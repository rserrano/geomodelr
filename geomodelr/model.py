
# Geomodelr query tool. Tools for using a geomodelr.com geological model.
# Copyright (C) 2016 Geomodelr, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import shared
import json
import datetime
import numpy as np
import sys
import cpp
from versions import upgrade_model

from numpy import linalg as la

# Validates that the feature is valid and can be used.
def validate_feature(f):
    # TODO
    pass

# Validate the transforms.
def validate_transform(transform, geology_type):
    
    # TODO
    if not 'type' in transform:
        raise shared.ModelException("transform needs a type")
    if transform['type'] not in ['topography', 'plane', 'identity']:
        raise shared.ModelException("transform types supported: topography, plane and identity")
    
    pass

# Validate a FeatureCollection.
def validate_feature_collection(fc):
    if not 'type' in fc:
        raise shared.ModelException("type is necessary for feature collection")
    if fc['type'] != 'FeatureCollection':
        raise shared.ModelException("type needs to be FeatureCollection")
    if not 'geology_type' in fc:
        raise shared.ModelException("geology type is necessary for feature collection")
    if not fc['geology_type'] in ['map', 'section', 'boreholes', 'dips', 'faults']:
        raise shared.ModelException("not recognized geology type of feature collection")
    if not 'transform' in fc:
        raise shared.ModelException("transform is necessary for feature collection")
    validate_transform(fc['transform'], fc['geology_type'])
    if not 'features' in fc:
        raise shared.ModelException("features are necessary in feature collection")
    for f in fc:
        validate_feature(f)

class GeologicalModel(cpp.Model):
    """
    Interface to query a Geological model from Geomodelr.com. The models in Geomodelr.com
    are saved in Geological JSON. A Geological JSON is a set of GeoJSON FeatureCollections 
    with a transformation. Go to Geomodelr.com, create a new model and use it with this 
    tool.
    """
    def __init__( self, geolojson, delete=True, params={'faults': 'basic'} ):
        """ 
        Initializes the geological model from a Geological JSON 
        file created in www.geomodelr.com.
        
        You can create a free user at https://geomodelr.com and it will allow you to create the 
        Geological Model. After you are finished, create a version and it will allow you
        to download it as a Geological JSON. You can use this constructor by loading the 
        json, like this::
        
            import json
            import geomodelr
            mfile = open('/path/to/your/version.json')
            geomodel = geomodelr.GeologicalModel(json.loads(mfile.read()))
        
        Args:
            (dict) geolojson: The Geological JSON.
        """
        
        self.geojson = geolojson
        
        upgrade_model( self.geojson )

        sections = []
        
        geomap = []
        topography = {}
        # Search map.
        for feature in self.geojson['features']:
            if feature['geology_type'] == 'map':
                gm = shared.points_index_repr(feature)
                geomap = [gm['points'], gm['polygons'], gm['units'], gm['lines'], gm['lnames']]
                topography = feature['transform']
                break
        else:
            self.geomap = None
        
        # First get the base section, which will locate all other sections.
        base_section = None
        for feature in self.geojson['features']:
            if feature['geology_type'] == 'section' and 'base' in feature['properties'] and feature['properties']['base']:
                base_section = feature
                break
        
        # Calculate direction and base point.
        try:
            if base_section['transform']['type'] == 'plane':
                if not 'orientation' in base_section['transform'] or base_section['transform']['orientation'] == 'vertical': 
                    base_line = base_section['transform']['line']
                    orientation = 'vertical'
                elif base_section['transform']['orientation'] == 'horizontal':
                    base_line = None
                    orientation = 'horizontal'
                else:
                    raise shared.ModelException("The base section should be a plane")
            else:
                raise shared.ModelException("The base section should be a plane")
        except TypeError:
            raise shared.ModelException("The model needs minimum a base cross section.")
        
        if orientation == 'vertical':
            base_point = np.array(base_line[0][:2])
            direction = np.array(base_line[1][:2])-base_point
            direction = direction/la.norm(direction)
        else:
            base_point = None
            direction = None
        abbox = [ float('inf'), float('inf'), float('inf'), -float('inf'), -float('inf'), -float('inf') ]
        for idx, feature in enumerate(self.geojson['features']):
            if feature['geology_type'] == 'section' and 'interpolation' in feature['properties'] and feature['properties']['interpolation']:
                if orientation == 'vertical':
                    cut, cs = shared.cross_idx_repr( feature, base_line )
                    
                    # To calculate aligned bounding box.
                    bds = shared.calc_line_bounds( feature['transform']['line'], base_line )
                    abbox[0] = min( bds[0], abbox[0] )
                    abbox[3] = max( bds[1], abbox[3] )
                    abbox[2] = min( cut, abbox[2] )
                    abbox[5] = max( cut, abbox[5] )

                    sect = [feature['name'], cut, cs['points'], cs['polygons'], cs['units'], cs['lines'], cs['lnames'], cs['anchored_lines']]
                    sections.append(sect)
                else:
                    cs = shared.points_index_repr(feature)
                    sect = [feature['name'], feature['transform']['height'], cs['points'], cs['polygons'], cs['units'], cs['lines'], cs['lnames'], cs['anchored_lines']]
                    sections.append(sect) 
        
        # Obtain the possible farthest cuts to add triangles towards them.
        bbox = self.geojson['bbox']
        if orientation == 'horizontal':
            abbox = bbox
        else:
            abbox[1] = bbox[2]
            abbox[4] = bbox[5]
        
        lines = self.geojson['properties']['lines']
        if orientation == 'horizontal':
            super(GeologicalModel, self).__init__(bbox, abbox, geomap, topography, sections, lines, params)
        else:
            super(GeologicalModel, self).__init__(bbox, abbox, list(base_point), list(direction), geomap, topography, sections, lines, params)
        
        self.make_matches()
        
        # Add units to model before deleting geojson.
        units = self.geojson['properties']['units'].keys()
        self.units = units
        
        # Save space.
        if delete:
            del self.geojson
        

    def make_matches(self):
        """ Prepares the model to query by matching polygons and lines.
            It finds which polygons, when projected to the next cross section,
            intersect. After that, it tries to match faults with the same name
            by triangulating them and trying to find a continuous set of triangles
            between the two lines that go from the ends to the other side.
        """
        super(GeologicalModel, self).make_matches()
    
    def print_information( self, verbose=False ):
        """
        Prints the information of the geological model just loaded.
        
        Prints the version, coordinate system and valid coordinates 
        that the geological model takes.
        
        Args:
            (boolean) verbose: You can print more information with verbose=True.
        """
        # Get name of the study.
        if 'name' in self.geojson:
            print "Geological Model Name:", self.geojson['name']
        
        else:
            print "No name"
        
        # Get the version of the study.
        if 'version' in self.geojson['properties']:
            version = self.geojson['properties']['version']
            timestamp = datetime.datetime.strptime(self.geojson['properties']['timestamp'], "%Y-%m-%dT%H:%M:%S.%fZ").isoformat()
            print "\tVersion %(version)s created at %(created)s" % {'version': version, 'created': timestamp} 
        
        # Bounding box.
        if 'crs' in self.geojson:
            crs = self.geojson['crs']
            print "Coordinate System:\n\t%s" % crs['properties']['name']
        
        
        # Bounding box.
        if 'bbox' in self.geojson:
            bbox = self.geojson['bbox']
            print "Bounding Box:\n\tMin. X: %s Min. Y: %s Min. Z: %s Max. X: %s Max. Y: %s Max. Z: %s" % tuple(bbox)
        
        units  = set()
        pprops = set()
        lnames = set()
        lprops = set()
        
        def verbose_info(collection):
            polys = filter(lambda f: f['geometry']['type'] == 'Polygon', collection['features'])
            for p in polys:
                if 'unit' in p['properties']:
                    units.add(p['properties']['unit'])
                pprops.update(set(p['properties'].keys()))
            
            flts = filter(lambda f: f['geometry']['type'] == 'LineString', collection['features'])
            for f in flts:
                if 'name' in f['properties']:
                    lnames.add(f['properties']['name'])
                lprops.update(set(f['properties'].keys()))
            
            return { 'polygons': len(polys), 'faults': len(flts) }
        
        maps = filter(lambda f: f['geology_type'] == 'map', self.geojson['features'])
        
        print "Content:"
        # Has map.
        if len(maps) == 0:
            print "\tGeological map not present"
        elif len(maps) == 1:
            print "\tGeological map present"
            if verbose:
                print "\t\tpolygons: %(polygons)s, faults: %(faults)s" % verbose_info(maps[0])
        else:
            print "\tMore than one map"
        
        # Has sections.
        sections = filter(lambda f: f['geology_type'] == 'section', self.geojson['features'])
        snames = map( lambda s: s['name'], sections )
        if len(sections) > 0:
            # Number of sections and names of sections.
            print "\tNumber of geological cross sections present: %s" % len(sections)
            print "\tNames of the cross sections: %s" % ", ".join(snames)
            if verbose:
                for s in sections:
                    print "\tSection %s" % s['name']
                    print "\t\tpolygons: %(polygons)s, faults: %(faults)s" % verbose_info(s)
        
        if verbose:
            # Polygons
            print "\tPolygon units present: %s" % ", ".join(units)
            print "\tPolygon properties present: %s" % ", ".join(pprops)
            # Fault names.
            print "\tFault names present: %s" % ", ".join(lnames)
            print "\tFault properties present: %s" % ", ".join(lprops)
    
    def validate( self ):
        """
        Validates that the Geological JSON has correct information.
        """
        # Validate basic GeoJSON contents.
        if type(self.geojson) != dict:
            raise shared.ModelException("You did not load a dictionary in geojson format")
        if not 'type' in self.geojson:
            raise shared.ModelException("No type in json")
        if self.geojson['type'] != 'GeologyCollection':
            raise shared.ModelException("type is not GeologyCollection")
        if not 'bbox' in self.geojson:
            raise shared.ModelException("bounding box (bbox) is required in model")
        if not 'crs' in self.geojson:
            raise shared.ModelException("coordinate system (crs) is required in model")
        if not 'type' in self.geojson['crs']:
            raise shared.ModelException("coordinate system (crs) type is required in model")
        if self.geojson['crs']['type'] != 'name':
            raise shared.ModelException("only name is supported for coordinate system (crs) type")
        if not 'properties' in self.geojson['crs']:
            raise shared.ModelException("properties in coordinate system (crs) is required in model")
        if not 'name' in self.geojson['crs']['properties']:
            raise shared.ModelException("name property in coordinate system (crs) is required in model")
        if self.geojson['crs']['properties']['name'].split(":")[0] != "EPSG":
            raise shared.ModelException("only EPSG codes are supported as coordinate systems")
        try:
            int(self.geojson['crs']['properties']['name'].split(":")[1])
        except ValueError:
            raise shared.ModelException("EPSG should be integer code")
        
        # Validate the feature collections of the GeoJSON.
        if not 'features' in self.geojson:
            raise shared.ModelException("features is necessary in model")
        for fc in self.geojson['features']:
            validate_feature_collection(fc)

def model_from_file(filename):
    """
    Entry point for the API. It creates the geological model 
    from the file path. The geological model is a model of 
    geomodelr.com, downloaded as a version.
    
    Args:
        (str) filename: The path to the Geological JSON file downloaded from 
        Geomodelr.com.

    Returns:
        (GeologicalModel): The output Geological model to query the geological
        units freely.
    """
    with open(filename) as f:
        m = GeologicalModel(json.load(f))
        return m

