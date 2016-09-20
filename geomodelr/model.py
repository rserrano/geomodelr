import matching
import shared
import faults
import query
import json
import datetime
import numpy as np
import sys


from numpy import linalg as la

def validate_feature(f):
    """
    Validates that the feature is valid and can be used.
    """
    # TODO
    pass

def validate_transform(transform, geology_type):
    """
    Validate the transforms.
    """
    # TODO
    if not 'type' in transform:
        raise shared.ModelException("transform needs a type")
    if transform['type'] not in ['topography', 'plane', 'identity']:
        raise shared.ModelException("transform types supported: topography, plane and identity")
    
    pass

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

class GeologicalModel(object):
    def __init__( self, geojson ):
        """ 
        Initializes the geological model from a geojson file.
        """
        self.geojson = geojson
    
    def validate( self ):
        """
        Validates that the geojson is correct. 
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
    
    def print_information( self, verbose=False ):
        """
        Validates the information of the geological model.
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
    
    def calc_cache( self ):
        """
        Calculates objects which are easier to work with and faster
        but can't be serialized, like shapely polygons.
        """
        
        self.sections = []
        self.section_idxs = []
        
        # Search map.
        for feature in self.geojson['features']:
            if feature['geology_type'] == 'map':
                self.topography = feature['transform']
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
        base_line = base_section['transform']['line']
        self.base_point = np.array(base_line[0][:2])
        self.direction = np.array(base_line[1][:2])-self.base_point
        self.direction = self.direction/la.norm(self.direction)
        
        # Then add non serializable cross sections, (TODO: Kill this some day.)
        for idx, feature in enumerate(self.geojson['features']):
            if feature['geology_type'] == 'section' and 'interpolation' in feature['properties'] and feature['properties']['interpolation']:
                self.sections.append(shared.cross_from_geojson(feature, base_line))
                self.section_idxs.append(idx)
        
        # Then reorder all the sections.
        reorder = sorted(range(len(self.sections)), key=lambda o: self.sections[o].cut)
        self.sections = [self.sections[i] for i in reorder]
        self.section_idxs = [self.section_idxs[i] for i in reorder]
        
        # Check that interpolation exists, if so, add a shortcut to the matching.
        if 'interpolation' in self.geojson:
            # TODO: Check that it's in the correct order or has not changed.
            self.matching = self.geojson['interpolation']['matching']
        
        # Check that faults exist. If so, make a side to side representation.
        fexist = False
        for idx, feature in enumerate(self.geojson['features']):
            # TODO: Check that it's in the correct order or has not changed.
            if feature['geology_type'] == 'faults':
                if not fexist:
                    fexist = True
                    self.faults = []
                    self.joined_faults = []
                self.faults.append(faults.align_fault_with(feature, self.model_point))
                for surface in feature['features']:
                    name = surface['properties']['name']
                    fplane = surface['geometry']['coordinates']
                    if name in self.joined_faults:
                        self.joined_faults[name] += fplane
                    else:
                        self.joined_faults[name] = fplane
    
    def has_interpolation( self ):
        return 'interpolation' in self.geojson

    def calc_interpolation( self ):
        """
        Adds a property to self.geojson called interpolation 
        with a matching that defines which polygons match each
        other.
        """
        self.geojson['interpolation'] = { 'type': 'parallel', 'matching': [] }
        mtch = self.geojson['interpolation']['matching']
        for idx in xrange(len(self.sections)-1):
            mtch.append(matching.create_match(self.sections[idx], self.sections[idx+1]))
            name1 = self.geojson['features'][self.section_idxs[idx]]['name']
            name2 = self.geojson['features'][self.section_idxs[idx+1]]['name']
            mtch[-1]['name'] = [name1, name2]
        self.matching = self.geojson['interpolation']['matching']
    
    def has_faults( self ):
        for feature in self.geojson['features']:
            if feature['geology_type'] == 'faults':
                return True
        return False
    
    def calc_faults( self ):
        """ 
        Adds a FeatureCollection to the GeoJSON with 
        surfaces containing triangles that
        represent a fault between two cross sections.
        """
        self.faults = []
        self.joined_faults = { }
        for idx in xrange(len(self.section_idxs)-1):
            mdx = self.section_idxs[idx]
            ndx = self.section_idxs[idx+1]
            fplanes = faults.faultplanes_between_sections(self.geojson['features'][mdx], self.geojson['features'][ndx])
            
            for name, fplane in fplanes.iteritems():
                if name in self.joined_faults:
                    self.joined_faults[name] += fplane
                else:
                    self.joined_faults[name] = fplane
            
            fplanes = faults.faultplanes_to_featurecollection(fplanes)
            name1 = self.geojson['features'][mdx]['name']
            name2 = self.geojson['features'][ndx]['name']
            fplanes['name'] = [name1, name2]
            self.geojson['features'].append(fplanes)
            self.faults.append(faults.align_fault_with(fplanes, self.model_point))

        
    def model_point(self, point):
        """
        Returns the point in geological model representation.
        parameters:
            point: point to query in spatial coordinates
        return:
            point: point to query in model coordinates.
        """
        ppoint = np.array([point[0], point[1]])
        v = ppoint-self.base_point
        n = la.norm(v)
        if n < 1e-10:
            return (0, point[2], 0)
        else:
            nv = v/n
            c = self.direction[0]*nv[1] - self.direction[1] * nv[0]
            d = self.direction[0]*nv[0] + self.direction[1] * nv[1]
            
            return [d*n, point[2], c*n]
    
    def inverse_point(self, point):
        """
        Given a point in model coordinates, returns its inverse.
        parameters:
            point: point to query in model coordinates
        return:
            point: point in normal coordinates.
        """
        dsq = self.direction[0]**2 + self.direction[1]**2
        v1 = (point[2]*self.direction[0] + point[0]*self.direction[1])/dsq
        if math.fabs(self.direction[1]) < 1e-9:
            v0 = point[0]/self.direction[0]
        else:
            v0 = (self.direction[0]*v1-point[2])/self.direction[1]
        return [ self.base_point[0] + v0, self.base_point[1] + v1, point[1] ]
    
    def query_point(self, point):
        """
        Returns the formation in which the point is located.
        parameters:
            point:  the 3 coordinates of the point to query.
        """
        point = self.model_point(point)
        return self.query_point_cut((point[0], point[1]), point[2])
    
    def query_point_topo(self, point):
        """
        Convenience function that also evaluates if the point is above
        the topography and returns none.
        parameters:
            point: the 3 coordinates of the point to query.
        """
        height = self.query_height((point[0], point[1]))
        if height < point[2] or height is None:
            return None
        else:
            return self.query_point(point)
        
    def query_point_cut(self, point, cut):
        
        if cut <= self.sections[0].cut:
            return query.query_cross_point(self.sections[0], point)
        elif cut >= self.sections[-1].cut:
            return query.query_cross_point(self.sections[-1], point)
        
        rng = [0, len(self.sections)-1]
        while rng[0] < (rng[1] - 1):
            mid = (rng[0] + rng[1])/2
            if cut > self.sections[mid].cut:
                rng[0] = mid
            elif cut < self.sections[mid].cut:
                rng[1] = mid
            else:
                rng = [mid, mid+1]
        
        return query.query_single_point(self.sections[rng[0]], self.sections[rng[1]], self.matching[rng[0]], self.faults[rng[0]], point, cut)
    
    def query_height(self, point):
        """
        Queries the height at a point in the plane.
        Returns the closest point if outside the bounds.
        """
        abspos = [(point[0]-self.topography['point'][0])/self.topography['sample'][0], 
                  (point[1]-self.topography['point'][1])/self.topography['sample'][1]]
        
        x = int(abspos[0])
        y = int(abspos[1])
        
        if x < 0:
            x = 0
        
        if y < 0:
            y = 0
        
        if x >= self.topography['dims'][0]:
            x = self.topography['dims'][0]-1
        
        if y >= self.topography['dims'][1]:
            y = self.topography['dims'][1]-1
        
        return self.topography['heights'][x][y]
    
    def intersect_plane( self, plane ):
        return faults.find_faults_plane_intersection( self.joined_faults, plane )

    def intersect_planes( self, planes ):
        return faults.find_faults_multiple_planes_intersection( self.joined_faults, planes )

