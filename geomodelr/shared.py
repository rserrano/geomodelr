from shapely.geometry import Polygon, Point, LineString
import numpy as np
from numpy import linalg as la

def line_side(surface_line, cs_line):
    """
    Distance between two surface lines.
    """
    section_dir = np.array(surface_line[1])-np.array(surface_line[0])
    nsect = la.norm(section_dir)
    nvsect = section_dir/nsect
    v = np.array(cs_line[0])-np.array(surface_line[0])
    n = la.norm(v)
    if n < 1e-10:
        return 0
    else:
        nv = v/n
        s = nvsect[0]*nv[1] - nvsect[1] * nv[0]
        return s*n

def line_offset(surface_line, cs_line):
    """
    Offset to the first line.
    """
    section_dir = np.array(surface_line[1])-np.array(surface_line[0])
    nsect = la.norm(section_dir)
    nvsect = section_dir/nsect
    v  = np.array(cs_line[0])-np.array(surface_line[0])
    n  = la.norm(v)
    if n < 1e-10:
        return 0
    else:
        nv = v/n
        s = nvsect[0]*nv[0] + nvsect[1] * nv[1]
        return s*n

def nodes_offset(nodes, offset):
    """
    Returns the set of nodes but transformed so that the begining of coordinates is 
    The first point of the base, with get_base_transform
    """
    for node in nodes:
        node[0] = node[0] + offset

class Cross:
    """ Cross is a representation for the cross section where points-polygons are separated
        and also shapely_polygons are calculated. """
    def __init__(self, cut, points, polygons, units, lines, lnames):
        self.cut = cut
        self.points = points
        # Import polygons with shapely versions.
        self.polygons = []
        self.shapely_polygons = []
        self.units = []
        for idx, poly in enumerate(polygons):
            ppoly = []
            for ring in poly:
                pring = [self.points[n] for n in ring]
                ppoly.append(pring)
            try:
                self.shapely_polygons.append(Polygon(ppoly[0], ppoly[1:]))
                self.polygons.append(poly)
                self.units.append(units[idx])
            except:
                # Ignore errors.
                pass
        reorder = sorted(range(len(self.polygons)), key=lambda o: abs(self.shapely_polygons[o].area))
        self.polygons = [self.polygons[i] for i in reorder]
        self.shapely_polygons = [self.shapely_polygons[i] for i in reorder]
        self.units = [self.units[i] for i in reorder]
        
        
        self.lines = []
        self.shapely_lines = []
        self.lnames = []
        
        # Import faults, with shapely versions.
        for idx, line in enumerate(lines):
            pline = [self.points[n] for n in ring]
            try:
                
                self.shapely_lines.append(LineString(pline))
                self.lines.append(line)
                self.lnames.append(lnames[idx])
            except Exception as e:
                # Ignore errors.
                print(e)
                pass

def points_index_repr(geojson):
    """ Returns the geojson data a set of points with indexed lines and polygons.  """
    lines = []
    lnames = []
    polygons = []
    units = []
    pointsd = {}
    
    def ipoint(point):
        """
        Gets the index of the point and increases the counter if it does not exist.
        """
        tpoint = tuple(point)
        if tpoint in pointsd:
            ipoint = pointsd[tpoint]
        else:
            ipoint = len(pointsd)
            pointsd[tpoint] = ipoint
        return ipoint
    
    for feature in geojson['features']:
        geom = feature['geometry']
        if geom['type'] == 'Polygon':
            ipolygon = []
            for ring in geom['coordinates']:
                iring = []
                for point in ring:
                    iring.append(ipoint(point))
                ipolygon.append(iring)
            polygons.append(ipolygon)
            if 'unit' in feature['properties']:
                units.append(feature['properties']['unit'])
            else:
                units.append("NONE")
        elif feature['geometry']['type'] == 'LineString':
            iline = []
            for point in geom['coordinates']:
                iline.append(ipoint(point))
            lines.append(iline)
            if 'name' in feature['properties']:
                lnames.append(feature['properties']['name'])
            else:
                lnames.append("")
    points = [ None for i in xrange(len(pointsd)) ]
    for p in pointsd:
        points[pointsd[p]] = list(p)
    
    return { 'points': points, 'lines': lines, 'polygons': polygons, 'units': units, 'lnames': lnames }

def cross_from_geojson(geojson, base_line):
    """
    Converts a GeoJSON cross section into a 
    """
    pi = points_index_repr(geojson)
    line = geojson['transform']['line']
    offset = line_offset(base_line, line)
    nodes_offset(pi['points'], offset)
    cut = line_side(base_line, line)
    return Cross(cut, pi['points'], pi['polygons'], pi['units'], pi['lines'], pi['lnames'])


