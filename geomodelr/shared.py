
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

from shapely.geometry import Polygon, Point, LineString
import numpy as np
from numpy import linalg as la
from scipy.spatial import Delaunay
import itertools
import math

class GeometryException(Exception):
    pass

class ModelException(Exception):
    pass

# Method converts a set of shapes, (MultiGeometry, MultiPolygon, etc.), into a list of simpler shapes.
def shape_list(shape, sh_type):
   try:
       shapes = []
       for sh in shape.geoms:
           if sh_type is None or sh.geometryType() == sh_type:
               shapes.append(sh)
       return shapes
   except:
       if sh_type is None or shape.geometryType() == sh_type:
           return [shape]
       return []

def over_point( points ):
    """
    Returns a point perpendicular to the plane and above all the others. 
    This is to avoid errors with the triangulation.
    """
    mxx = -float('inf')
    mxy = -float('inf')
    mxz = -float('inf')
    mnx =  float('inf')
    mny =  float('inf')
    mnz =  float('inf')
    
    for p in points:
        mxx = max(p[0], mxx)
        mxy = max(p[1], mxy)
        mxz = max(p[2], mxz)
        mnx = min(p[0], mnx)
        mny = min(p[1], mny)
        mnz = min(p[2], mnz)
    
    maxd = -1
    maxp = None
    for p in points:
        dx = max(mxx-p[0], p[0]-mnx)
        dy = max(mxy-p[0], p[0]-mny)
        dz = max(mxz-p[0], p[0]-mnz)
        sq = math.sqrt(dx*dx+dy*dy+dz*dz)
        if sq > maxd:
            maxd = sq
            maxp = p
    
    p0 = np.array(maxp)
    A = np.array( map( lambda p: [p[0], p[1], 1], points ) )
    y = np.array( map( lambda p: p[2], points ) )
    x = la.lstsq(A, y)[0]
    d = math.sqrt( 1 + x[0]*x[0] + x[1]*x[1] )
    x /= d
    n = np.array([-x[0], -x[1], 1.0/d])
    sqd = math.sqrt((mxx-mnx)*(mxx-mnx)+(mxy-mny)*(mxy-mny)+(mxz-mnz)*(mxz-mnz))
    
    return p0 + n*sqd
    
# Triangulates a set of points and returns the triangles filtered by the given filter.
def triangulate(points, fltr=lambda x: True):

    tetras = Delaunay(points)
    tris = set()
    for tet in tetras.simplices.copy():
        for tri in itertools.combinations(tet,3):
            srt = tuple(sorted(tri))
            if fltr(srt):
                tris.add(srt)
    return list(tris)

# Triangulates a set of points filtering them because they are on both sides of the line.
def opposite_triangles(points, na):
    n = len(points)
    pt = over_point(points)
    def is_opposite(tri):
        if n in tri:
            return False
        if not tri[0] < na or not tri[2] >= na:
            return False
        if tri[0]+1 == tri[1] and tri[1] < na:
            return True
        if tri[1]+1 == tri[2] and tri[1] >= na:
            return True
        return False

    return map( lambda tri: map(int, tri), triangulate(points+[pt], is_opposite) )

# Calculates the distance between two parallel lines.
def line_side(surface_line, cs_line):
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

# Calculates how much distance this section has to be offset to be aligned with the first.
def line_offset(surface_line, cs_line):
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

# Returns the set of nodes but transformed so that the begining of coordinates is 
# The first point of the base, with get_base_transform
def nodes_offset(nodes, offset):
    for node in nodes:
        node[0] = node[0] + offset

def fault_points_index_repr(fault):
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
    itris = []
    for triangle in fault:
        itris.append( map( ipoint, triangle ) )
    
    points = [ None for i in xrange(len(pointsd)) ]
    for p in pointsd:
        points[pointsd[p]] = list(p)
    return {'triangles': itris, 'points': points}    

def faults_points_index_repr(faults):
    return { k: fault_points_index_repr(v) for k, v in faults.iteritems() }

# Returns from the geojson data a set of points with indexed lines and polygons.  
def points_index_repr(geojson):
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

# Returns the points index representation of a cross section.
def cross_idx_repr(geojson, base_line):
    pi = points_index_repr(geojson)
    line = geojson['transform']['line']
    offset = line_offset(base_line, line)
    nodes_offset(pi['points'], offset)
    cut = line_side(base_line, line)
    return ( cut, pi )

# returns an encoded string from a unicode, forcing it
# reaaaally forcing it.
# This is probably the ugliest function ever.
# Probably someone will shame me for this. :(.
def force_encode(f):
    try:
        return f.encode()
    except:
        pass
    try:
        return f.encode('utf')
    except:
        pass
    try:
        return f.encode('utf8')
    except:
        pass
    try:
        return f.encode('utf16')
    except:
        pass
    return f.encode('utf32')

