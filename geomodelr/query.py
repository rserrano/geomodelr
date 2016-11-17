from shapely.geometry import Polygon, Point, LineString
from shapely import speedups
# if speedups.available:
#     speedups.enable()
import numpy as np
from numpy import linalg as la
import itertools
import math

# Functions to query points.
def crosses_triangles( faults, point, cut ):
    pt2 = Point(point)
    pt3 = np.array((point[0], point[1], cut))
    for name, triangles in faults.iteritems():
        for tr in triangles:
            if tr.triangle.contains( pt2 ):
                d = np.dot( pt3-tr.point, tr.normal )
                if d < 0:
                    return -1
                else:
                    return 1
    return 0
        
def query_cross_point( cross, point ):
    """ 
    Queries a point in a single cross. Useful for querying outside of the range of the crosss. 
    """
    min_dist = float('inf')
    min_match = None
    sh_point = Point(point[0], point[1])
    # Searches for the closest polygon at that point without a 'NONE' geological unit.
    for idx in xrange(len(cross.polygons)):
        if cross.units[idx] == 'NONE':
            continue
        if cross.shapely_polygons[idx].contains(sh_point):
            return cross.units[idx]
        dist = cross.shapely_polygons[idx].distance(sh_point)
        if dist < min_dist:
            min_dist = dist
            min_match = idx
    # Returns NONE if none was found, or the closest match.
    return 'NONE' if min_match is None else cross.units[min_match]

def query_single_point(cross_a, cross_b, match, faults, point, cut):
    """ 
    Queries a point between two crosss. 
    """
    # Check if it crosses a fault plane. If so, just query the cross section at the side of the cross.
    crosses = crosses_triangles( faults, point, cut )
    
    if crosses != 0:
        if crosses < 0:
            return query_cross_point(cross_a, point)
        else:
            return query_cross_point(cross_b, point)

    dist = (cross_b.cut-cross_a.cut)
    mult_a = (cut-cross_a.cut)/dist
    mult_b = (cross_b.cut-cut)/dist
    min_dist  = float('inf')
    min_match = None
    for idx in xrange(len(match['match'])):
        vect_a = (mult_a*match['rays'][idx][0], mult_a*match['rays'][idx][1])
        vect_b = (mult_b*match['rays'][idx][0], mult_b*match['rays'][idx][1])
        point_a = Point(point[0]-vect_a[0], point[1]-vect_a[1])
        point_b = Point(point[0]+vect_b[0], point[1]+vect_b[1])
        if cross_a.shapely_polygons[match['match'][idx][0]].contains(point_a) and \
           cross_b.shapely_polygons[match['match'][idx][1]].contains(point_b):
            return cross_a.units[match['match'][idx][0]]
        d_a = cross_a.shapely_polygons[match['match'][idx][0]].distance(point_a)
        d_b = cross_b.shapely_polygons[match['match'][idx][1]].distance(point_b)
        
        dist = (((1-abs(mult_a))*d_a + (1-abs(mult_b))*d_b))
        if dist < min_dist:
            min_match = idx
            min_dist = dist
    
    # If no matches, (when everything is NONE), just return the match at the closest cross section.
    if min_match is None:
        if mult_a < mult_b:
            return query_cross_point(cross_a, point)
        else:
            return query_cross_point(cross_b, point)
    
    # Return the given match.
    return cross_a.units[match["match"][min_match][0]]

