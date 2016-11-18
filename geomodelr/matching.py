"""
	Geomodelr query tool. Tools for using a geomodelr.com geological model.
	Copyright (C) 2016 Geomodelr, Inc.
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Affero General Public License for more details.
	
	You should have received a copy of the GNU Affero General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

from scipy.spatial import Delaunay
from scipy.spatial.distance import euclidean as dist
import numpy as np
from numpy import linalg as la
import itertools

def classify_triangle(tri, a_len):
    """
    Classifies classifies the triangle for the number of sides it has in each cross section.
    parameters:
        tri:    the triangle.
        a_len:  the separation between the points in each cross section.
    """
    side_a = []
    side_b = []
    for node in tri:
        if ( node < a_len ):
            side_a.append(node)
        else:
            side_b.append(node)
    n = len(side_a)
    if ( n == 0 ):
        return ("only_b", tuple(side_b))
    if ( n == 1 ):
        return ("b", tuple(side_b), side_a[0])
    if ( n == 2 ):
        return ("a", tuple(side_a), side_b[0])
    if ( n == 3 ):
        return ("only_a", tuple(side_a))

def polygons_dual(cross):
    """ 
    Creates a graph that represents the connection between points and polygons, with next 
    parameters:
        cross:      the cross section with the polygons and points.
    return:
        belongs_to: the set of (polygon, ring) that a point belongs to.
        next_node:  is the next node in the ring after this.
    """
    belongs_to = [[] for i in range(len(cross.points))]
    next_node  = [[] for i in range(len(cross.points))]
    
    for npol, poly in enumerate(cross.polygons):
        for nring, ring in enumerate(poly):
            n = len(ring)
            for nnod, node in enumerate(ring):
                belongs_to[node].append((npol, nring))
                next_node[node].append(ring[(nnod+1)%n])
    return {"belongs_to": belongs_to, "next_node": next_node}


def concat_duals(dual_a, dual_b, a_len):
    """ 
    Concatenates the duals created in polygons_dual re enumerating
    the points in the second cross section.
    parameters:
        dual_a:     The output for polygons_dual for cross_a
        dual_b:     The output for polygons_dual for cross_b
        a_len:      The number of points in cross_a, used to renumber.
    return:
        belongs_to: The (polygon, ring) that the node belongs to.
        next_node:  The next node in the polygon ring, with absolute indices.
    """
    belongs_to = dual_a["belongs_to"] + dual_b["belongs_to"]
    next_node = dual_a["next_node"] + dual_b["next_node"]
    
    for idx in xrange(a_len, len(next_node)):
        for jdx in xrange(len(next_node[idx])):
            next_node[idx][jdx] += a_len
    
    return {"belongs_to": belongs_to, "next_node": next_node}

def segment_in_polygons(dual, seg):
    """ 
    Returns all the pairs (polygon, ring) that are shared by the segment 
    """
    pols = set()
    for num, nnod in enumerate(dual["next_node"][seg[0]]):
        if ( nnod == seg[1] ):
            pols.add(dual["belongs_to"][seg[0]][num])
    for num, nnod in enumerate(dual["next_node"][seg[1]]):
        if ( nnod == seg[0] ):
            pols.add(dual["belongs_to"][seg[1]][num])
    return list(pols)

def connected_length(cross_a, cross_b):
    """ One of the most important parts of the algorithm.
        Its result is the segments that are connected between 
        two polygons of the same unit by the delaunay 
        triangulation. """
    # Create matrix of connected lengths between polygons.
    con_len = np.zeros((len(cross_a.polygons), len(cross_b.polygons)))
    a_len = len(cross_a.points)
    
    # Set of points from the cs a has coordinate z -1.
    points_a = np.zeros((len(cross_a.points), 3))
    points_a[:,:-1] = np.array(cross_a.points)
    points_a[:,2] = -1
    
    # Set of points from the cs b has coordinate z 1.
    points_b = np.zeros((len(cross_b.points), 3))
    points_b[:,:-1] = np.array(cross_b.points)
    points_b[:,2] = 1
    
    # Triangulate the points.
    points = np.concatenate((points_a, points_b), axis=0)
    tetras = Delaunay(points)
    
    # 
    dual_a = polygons_dual(cross_a)
    dual_b = polygons_dual(cross_b)
    dual = concat_duals(dual_a, dual_b, a_len)
    
    # From the tetrahedra get all the unique triangles. 
    tris = set()
    for tet in tetras.simplices.copy():
        for tri in itertools.combinations(tet,3):
            tris.add(tuple(sorted(tri)))
    
    # Get the segments that will be added to the connected length
    to_sum = []
    for tri in tris:
        clsd = classify_triangle(tri, a_len)
        if ( len(clsd[0]) == 1): # if the triangle goes from one side to the other.
            seg = segment_in_polygons(dual, clsd[1])
            for pol1 in seg:
                for pol2 in dual["belongs_to"][clsd[2]]:
                    if clsd[0] == "a":
                        unit_a = cross_a.units[pol1[0]]
                        unit_b = cross_b.units[pol2[0]]
                        if unit_a == unit_b:
                            tup = (clsd[0], clsd[1], (pol1[0], pol2[0]))
                            to_sum.append(tup)
                    else: 
                        # clsd[0] == "b"
                        unit_b = cross_b.units[pol1[0]]
                        unit_a = cross_a.units[pol2[0]]
                        if unit_a == unit_b:
                            tup = (clsd[0], clsd[1], (pol2[0], pol1[0]))
                            to_sum.append(tup)
    
    for tup in to_sum:
        con_len[tup[2][0],tup[2][1]] += la.norm(points[tup[1][0],:-1]-points[tup[1][1],:-1])
    
    return con_len

def match_polygons(cross_a, cross_b):
    """ 
    Returns the matches that can be made by the polygons using the connected length.
    """
    con_len = connected_length(cross_a, cross_b)
    (m, n) = con_len.shape
    sel_a = [False for i in range(m)]
    sel_b = [False for i in range(n)]
    pairs = []
    for i in xrange(m):
        for j in xrange(n):
            if cross_a.units[i] != 'NONE' and \
               cross_a.units[i] == cross_b.units[j] and \
               con_len[i,j] > 0.0:
                pairs.append((i, j))
    
    sorted_pairs = sorted(pairs, key=lambda x: con_len[x[0], x[1]], reverse=True)
    
    matched_pairs = []
    for p in sorted_pairs:
        if ( sel_a[p[0]] and sel_b[p[1]] ):
            continue
        matched_pairs.append(p)
        sel_a[p[0]] = True
        sel_b[p[1]] = True
    return matched_pairs

def create_match(cross_a, cross_b):
    """ 
    Returns the matches, with the rays that point into the direction between polygons
    And it's ordered by the sum of areas of the matched polygons.
    parameters:
        cross_a:    the first cross section.
        cross_b:    the second cross section.
    """
    match = match_polygons(cross_a, cross_b)
    
    rays = []
    
    for c in match:
        ctr_a = cross_a.shapely_polygons[c[0]].centroid
        ctr_b = cross_b.shapely_polygons[c[1]].centroid
        # ray   = (ctr_b.x-ctr_a.x, ctr_b.y-ctr_a.y)
        ray = (0.0, 0.0)
        rays.append(ray)
    
    areasum = lambda o: abs(cross_a.shapely_polygons[match[o][0]].area)+abs(cross_b.shapely_polygons[match[o][1]].area)
    reorder = sorted(range(len(match)), key=areasum)
    reorder_match = [match[i] for i in reorder]
    reorder_rays  = [rays[i]  for i in reorder]
    areas = [ areasum(i) for i in reorder]
    
    return { 'match': reorder_match, 'rays': reorder_rays, 'areas': areas }
