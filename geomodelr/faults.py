import shared
from scipy.spatial.distance import euclidean as dist
import math
import numpy as np
from numpy import linalg as la
from scipy.spatial import Delaunay
import sortedcontainers

def angles(t):
    """
    Given a triangle it finds the angles of the triangle.
    """
    st = map(lambda a: a**2, t)
    cos = [(st[1] + st[2] - st[0])/(2*t[1]*t[2]), (st[0] + st[1] - st[2])/(2*t[0]*t[1])]
    ang = map(math.acos, cos)
    ang.append(math.pi - (ang[0] + ang[1]))
    return ang

def angle_between(p1, p2, p3, p4):
    """
    Finds the angle between two triangles with 
    shared edge (p1, p2) in 3D.
    """
    v1 = np.array(p2)-np.array(p1)
    v2 = np.array(p3)-np.array(p1)
    v3 = np.array(p4)-np.array(p1)
    c1 = np.cross(v2,v1)
    nc1 = c1/la.norm(c1)
    c2 = np.cross(v3,v1)
    nc2 = c2/la.norm(c2)
    dt = np.dot(nc1, nc2)
    if dt > 1.0:
        dt = 1.0
    if dt < -1.0:
        dt = -1.0
    ang = math.acos(dt)
    return ang

def next_and_check(tri,edg,na):
    """
    Gets the next edge and also the edge that now has been finished.
    """
    if tri[1] < na:
        # Case tri[0], tri[1] are in the same edge tri[2] the contrary.
        # (tri[0], tri[1]) need to be removed.
        if edg[0] == tri[0]:
            # Case edg[0] is tri[0], then next should be (tri[1], tri[2])
            return ( (tri[1], tri[2]), (tri[0], tri[1]) )
        else:
            # Case edg[0] is tri[1], then next should be (tri[0], tri[2])
            return ( (tri[0], tri[2]), (tri[0], tri[1]) )
    else:
        # Case tri[1], tri[2] are in the same edge, tri[0] the contrary.
        # (tri[1], tri[2]) need to be removed.
        if edg[1] == tri[1]:
            # Case edg[0] is tri[0], then next should be (tri[1], tri[2])
            return ( (tri[0], tri[2]), (tri[1], tri[2]) )
        else:
            # Case edg[0] is tri[1], then next should be (tri[0], tri[2])
            return ( (tri[0], tri[1]), (tri[1], tri[2]) )
 
def fault_plane_for_lines(la, lb):
    """
    Get the faults plane between lines la, lb.
    """
    na = len(la)
    nb = len(lb)

    def is_between_faults(tri):
        """
        Filter the triangles so that only triangles in the edges 
        and that are tied between faults are accepted.
        It's required that tri is ordered.
        """
        if not tri[0] < na or not tri[2] >= na:
            return False
        if tri[0]+1 == tri[1] and tri[1] < na:
            return True
        if tri[1]+1 == tri[2] and tri[1] >= na:
            return True
        return False
    
    # Obtain the triangles that fulfil oposite
    tris = shared.triangulate(la+lb, is_between_faults)
    
    
    # pt is the point depending on the index.
    pt = lambda i: lb[i-na] if i >= na else la[i]
    pt_dist=lambda e: dist(pt(e[0]), pt(e[1]))
    pos_start = sortedcontainers.SortedListWithKey(key=pt_dist)
    
    # Obtain the edges count in point distance order.
    for tri in tris:
        edgs = []
        edgs.append((tri[0], tri[2]))
        if tri[1] < na:
            edgs.append((tri[1], tri[2]))
        else:
            edgs.append((tri[0], tri[1]))
        for edg in edgs:
            if (edg[0] == 0 or edg[0] == na-1) and (edg[1] == na or edg[1] == na+nb-1):
                pos_start.add(edg)
    
    if not len(pos_start):
        raise shared.GeometryException("Could not find a configuration to start the triangulation.")
    
    # Create a set of triangles and order by weight.
    def weight( tri ):
        """
        Minimum angle in triangle denotes the weight.
        """
        dt = [dist(pt(tri[0]), pt(tri[1])), 
              dist(pt(tri[1]), pt(tri[2])), 
              dist(pt(tri[2]), pt(tri[0]))]
        ang = angles(dt)
        return min(ang)
    
    sorted_triangles = sortedcontainers.SortedListWithKey( tris, weight )
    
    # Obtain the minimum count for the edges, (hopefully is one, but, well, weird cases).
    mincnt = min( map(pos_start.count, pos_start) )
    
    # Get the starting edge 
    for edg in pos_start:
        if pos_start.count(edg) == mincnt:
            start = edg
            break       
    nxtedg = start
    
    output_triangles = []
    
    has_edge = lambda tri, edg: (edg[0] in tri) and (edg[1] in tri)
    while nxtedg is not None:
        for tri in sorted_triangles:
            if has_edge(tri, nxtedg):
                sorted_triangles.remove(tri)
                output_triangles.append(tri)
                edgs = next_and_check(tri,nxtedg,na)
                burnedg = edgs[1]
                toremove = filter( lambda trr: has_edge(trr, burnedg) or has_edge(trr, nxtedg), sorted_triangles )
                for trr in toremove:
                    sorted_triangles.remove(trr)
                nxtedg = edgs[0]
                break
        else:
            nxtedg = None
    
    if len(sorted_triangles) != 0:
        raise shared.GeometryException("The faults were not completely triangulated")
    
    return output_triangles

def fault_planes_between_sections(cross_a, cross_b):
    """
    Creates the fault planes given the cross sections with faults with the same name.
    """
    # Create map of related faults.
    rel_faults = {  }
    for idx, nam in enumerate(cross_a.lnames):
        # Avoid two with the same name.
        if not nam in rel_faults:
            rel_faults[nam] = [idx]
    
    for idx, nam in enumerate(cross_b.lnames):
        # Check if present and avoid two of the same name.
        if nam in rel_faults and len(rel_faults[nam]) == 1:
            rel_faults[nam].append(idx)
        
    val_faults = filter( lambda f: len(rel_faults[f]) == 2, rel_faults )
    fault_planes = []
    for nam in val_faults:
        # Find which of the four points relate to each other.
        # First find the triangulation of the points at both sides.
        la = map( lambda n: cross_a.points[n]+[cross_a.cut], cross_a.lines[idx_a])
        lb = map( lambda n: cross_b.points[n]+[cross_b.cut], cross_b.lines[idx_b])
        
        fplane = fault_plane_for_lines(la, lb)
        na = len(la)
        idx_fplane = []
        for tri in fplane:
            idx_tri = []
            for n in tri:
                if n < na:
                    idx_tri.append((0, cross_a.lines[idx_a][n]))
                else:
                    idx_tri.append((1, cross_b.lines[idx_b][n-na]))
            idx_fplane.append(tri)
        fault_planes.append( [nam, fplane] )
    return fault_planes

