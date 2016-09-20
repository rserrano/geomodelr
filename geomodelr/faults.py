import shared
from scipy.spatial.distance import euclidean as dist
import math
import numpy as np
from numpy import linalg as la
from scipy.spatial import Delaunay
import sortedcontainers
from shapely.geometry import Polygon, Point, LineString
from copy import deepcopy

def angles(t):
    """
    Given a triangle it finds the angles of the triangle.
    """
    st = map(lambda a: a**2, t)
    cos = [(st[1] + st[2] - st[0])/(2*t[1]*t[2]), (st[0] + st[1] - st[2])/(2*t[0]*t[1])]
    ang = map(math.acos, cos)
    ang.append(math.pi - (ang[0] + ang[1]))
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
 
def faultplane_for_lines(l_a, l_b):
    """
    Get the faults plane between lines la, lb.
    """
    na = len(l_a)
    nb = len(l_b)

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
    tris = shared.triangulate(l_a+l_b, is_between_faults)
    
    
    # pt is the point depending on the index.
    pt = lambda i: l_b[i-na] if i >= na else l_a[i]
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
    
    # # Obtain the minimum count for the edges, (hopefully is one, but, well, weird cases).
    # mincnt = min( map(pos_start.count, pos_start) )
    # 
    # # Get the starting edge 
    # for edg in pos_start:
    #     if pos_start.count(edg) == mincnt:
    #         start = edg
    #         break       
    start = pos_start[0]

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

def faultplanes_between_sections(cross_a, cross_b):
    """
    Creates the fault planes given the cross sections with faults with the same name.
    """
    # Create map of related faults.
    rel_faults = {  }
    for idx, feature in enumerate(cross_a['features']):
        if 'name' in feature['properties']:
            nam = feature['properties']['name']
            # Avoid two with the same name.
            if not nam in rel_faults:
                rel_faults[nam] = [idx]
    
    for idx, feature in enumerate(cross_b['features']):
        if 'name' in feature['properties']:
            nam = feature['properties']['name']
            # Avoid two with the same name.
            if nam in rel_faults and len(rel_faults[nam]) == 1:
                rel_faults[nam].append(idx)
    
    pa = np.array(cross_a['transform']['line'][0]) 
    va = np.array(cross_a['transform']['line'][1]) - pa
    va = va/la.norm(va)
    def transform_a(point):
        """
        Only supports transversal cs for now.
        """
        return ( pa[0] + va[0]*point[0], pa[1] + va[1]*point[0], point[1] )
    
    pb = np.array(cross_b['transform']['line'][0]) 
    vb = np.array(cross_b['transform']['line'][1]) - pb
    vb = vb/la.norm(vb)
    def transform_b(point):
        """
        Only supports transversal cs for now.
        """
        return ( pb[0] + vb[0]*point[0], pb[1] + vb[1]*point[0], point[1] )
    
    faultplanes = {}
    for nam, idx in rel_faults.iteritems():
        # Find which of the four points relate to each other.
        # First find the triangulation of the points at both sides.
        if len(idx) != 2:
            continue
        
        l_a = map( transform_a, cross_a['features'][idx[0]]['geometry']['coordinates'])
        l_b = map( transform_b, cross_b['features'][idx[1]]['geometry']['coordinates'])
        try:
            fplane = faultplane_for_lines(l_a, l_b)
            na = len(l_a)
            point_fplane = []
            for tri in fplane:
                point_tri = []
                for n in tri:
                    if n < na:
                        point_tri.append(l_a[n])
                    else:
                        point_tri.append(l_b[n-na])
                point_fplane.append(point_tri)
            faultplanes[nam] = point_fplane
        except shared.GeometryException as e:
            print "could not interpolate fault %s between %s and %s: " % (nam, cross_a['name'], cross_b['name']), e
    return faultplanes

def faultplane_to_feature(fplane, name):
    """
    Converts a fault plane to a Feature.
    """
    return { 'type': 'Feature',
             'geometry': { 'type': 'Surface', 
                           'coordinates': fplane },
             'properties': { 'name': name } }

def faultplanes_to_featurecollection(fplanes):
    """
    Converts a set of faultplanes to a FeatureCollection
    """
    features = [ faultplane_to_feature(fplane, name) for name, fplane in fplanes.iteritems() ]
    return  { 
        'type': 'FeatureCollection',
        'geology_type': 'faults',
        'features': features,
        'transform': {
                    'type': 'identity'
                }
            }

class AlignedTriangle(object):
    def __init__(self, triangle):
        # First convert the triangle to model coordinate system.
        tr = map(np.array, triangle)
        # Then calculate the normal, but looking in the direction of the model.
        vct = np.cross(tr[1]-tr[0], tr[2]-tr[0])
        
        if vct[2] < 0:
            vct = -vct
        vct /= la.norm(vct)
        
        # Set it in the AlignedTriangle
        self.normal = vct
        self.point  = tr[0]
        
        # Add a polygon.
        self.triangle = Polygon( map( lambda pt: pt[:2], triangle ) )
    
    def __repr__(self):
        return ",".join((str(self.triangle), str(self.point), str(self.normal)))

def align_fault_with(faults, model_point):
    """
    It returns a structure that can be used to know if the triangle is before of after the point.
    """
    ret = {}
    for feature in faults['features']:
        triangles = []
        for triangle in feature['geometry']['coordinates']:
            triangles.append(AlignedTriangle(map(model_point, triangle)))
        ret[feature['properties']['name']] = triangles 
    return ret

def find_fault_plane_intersection(fplane, x0, v1, v2, nv, intersect_line):
    """
    Finds the lines that result after intersecting a fault plane with a plane.
    """
    lines = []
    for tp in fplane:
        ints = []
        for i in range(3):
            # Find the ray parameters.
            ni = (i+1)%3
            p0 = np.array(tp[i])
            p1 = np.array(tp[ni])
            nt = p1-p0
            d  = la.norm(nt)
            nt = nt/d
            
            # Find if the ray and the plane are parallel. (plane normal perpendicular).
            div = np.dot(nv, nt)
            if math.fabs(div) < 1e-9:
                continue
            
            # Find the intersection between the ray and the plane.
            a = np.dot(x0 - p0, nv)/div
            
            # Find if it lies within the line.
            if a > (d + 1e-9) or (a < -1e-9):
                continue
            
            # Find if the point is already in the list.
            p = p0 + nt*a
            for pp in ints:
                if la.norm(pp-p) < 1e-9:
                    break
            else:
                ints.append(p)
        
        # Size of ints must be 2 to append it as lines.
        if len(ints) == 2:
            lines.append(ints)
    
    filt_lines = []
    # Cut the lines to the plane.
    for line in lines:
        res = intersect_line(line)
        filt_lines += res
    
    # Join the lines.
    j_lines = []
    for l in filt_lines:
        # Check that the resulting lines have two coordinates.
        assert len(l) == 2
        for i in range(2):
            p = l[i]
            pc = l[(i+1)%2]
            br = False
            for pl in j_lines:
                fp = pl[0]
                if dist(fp, p) < 1e-9:
                    pl.insert(0, pc)
                    br = True
                    break
                lp = pl[-1]
                if dist(lp, p) < 1e-9:
                    pl.append(pc)
                    br = True
                    break
            if br:
                break
        else:
            j_lines.append(deepcopy(l))
    return j_lines
    
def find_faults_plane_intersection(fplanes, plane):
    """
    Finds the lines that intersect a plane with the fault planes.
    fplanes:    the set of fault planes to intersect with the plane.
    plane:      plane to intersect. It's the four corners of the plane.
    """
    # Find the vectors in the plane that generate the space.
    pn = map(np.array, plane)
    # First is the vector with two corners in the top.
    v1 = pn[1] - pn[0]
    # Second is last point and first point.
    v2 = pn[-1] - pn[0]
    v1 = v1/la.norm(v1)
    
    # Gram-Schmidt is used to make them perpendicular.
    v2 = v2 - np.dot(v1, v2)*v1
    v2 = v2/la.norm(v2)
    
    # Normal orthogonal to both.
    nv = np.cross(v1, v2)
    nv = nv/la.norm(nv)
    x0 = pn[0]
    
    # planep is a function that returns the plane coordinate in terms of v1, v2.
    planep = lambda p: [np.dot(v1, p-pn[0]), np.dot(v2, p-pn[0])]
    
    # Get the polygon terms of v1, v2.
    poly = map(planep, pn)
    shpoly = Polygon(poly)
    
    # Corverts line to point coordinates.
    to_coords = lambda line: map(list, line.coords)
    # Converts returned multigeometry to coordinate lines.
    sep_lines = lambda lg: map( to_coords, shared.shape_list(lg, 'LineString') )
    def intersect_line(line):
        """
        Returns line after cutting it with a polygon.
        """
        p0 = planep(line[0])
        p1 = planep(line[1])
        l = LineString([p0, p1])
        inter = shpoly.intersection(l)
        return sep_lines(inter)
     
    joined_lines = {}
    for nam, fplane in fplanes.iteritems():
        joined_lines[nam] = find_fault_plane_intersection(fplane, x0, v1, v2, nv, intersect_line)
    
    return joined_lines

def find_faults_multiple_planes_intersection(fplanes, planes):
    # All joined lines.
    joined_lines = {}
    
    # Coordinate in x to add to the plane lines.
    startx = 0.0
    for plane in planes:
        lines = find_faults_plane_intersection(fplanes, plane)
        for nam, lines in lines.iteritems():
            for line in lines:
                for point in line:
                    point[0] += startx
            if nam in joined_lines:
                joined_lines[nam] += lines
            else:
                joined_lines[nam] = lines
        # Find the vectors in the plane that generate the space.
        pn = map(np.array, plane)
        # First is the vector with two corners in the top.
        v1 = pn[1] - pn[0]
        # Update the coordinate in x to start.
        startx += la.norm(v1)
