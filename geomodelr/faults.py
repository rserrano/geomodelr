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
from scipy.spatial.distance import euclidean as dist
import math
import numpy as np
from numpy import linalg as la
from shapely.geometry import Polygon, Point, LineString
from copy import deepcopy

# Finds the lines that result after intersecting a fault plane with a plane.
def find_fault_plane_intersection(fplane, x0, v1, v2, nv, intersect_line):
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

# Finds the lines that intersect a plane with the fault planes.
# fplanes:    the set of fault planes to intersect with the plane.
# plane:      plane to intersect. It's the four corners of the plane.
def find_faults_plane_intersection(fplanes, plane):
    
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

# Finds the intersection but for multiple planes.
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
    return joined_lines
