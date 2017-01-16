
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

# The utils file contains scripts that can be used by users to 
# make calculations of their models.

from model import GeologicalModel
from shared import ModelException
import random
import numpy as np
import itertools

def generate_simple_grid(query_func, bbox, grid_divisions):
    """
    Returns a uniform grid of sizes grid_divisions x grid_divisions x grid_divisions 
    that covers the given bbox evaluated with the query function.
    """
    
    # Obtain the distance ranges.
    dists = [ bbox[i+3] - bbox[i] for i in range(3) ]
    
    # Obtain the sizes of the grid.
    grid_sizes = [dists[i]/grid_divisions for i in range(3)]
    
    # Generate empty grids and formations arrays.
    grid = np.zeros( ( grid_divisions+1, grid_divisions+1, grid_divisions+1, 3 ) )
    units = np.zeros( ( grid_divisions+1, grid_divisions+1, grid_divisions+1 ), dtype=int)
    
    # Now fill the grid with positions and the forms with formations.
    for i in xrange(grid_divisions+1):
        for j in xrange(grid_divisions+1):
            for k in xrange(grid_divisions+1):
                idx = [i,j,k]
                grid[i,j,k,:] = [ bbox[l] + grid_sizes[l]*idx[l] for l in range(3) ]
                units[i,j,k] = query_func(grid[i,j,k,:])
    
    # Return the generated grid and formations.
    return { 'grid': grid, 'units': units }

# Returns a random number between two points, skewed to the center.
def triangular(a, b):
    n = random.triangular(a,b)
    # Return mid if it's outside or equal to a or b.
    return n

# Returns the a random point between two corners.
def triangular_pt(p0, p1):
    return (triangular(p0[0], p1[0]), triangular(p0[1], p1[1]), triangular(p0[2], p1[2]))

# Simple sum of points.
def ts(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def generate_fdm_grid(query_func, bbox, grid_divisions, max_refinements):
    """
    Generates a grid of points with a FDM like
	refinment method. It first generates a simple grid.
    then it checks if a cell needs refinement. If it does,
    it marks it as a cell to refine.
    Then it goes through every axis, creating planes where the
    cell needs refinements, plus marking the cells as not needing
    refinement.
    """
    
    # First, get the initial grid.
    simple = generate_simple_grid(query_func, bbox, grid_divisions)
    units = simple['units']
    grid = simple['grid']
    
    random.seed()
    
    def should_refine_cell(cell):
        # Finally find a random point and check if the formation is not in paired.
        a = grid[cell]
        b = grid[ts(cell, (1,1,1))]
        g = triangular_pt(a, b) #(a+b)/2.0# 
        ug = query_func(g)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if ug == units[ts(cell, (i,j,k))]:
                        return False
        return True
    
    # Function that evaluates the array in place
    def vquery_func( arr ):
        sh = arr.shape[:2]
        unitm = np.zeros(sh, dtype=int)
        for i in xrange(sh[0]):
            for j in xrange(sh[1]):
                unitm[i,j] = query_func(arr[i,j,:])
        return unitm
    
    # First check which rows or columns have to be refined.
    # last_checked = set()
    for r in range(max_refinements):
        rsiz = ts(units.shape, (-1,-1,-1))
        refine = np.zeros(rsiz, dtype=bool)
        
        for i,j,k in itertools.product(xrange(rsiz[0]), xrange(rsiz[1]), xrange(rsiz[2])):
            refine[i,j,k] = should_refine_cell((i,j,k))
        
        # Check which need refinement in the given axis.
        to_refine = []
        for i in (2,1,0):
            for j in xrange(rsiz[i]-1, -1, -1):
                sl = tuple( slice(None) if k != i else j for k in range(3) )
                sm = refine[sl].sum()
                if sm > 0:
                    to_refine.append((i,j,sm))
                
        if not len(to_refine):
            break
        # Sort by the largest number to refine.
        srt_ref = sorted(to_refine, key=lambda a: a[2], reverse=True)
        
        # Go again through every range, but in the given order by sum.
        to_refine = []
        for axis, idx, osm in srt_ref:
            sl = tuple( slice(None) if k != axis else idx for k in range(3) )
            sm = refine[sl].sum()
            if sm > 0:
                to_refine.append((axis, idx, sm))
                refine[sl] = False
        # Order so insert indexes are correct.
        to_refine.sort(reverse=True) 
        for axis, idx, sm in to_refine:
            gridm =  (np.take(grid, idx, axis=axis)+np.take(grid, idx+1, axis=axis))/2.0
            unitsm = vquery_func(gridm)
            grid =   np.insert(grid, idx+1, gridm, axis=axis)
            units =  np.insert(units, idx+1, unitsm, axis=axis)
        
    return { 'units': units, 'grid': grid }

def generate_octtree_grid(query_func, bbox, grid_divisions, fdm_refine, oct_refine):
    
    simple = generate_fdm_grid(query_func, bbox, grid_divisions, fdm_refine)
    grid  = simple['grid']
    prev_units = simple['units']
    # Convert it to mesh.
    pdims = (grid.shape[0], grid.shape[1], grid.shape[2])
    points = np.zeros((pdims[0]*pdims[1]*pdims[2], 3))
    mdims = ts(pdims, (-1,-1,-1))
    elems = np.zeros((mdims[0]*mdims[1]*mdims[2], 8), dtype=int)
    units = np.zeros((pdims[0]*pdims[1]*pdims[2],), dtype=object)
    for i in xrange(pdims[0]):
        for j in xrange(pdims[1]):
            for k in xrange(pdims[2]):
                pidx = i*pdims[1]*pdims[2] + j*pdims[2] + k
                points[pidx,:] = grid[i,j,k,:]
                units[pidx] = prev_units[i,j,k]
    
    for i in xrange(mdims[0]):
        for j in xrange(mdims[1]):
            for k in xrange(mdims[2]):
                pidx = i*pdims[1]*pdims[2] + j*pdims[2] + k
                elem = []
                for ix in range(2):
                    for jx in range(2):
                        for kx in range(2):
                            elem.append(pidx+ix*pdims[1]*pdims[2]+jx*pdims[2]+kx)
                elems[i*mdims[1]*mdims[2] + j*mdims[2] + k,:] = elem
    
    def should_refine_cell(idx):
        # Check that all the points are equal.
        u = None
        for i in elems[idx,:]:
            if u == None:
                u = units[i]
            elif u != units[i]:
                return True
        
        # Find a random point and check if the formation is different to it.
        a = points[elems[idx][0]]
        b = points[elems[idx][1]]
        g = triangular_pt(a, b)
        ug = query_func(g)
        if ug != u:
            return True
        return False
    
    midpairs = [(0,1),(2,3),(4,5),(6,7),(0,2),(1,3),(4,6),(5,7),(0,4),(1,5),(2,6),(3,7),(0,3),(4,7),(0,5),(2,7),(0,6),(1,7),(0,7)]
    elemsidx = [(0,8,12,20,16,22,24,26),
                (8,1,20,13,22,17,26,25),
                (12,20,2,9,24,26,18,23),
                (20,13,9,3,26,25,23,19),
                (16,22,24,26,4,10,14,21),
                (22,17,26,25,10,5,21,15),
                (24,26,18,23,14,21,6,11),
                (26,25,23,19,21,15,11,7)]
    
    checked = np.zeros((elems.shape[0],), dtype=bool)
    for r in xrange(oct_refine):
        elems_ti = []
        points_ti = []
        units_ti = []
        creamids = {}
        for idx in xrange(elems.shape[0]):
            if not checked[idx] and should_refine_cell(idx):
                newelems = []
                newpoints = []
                currelem = list(elems[idx,:])
                for i,j in midpairs:
                    ie, je = (elems[idx][i],elems[idx][j])
                    if ( ie,je ) in creamids:
                        currelem.append(creamids[ie,je])
                    else:
                        l = points.shape[0]+len(points_ti)+len(newpoints)
                        newpoints.append((points[ie]+points[je])/2.0)
                        currelem.append(l)
                for e in elemsidx:
                    newelems.append([])
                    for ie in e:
                        newelems[-1].append(currelem[ie])
                units_ti += map( query_func, newpoints )
                elems[idx] = newelems[0]
                elems_ti += newelems[1:]
                points_ti += newpoints
            else:
                checked[idx] = True

        units = np.insert(units, units.shape[0], units_ti)
        points = np.insert(points, points.shape[0], np.array(points_ti), axis=0)
        elems = np.insert(elems, elems.shape[0], elems_ti, axis=0)
        ccopy = np.zeros((elems.shape[0],),dtype=bool)
        ccopy[:checked.shape[0]] = checked
        checked = ccopy
    
    return {'points': points, 'elems': elems, 'units': units}
