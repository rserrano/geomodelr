from model import GeologicalModel
from shared import ModelException, TaskException
import random
import numpy as np
import pyopenvdb as vdb
from math import ceil

MESH_AVAILABLE = False

try:
    from scipy.spatial import cKDTree
    from stl import mesh
    from skimage import measure
    MESH_AVAILABLE = True
except ImportError:
    pass

PLOT_AVAILABLE = False
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    if MESH_AVAILABLE:
        PLOT_AVAILABLE = True
except ImportError:
    pass

 
import itertools
import gc

def set_SDFvoxels(grid, vals, iter, step, active_bools):
    """
    Assigns the values of the signed distance function to the LeafNodes of a 8x8x8 block.

    Args:
        (pyopenvdb.FloatGrid) grid: The tree of pyopenvdb that represents the unit surface.

        (numpy.array) vals: values of the signed distance function

        (pyopenvdb.Value) iter: Iterator over active voxels of the tree.

        (int) step: Step to go through the 8x8x8 block.

        (numpy.array) active_bools: the true values of the voxels.

        (float) dy: size of the openvdb voxel in the Y direction.

        (float) dz: size of the openvdb voxel in the Z direction.
    """
        
    min_v = iter.min; max_v = iter.max
    dAccessor = grid.getAccessor()
    for i in range(min_v[0],max_v[0]+1,step):
        I = (i-min_v[0])/step
        for j in range(min_v[1],max_v[1]+1,step):
            J = (j-min_v[1])/step
            for k in range(min_v[2],max_v[2]+1,step):
                K = (k-min_v[2])/step
                if active_bools[I,J,K]:
                    dAccessor.setValueOn((i,j,k), vals[I,J,K])
                else:
                    dAccessor.setActiveState((i,j,k), False)
                    # dAccessor.setValueOff((i,j,k), vals[I,J,K])

    del dAccessor

def sdf_voxels(signed_distance,dx,dy,dz,xi,yi,zi,min_v,max_v,step):
    """
    Calculate the values of the signed distance function of the LeafNodes of a 8x8x8 block.

    Args:
        (function) signed_distance: Signed distance function.

        (float) dx: size of the openvdb voxel in the X direction.

        (float) dy: size of the openvdb voxel in the Y direction.

        (float) dz: size of the openvdb voxel in the Z direction.

        (float) xi: X coordinate of the lower-left corner of the model.

        (float) yi: Y coordinate of the lower-left corner of the model.

        (float) zi: Z coordinate of the lower-left corner of the model.

        (tuple) min_v: (i,j,k) lower-left corner of the 8x8x8 block.

        (tuple) max_v: (i,j,k) upper-right corner of the 8x8x8 block.

        (int) step: Step to go through the 8x8x8 block.
    """
    ni = int(ceil((max_v[0]-min_v[0]+1.0)/step))
    nj = int(ceil((max_v[1]-min_v[1]+1.0)/step))
    nk = int(ceil((max_v[2]-min_v[2]+1.0)/step))
    vals=np.zeros((ni,nj,nk))

    for i in range(min_v[0],max_v[0]+1,step):
        I = (i-min_v[0])/step
        for j in range(min_v[1],max_v[1]+1,step):
            J = (j-min_v[1])/step
            for k in range(min_v[2],max_v[2]+1,step):
                K = (k-min_v[2])/step
                vals[I,J,K] = signed_distance(dx*i+xi,dy*j+yi,dz*k+zi)
    return vals

def grid_openvdb(signed_distance, grid_divisions, bbox, ndelta):
    """
    Creates a pyopenvdb.FloatGrid of the unit surface.

    Args:
        (function) signed_distance: Signed distance function.

        (int) grid_divisions: The divisions of the grid in the X, Y and Z directions.

        (list) bbox: the values of [minx, miny, minz, maxx, maxy, maxz].

        (int) ndelta: number of voxels to use to increase the bounding box.
    """

    # Fixed parameters

    nx, ny, nz = grid_divisions
    steps = 2

    dx = (bbox[3] - bbox[0])/(nx-1)
    dy = (bbox[4] - bbox[1])/(ny-1)
    dz = (bbox[5] - bbox[2])/(nz-1)

    xi = bbox[0]-ndelta*dx; yi = bbox[1]-ndelta*dy; zi = bbox[2]-ndelta*dz
    nx+=2*ndelta; ny+=2*ndelta; nz+=2*ndelta

    bg = 1.5*np.sqrt(dx**2 + dy**2 + dz**2)

    grid = vdb.FloatGrid(background=bg)
    grid.fill(min=(0,0,0), max=(nx-1, ny-1, nz-1), value=bg)
    grid.gridClass = vdb.GridClass.LEVEL_SET

    bg = grid.background

    change = True
    while change:
        change = False
        for iter in grid.iterOnValues():

            if iter.depth==3:
                if iter.value==bg:
                    i,j,k = iter.min
                    dist = signed_distance(dx*i+xi,dy*j+yi,dz*k+zi)                        
                    if abs(dist)<bg:
                        iter.value = dist
                    else:
                        iter.active = False

            elif iter.depth==2:
                vals = sdf_voxels(signed_distance,dx,dy,dz,xi,yi,zi,iter.min,iter.max,steps)
                active_bools = np.abs(vals)<bg
                if np.any(active_bools):
                    set_SDFvoxels(grid,vals,iter,steps,active_bools)
                    change = True
                else:
                    iter.active = False
                    # iter.value = np.sign(vals[2,2,2])*bg
            else:
                i,j,k = iter.max
                dist = signed_distance(dx*i+xi,dy*j+yi,dz*k+zi)
                dAccessor = grid.getAccessor()
                if abs(dist)<bg:
                    dAccessor.setValueOn((i,j,k), dist)
                else:
                    dAccessor.setActiveState((i,j,k), False)
                    # dAccessor.setValueOff((i,j,k), dist)
                del dAccessor   
                change = True

    return (grid,dx,dy,dz,xi,yi,zi,nx,ny,nz)

def check_bbox(grid,ndelta,grid_size):
    """
    Checks if the bounding box of the object modeled is very small. When a geological unit covers a very
    small part of the model, it needs to be refined. The new bbox of the unit is returned to check that
    case.

    Args:
        (pyopenvdb.FloatGrid) grid: The tree of pyopenvdb that represents the unit surface.

    """

    mn = [ float("inf"),  float("inf"),  float("inf")]
    mx = [-float("inf"), -float("inf"), -float("inf")]

    for iter in grid.citerOnValues():
        if iter.value<0:
            i,j,k = iter.min
            mn[0] = min( float(i), mn[0] )
            mn[1] = min( float(j), mn[1] )
            mn[2] = min( float(k), mn[2] )

            mx[0] = max( float(i), mx[0] )
            mx[1] = max( float(j), mx[1] )
            mx[2] = max( float(k), mx[2] )

    if mn[0] > mx[0]:
        raise TaskException("This model does not contain the unit or the sample is too coarse")

    mn = [ max(mn[i]-2,ndelta) for i in xrange(3) ]
    mx = [ min(mx[i]+2,ndelta+grid_size[i]-1) for i in xrange(3) ]
    return mn + mx

def grid_to_mesh(grid,xi,yi,zi,dx,dy,dz,isovalue,adaptive):

    """
    Creates a triangulation of the unit using a tree of pyopenvdb.
    
    Args:
        (pyopenvdb.FloatGrid) grid: The tree of pyopenvdb that represents the unit surface.

        (float) xi: X coordinate of the lower-left corner of the model.

        (float) yi: Y coordinate of the lower-left corner of the model.

        (float) zi: Z coordinate of the lower-left corner of the model.

        (float) dx: size of the openvdb voxel in the X direction.

        (float) dy: size of the openvdb voxel in the Y direction.

        (float) dz: size of the openvdb voxel in the Z direction.
    """

    grid.transform.scale((dx,dy,dz))
    grid.transform.translate((xi,yi,zi))

    points, trians, quads = grid.convertToPolygons(isovalue=isovalue,adaptivity=adaptive)
    del grid
    nquads = len(quads)
    triangles = np.zeros((2*nquads,3),dtype='int32')

    for k in range(nquads):
        triangles[2*k] = quads[k,[0,1,2]]
        triangles[2*k+1] = quads[k,[2,3,0]]

    if len(trians)>0:
        triangles = np.concatenate((triangles, trians.astype('int32')))

    return (points,triangles)

# ==================================================
#                   Begin: Resample
# ===================================================

def eval_linear(vals,xc,yc,zc,dx,dy,dz,x,y,z):

    x = (x-xc)*2.0/dx -1.0; y = (y-yc)*2.0/dy -1.0; z = (z-zc)*2.0/dz -1.0
    subs_y = 1.0 - y; sum_y = 1.0 + y
    subs_z = 1.0 - z; sum_z = 1.0 + z

    v1 = subs_y*(vals[1] + vals[0] + x*(vals[1] - vals[0]))
    v2 = sum_y *(vals[3] + vals[2] + x*(vals[3] - vals[2]))
    v3 = subs_y*(vals[5] + vals[4] + x*(vals[5] - vals[4]))
    v4 = sum_y *(vals[7] + vals[6] + x*(vals[7] - vals[6]))

    return (subs_z*(v1 + v2) + sum_z*(v3 + v4))/8.0

def i_checker(vals,dl):
    eps = -1e-5
    b1 = (vals[0]*vals[1]>eps) & (abs(vals[0]+vals[1])<dl)
    b2 = (vals[2]*vals[3]>eps) & (abs(vals[2]+vals[3])<dl)
    b3 = (vals[4]*vals[5]>eps) & (abs(vals[4]+vals[5])<dl)
    b4 = (vals[6]*vals[7]>eps) & (abs(vals[6]+vals[7])<dl)
    return b1|b2|b3|b4

def j_checker(vals,dl):
    eps = -1e-5
    b1 = (vals[0]*vals[2]>eps) & (abs(vals[0]+vals[2])<dl)
    b2 = (vals[1]*vals[3]>eps) & (abs(vals[1]+vals[3])<dl)
    b3 = (vals[4]*vals[6]>eps) & (abs(vals[4]+vals[6])<dl)
    b4 = (vals[5]*vals[7]>eps) & (abs(vals[5]+vals[7])<dl)
    return b1|b2|b3|b4

def k_checker(vals,dl):
    eps = -1e-5
    b1 = (vals[0]*vals[4]>eps) & (abs(vals[0]+vals[4])<dl)
    b2 = (vals[1]*vals[5]>eps) & (abs(vals[1]+vals[5])<dl)
    b3 = (vals[2]*vals[6]>eps) & (abs(vals[2]+vals[6])<dl)
    b4 = (vals[3]*vals[7]>eps) & (abs(vals[3]+vals[7])<dl)
    return b1|b2|b3|b4

def get_ijk_values(dAccessor,i,j,k,dx,dy,dz,bg_resample):

    vals = np.zeros(8)
    states = np.zeros(8,dtype='bool')
    vals[0],states[0] = dAccessor.probeValue((i,j,k))
    vals[1],states[1] = dAccessor.probeValue((i+1,j,k))
    vals[2],states[2] = dAccessor.probeValue((i,j+1,k))
    vals[3],states[3] = dAccessor.probeValue((i+1,j+1,k))
    vals[4],states[4] = dAccessor.probeValue((i,j,k+1))
    vals[5],states[5] = dAccessor.probeValue((i+1,j,k+1))
    vals[6],states[6] = dAccessor.probeValue((i,j+1,k+1))
    vals[7],states[7] = dAccessor.probeValue((i+1,j+1,k+1))

    active_state = np.all(states)
    if active_state:
        if i_checker(vals,dx) | j_checker(vals,dy) | k_checker(vals,dz):
            return (vals,active_state,True)

    if np.all(vals>bg_resample) | np.all(vals<-bg_resample):
        active_state = False

    return (vals,active_state,False)

def fill_block(signed_distance,dAcc_resample,dAcc,bg_resample,rx,ry,rz,i,j,k,vals,dx,dy,dz,xi,yi,zi,bool_resample):

    xc = xi + i*dx; yc = yi + j*dy; zc = zi + k*dz

    if bool_resample:
        resample_distance = lambda x,y,z: signed_distance(x,y,z)
    else:
        resample_distance = lambda x,y,z: eval_linear(vals,xc,yc,zc,dx,dy,dz,x,y,z)

    def is_corner(iterV,res):
        if np.mod(iterV,res)==0:
            return True, int(iterV/res)
        return False, -1
    
    for k_block in xrange(k*rz,(k+1)*rz + 1):
        z_pt = zi + k_block*dz/rz
        z_corner, kc = is_corner(k_block,rz)
        
        for j_block in xrange(j*ry,(j+1)*ry + 1):
            y_pt = yi + j_block*dy/ry
            y_corner, jc = is_corner(j_block,ry)
            
            for i_block in xrange(i*rx,(i+1)*rx + 1):
                x_pt = xi + i_block*dx/rx
                x_corner, ic = is_corner(i_block,rx)

                ijk = (i_block,j_block,k_block)
                if dAcc_resample.getValue(ijk)==bg_resample:

                    if x_corner & y_corner & z_corner:
                        dist = dAcc.getValue((ic,jc,kc))
                    else:
                        dist = resample_distance(x_pt,y_pt,z_pt)

                    if abs(dist)<bg_resample:
                        dAcc_resample.setValueOn(ijk, dist)
                    else:
                        dAcc_resample.setValueOff(ijk, dist)

def resample_openvdb_grid(signed_distance,grid,grid_sample,grid_size,xyz_corner,num_res):

    dx,dy,dz = grid_sample
    nx,ny,nz = grid_size
    xi,yi,zi = xyz_corner
    rx,ry,rz = num_res

    bg_resample = 1.5*np.sqrt((dx/rx)**2 + (dy/ry)**2 + (dz/rz)**2)

    resam_grid = vdb.FloatGrid(background=bg_resample)
    resam_grid.gridClass = vdb.GridClass.LEVEL_SET
    bg_resample = resam_grid.background

    dAcc = grid.getConstAccessor()
    dAcc_resample = resam_grid.getAccessor()

    for iter in grid.citerOnValues():
        i,j,k = iter.min
        vals,state,bool_resample = get_ijk_values(dAcc,i,j,k,dx,dy,dz,bg_resample)
        if state:
            fill_block(signed_distance,dAcc_resample,dAcc,bg_resample,rx,ry,rz,i,j,k,vals,dx,dy,dz,xi,yi,zi,bool_resample)
    
    del dAcc, dAcc_resample, grid

    dx,dy,dz = (dx/rx, dy/ry, dz/rz)
    nx,ny,nz = (nx*rx-1, ny*ry-1, nz*rz-1)
    return (resam_grid,dx,dy,dz,xi,yi,zi,nx,ny,nz)

# ==================================================
#                   End: Resample
# ===================================================

def calculate_isovalues( signed_distance, unit, grid_divisions, bbox ):
    """
    Calculates a grid of isovalues in a bounding box to be used by method 
    skimage.measure to create a triangulation of the unit.
    
    Args:
        (geomodelr.model.GeologicalModel) geolojson: The Geological model created form a Geological JSON object.

        (int) grid_divisions: The divisions of the grid in the X, Y and Z directions.

        (list) bbox: the values of [minx, miny, minz, maxx, maxy, maxz].

        (bool) bounded: whether the object will be a solid, (the bounded signed distance will be used), or it will be an open surface.
    """

    dx = (bbox[3]-bbox[0])/grid_divisions
    dy = (bbox[4]-bbox[1])/grid_divisions
    dz = (bbox[5]-bbox[2])/grid_divisions
    
    bbox[0] -= 2*dx
    bbox[1] -= 2*dy
    bbox[2] -= 2*dz
    
    bbox[3] += 2*dx
    bbox[4] += 2*dy
    bbox[5] += 2*dz
    
    X,Y,Z = np.mgrid[bbox[0]:bbox[3]:dx, bbox[1]:bbox[4]:dy, bbox[2]:bbox[5]:dz]
    
    vsigned_distance = np.vectorize(signed_distance, otypes=[np.float])
    sd = vsigned_distance( X, Y, Z )
    return (X, Y, Z, sd)

def check_bbox_surface( sd ):
    """
    Checks if the bounding box of the object modeled is very small. When a geological unit covers a very
    small part of the model, it needs to be refined. The new bbox of the unit is returned to check that
    case.
    
    Args:
        (numpy.array) sd: The grid of signed distances.
    """
    mn = [ float("inf"),  float("inf"),  float("inf")]
    mx = [-float("inf"), -float("inf"), -float("inf")]
    
    for i in xrange(sd.shape[0]):
        for j in xrange(sd.shape[1]):
            for k in xrange(sd.shape[2]):
                if sd[i,j,k] < 0:
                    mn[0] = min( i, mn[0] )
                    mn[1] = min( j, mn[1] )
                    mn[2] = min( k, mn[2] )
                    
                    mx[0] = max( i, mx[0] )
                    mx[1] = max( j, mx[1] )
                    mx[2] = max( k, mx[2] )
    
    if mn[0] > mx[0]:
        raise TaskException("This model does not contain the unit or the sample is too coarse")
    # The bbox is the first possitive.
    mn = [ max(mn[i]-1, 0) for i in xrange(3) ]
    mx = [ min(mx[i]+1, sd.shape[i]-1) for i in xrange(3) ]
    return mn + mx

def calculate_normals(vertices, simplices):
    """
    Calculates the normals for all the simplices, (skimage returns normals per vertex).
    
    Args:
        (numpy.array) vertices: The vertices returned by the marching cubes algorithm.

        (numpy.array) simplices: The simplices (triangles) returned by the marching cubes algorithm.

    """
    normals = np.zeros(simplices.shape, dtype=float)
    v1 = vertices[simplices[:,1]] - vertices[simplices[:,0]]
    v2 = vertices[simplices[:,2]] - vertices[simplices[:,0]]
    return np.cross(v1, v2)

def calculate_isosurface(model, unit, grid_divisions, bounded=True, filter_by_normal=False, normal_upwards=True, aligned=False ):
    """
    Calculates an isosurface of a unit. It uses a signed distance and an isosurface algorithm present in skimage.measure.
    
    Args:
        (geomodelr.model.GeologicalModel) model: The geomodelr geological model.

        (unicode) unit: The unit to calculate the isosurface for.

        (int) grid_divisions: The number of divisions for all the axes.

        (bool) bounded: calculates the surface using the bounding box. This will result in a solid.

        (bool) filter_by_normal: filters by the normal of each triangle, depending on the normal_upwards argument.

        (bool) normal_upwards: if filter_by_normal is True, filters the triangles depending on its normal. It returns only
        the triangles pointing upwards if it's True, otherwise it returns the triangles pointing downwards.
    
    Returns:
        (list) vertices: The list of vertices.

        (list) triangles: The list of triangles indexes to vertexes.
    """
    bbox = list(model.bbox) if not aligned else list(model.abbox) # Copy so original is not modified.
    if aligned:
        if bounded:
            signed_distance = lambda x,y,z: model.signed_distance_bounded_aligned(unit, (x,y,z))
        else:
            signed_distance = lambda x,y,z: model.signed_distance_unbounded_aligned(unit, (x,y,z))
    else:
        if bounded:
            signed_distance = lambda x,y,z: model.signed_distance_bounded(unit, (x,y,z))
        else:
            signed_distance = lambda x,y,z: model.signed_distance_unbounded(unit, (x,y,z))
    
    # X, Y, Z, sd = calculate_isovalues( signed_distance, unit, grid_divisions, bbox )
    
    # Check if the surface covers a very small volume, then reduce the bbox to calculate surface.
    # bb = check_bbox_surface( sd )
    # total_cells = grid_divisions ** 3
    # obj_cells = (bb[3]-bb[0])*(bb[4]-bb[1])*(bb[5]-bb[2])
    
    # If the object is at least 8 times smaller than the full bbox, it will benefit lots from a thinner sample.
    # if ( float(total_cells) / float(obj_cells) ) > 8.0:
    #     bbox = [X[bb[0], 0, 0], Y[0, bb[1], 0], Z[0, 0, bb[2]], X[bb[3], 0, 0], Y[0, bb[4], 0], Z[0, 0, bb[5]]]
    #     X, Y, Z, sd = calculate_isovalues( signed_distance, unit, grid_divisions, bbox )
    
    ndelta = 4
    grid_size = (grid_divisions,grid_divisions,grid_divisions)
    grid,dx,dy,dz,xi,yi,zi,nx,ny,nz = grid_openvdb(signed_distance, grid_size, bbox, ndelta)
    bb = check_bbox(grid,ndelta,grid_size)
    obj_cells = (bb[3]-bb[0]-2*ndelta)*(bb[4]-bb[1]-2*ndelta)*(bb[5]-bb[2]-2*ndelta)

    # If the object is at least 8 times smaller than the full bbox, it will benefit lots from a thinner sample.
    if float(nx*ny*nz)/obj_cells>8:
        del grid
        bbox = [bb[0]*dx+xi, bb[1]*dy+yi, bb[2]*dz+zi, bb[3]*dx+xi, bb[4]*dy+yi, bb[5]*dz+zi]
        grid,dx,dy,dz,xi,yi,zi,nx,ny,nz = grid_openvdb(signed_distance, grid_size, bbox, ndelta)
    try:
        # if unit==u'Cardium Sand 3 top':
        # # vertices, simplices, normals, values = measure.marching_cubes(sd, 0)
        # # del normals
        #     grid_sample = (dx,dy,dz)
        #     grid_size = (nx,ny,nz)
        #     xyz_corner = (xi,yi,zi)
        #     num_res= (6,2,4)
        #     grid,dx,dy,dz,xi,yi,zi,nx,ny,nz = resample_openvdb_grid(signed_distance,grid,grid_sample,grid_size,xyz_corner,num_res)
        vertices, simplices = grid_to_mesh(grid,xi,yi,zi,dx,dy,dz,1e-10,0.12)

        if not(bounded):
            bbox[0] -= ndelta*dx
            bbox[1] -= ndelta*dy
            bbox[2] -= ndelta*dz

            bbox[3] += ndelta*dx
            bbox[4] += ndelta*dy
            bbox[5] += ndelta*dz

            def check_outsidePT(pt):
                eps = 1e-10
                x_out = (pt[0]<=(bbox[0]+eps+dx)) | (pt[0]>=(bbox[3]-eps-dx))
                y_out = (pt[1]<=(bbox[1]+eps+dy)) | (pt[1]>=(bbox[4]-eps-dy))
                z_out = (pt[2]<=(bbox[2]+eps+dz)) | (pt[2]>=(bbox[5]-eps-dz))
                return x_out|y_out|z_out

            def check_outsideTR(tri):
                pts = vertices[tri]
                return not(np.all(map(check_outsidePT,pts)))

            simplices = simplices[ map(check_outsideTR,simplices) ] 

        # vertices = vertices.tolist()

    except ValueError:
        raise TaskException("This model does not contain the unit or the sample is too coarse")
    
    if filter_by_normal:
        normals = calculate_normals( vertices, simplices )
        # First, remove triangles with the normal incorrect.
        if normal_upwards:
            fltr = normals[:,2] >= 0
        else:
            fltr = normals[:,2] <= 0
        simplices = simplices[fltr,:]
        
        # Remove vertices that don't belong to any triangle.
        lverts = vertices.shape[0]
        belongs = np.zeros(lverts, dtype=bool)
        for i in xrange(3):
            belongs[simplices[:,i]] = True
        vertices = vertices[belongs]
        
        # Renum
        idx = 0
        renum = np.zeros(lverts, dtype=int)
        for i in xrange(renum.shape[0]):
            if belongs[i]:
                renum[i] = idx
                idx += 1
        
        for i in xrange(simplices.shape[0]):
            for j in xrange(3):
                simplices[i,j] = renum[simplices[i,j]]
    
    # ranges = [ X[:,0,0], Y[0,:,0], Z[0,0,:] ]
    
    def real_pt_simple( pt ):
        gr = map( lambda c: ( int(c), c-int(c) ), pt )
        outp = []
        for i in range(3):
            assert gr[i][0] < len(ranges[i])
            if gr[i][0]+1 == len(ranges[i]):
                c = ranges[i][gr[i][0]]
            else:
                c = ranges[i][gr[i][0]]*(1-gr[i][1]) + ranges[i][gr[i][0]+1]*gr[i][1]
            outp.append( c )
        return outp
    
    if aligned:
        # real_pt = lambda p: model.inverse_point( real_pt_simple( p ) )
        real_pt = lambda p: model.inverse_point( p )
    else:
        real_pt = lambda p: p
    
    vertices = map(real_pt, vertices.tolist())
    
    return vertices, simplices.tolist()

def plot_unit( model, unit, grid_divisions, bounded=True, filter_by_normal=False, normal_upwards=False ):
    """
    Plots a unit previously modeled with calculate_isosurface.
    
    Args:
        (geomodelr.model.GeologicalModel) model: The geomodelr geological model.

        (unicode) unit: The unit to calculate the isosurface for.

        (int) grid_divisions: The number of divisions for all the axes.

        (bool) bounded: calculates the surface using the bounding box. This will result in a solid.

        (bool) filter_by_normal: filters by the normal of each triangle, depending on the normal_upwards argument.

        (bool) normal_upwards: if filter_by_normal is True, filters the triangles depending on its normal. It returns only

        the triangles pointing upwards if it's True, otherwise it returns the triangles pointing downwards.
    
    """
    assert MESH_AVAILABLE, "To be able to plot units, you need the following packages installed: scipy, numpy-stl, scikit-image."
    assert PLOT_AVAILABLE, "To be able to plot units, you need the following packages installed: matplotlib."

    vertices, simplices = calculate_isosurface(model, unit, grid_divisions, bounded, filter_by_normal, normal_upwards)
    x,y,z = zip(*vertices)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(x, y, z, triangles=simplices)
    
    fig.show()
    raw_input("Enter to close...")

def stl_mesh( vertices, simplices ):
    """
    Creates a numpy-stl mesh from a sets of vertices and triangles.

    Args:
        (list) vertices: vertices of the mesh.

        (list) simplices: triangles of the mesh.
    Returns:
        (stl.mesh.Mesh): a numpy-stl mesh.
    """
    # fix_solid(vertices, simplices)
    m = mesh.Mesh(np.zeros(len(simplices), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(simplices):
        for j in range(3):
            m.vectors[i][j] = vertices[f[j]]
    return m

def save_unit( name, model, unit, grid_divisions, bounded=True, filter_by_normal=False, normal_upwards=False, aligned=False ):
    """
    Saves the wireframe of a geological unit to the disk. It uses a marching cubes and a signed distance from
    the model.

    Args:
        (str) name: the name to save the STL file.

        (geomodelr.model.GeologicalModel) model: the model to be queried.

        (unicode) unit: the unit to be meshed.

        (int) grid_divisions: the number of divisions of the grid to mesh the object.

        (bool) bounded: whether this surface is bounded by the bbox or only by the topography.

        (bool) filter_by_normal: whether to filter this mesh by normal to the surface. Useful 
        if you want to see the top or bottom of your formation.

        (bool) normal_upwards: if filter_by_normal is True, whether you want the triangles that look up or
        the triangles that look down.
    """
    assert MESH_AVAILABLE, "To be able to save units, you need the following packages installed: scipy, numpy-stl, scikit-image."
    v, s = calculate_isosurface( model, unit, grid_divisions, bounded, filter_by_normal, normal_upwards, aligned=False )
    m = stl_mesh( v, s )
    del v
    del s
    m.save(name)

def triangulate_unit( model, unit, grid_divisions, bounded=True, filter_by_normal=False, normal_upwards=False, aligned=False ):
    """
    Triangulates a geological unit and returns it for further processing, (or saving it to the database).
    It uses a marching cubes and a signed distance from the model.

    Args:
        (geomodelr.model.GeologicalModel) model: the model to be queried.

        (unicode) unit: the unit to be meshed.

        (int) grid_divisions: the number of divisions of the grid to mesh the object.

        (bool) bounded: whether this surface is bounded by the bbox or only by the topography.

        (bool) filter_by_normal: whether to filter this mesh by normal to the surface. Useful 
        if you want to see the top or bottom of your formation.

        (bool) normal_upwards: if filter_by_normal is True, whether you want the triangles that look up or
        the triangles that look down.
    Returns:
        (dict): the triangulated unit with a few useful properties.
    """
    
    assert MESH_AVAILABLE, "To be able to triangulate units, you need the following packages installed: scipy, numpy-stl, scikit-image."
    vertices, triangles = calculate_isosurface( model, unit, grid_divisions, bounded, filter_by_normal, normal_upwards, aligned )
    m = stl_mesh( vertices, triangles )
    volume, cog, inertia = m.get_mass_properties()
    del m
    
    mins = [ float('inf')] * 3
    maxs = [-float('inf')] * 3
    
    for v in vertices:
        for i in range(3):
            mins[i] = min( mins[i], v[i] )
            maxs[i] = max( maxs[i], v[i] )
    
    area = measure.mesh_surface_area( np.array(vertices), triangles )
    props = { 'properties': { 'center_of_gravity': cog.tolist(), 'bbox': mins + maxs, 'surface_area': area }, 
              'mesh': { 'vertices': vertices, 'triangles': triangles } }
    # Volume is incorrect and misleading when the feature is unbounded.
    if bounded:
        props['properties']['volume'] = volume
    return props


