from model import GeologicalModel
from shared import ModelException, TaskException
import random
import numpy as np

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

def calculate_isovalues(model, unit, grid_divisions, bbox, bounded=True):
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
    if bounded:
        signed_distance = lambda x,y,z: model.signed_distance_bounded(unit, (x,y,z))
    else:
        signed_distance = lambda x,y,z: model.signed_distance_unbounded(unit, (x,y,z))

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

def calculate_isosurface(model, unit, grid_divisions, bounded=True, filter_by_normal=False, normal_upwards=True):
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

    bbox = list(model.bbox) # Copy so original is not modified.
    X, Y, Z, sd = calculate_isovalues( model, unit, grid_divisions, bbox, bounded=bounded )
    
    # Check if the surface covers a very small volume, then reduce the bbox to calculate surface.
    bb = check_bbox_surface( sd )
    total_cells = grid_divisions ** 3
    obj_cells = (bb[3]-bb[0])*(bb[4]-bb[1])*(bb[5]-bb[2])
    
    # If the object is at least 8 times smaller than the full bbox, it will benefit lots from a thinner sample.
    if ( float(total_cells) / float(obj_cells) ) > 8.0:
        bbox = [X[bb[0], 0, 0], Y[0, bb[1], 0], Z[0, 0, bb[2]], X[bb[3], 0, 0], Y[0, bb[4], 0], Z[0, 0, bb[5]]]
        X, Y, Z, sd = calculate_isovalues( model, unit, grid_divisions, bbox, bounded=bounded )
    try:
        vertices, simplices, normals, values = measure.marching_cubes(sd, 0)
        del normals
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
    
    ranges = [ X[:,0,0], Y[0,:,0], Z[0,0,:] ]
    
    def real_pt( pt ):
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
    
    vertices = map(real_pt, vertices)
    
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

def save_unit( name, model, unit, grid_divisions, bounded=True, filter_by_normal=False, normal_upwards=False ):
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
    v, s = calculate_isosurface( model, unit, grid_divisions, bounded, filter_by_normal, normal_upwards )
    m = stl_mesh( v, s )
    del v
    del s
    m.save(name)

def triangulate_unit( model, unit, grid_divisions, bounded=True, filter_by_normal=False, normal_upwards=False ):
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
    vertices, triangles = calculate_isosurface( model, unit, grid_divisions, bounded, filter_by_normal, normal_upwards )
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


