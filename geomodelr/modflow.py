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

import numpy as np
import flopy as fp
from math import ceil,floor

from datetime import datetime
import time
from copy import deepcopy

class LENGTH_UNIT:
    UNDEFINED = 0
    FEET =1
    METERS = 2
    CENTIMETERS = 3

class TIME_UNIT:
    UNDEFINED = 0
    SECONDS = 1
    MINUTES = 2
    HOURS = 3
    DAYS = 4
    YEARS = 5

class ALGORITHM:
    REGULAR = u'regular'
    ADAPTIVE = u'adaptive'

def create_modflow_inputs(name, model, units_data,
    length_units=LENGTH_UNIT.METERS, rows=100, cols=100, layers=100,
    bbox=None, angle=20, dz_min = 1.0, time_units=TIME_UNIT.SECONDS, 
    algorithm=ALGORITHM.REGULAR,faults_data={},faults_method=u'regular'):
    """
    Generates the DIS, BAS, LPF and NAM files, which are used by classical
    MODFLOW processors. The user has to import the NAM file from his MODFLOW
    program. 
    
    Args:
        (string) name: name of generated files.

        (GeologicalModel) model: geomodlr model to work on.

        (dict) units_data: the dictionary keys correpond to the unit names
        and the items are tuples with 4 values: (Kh_x,ani,Kv,i_bound).
        Kh_x and Kv are the horizontal and vertical hydraulic conductivy
        respectively and ani variable is the horizontal anysotropy (rate between
        Kh_y and Kh_x, i.e., ani=Kh_y/Kk_x. Finally, i_bound is equal to 1 
        if the respective unit is active and 0 otherwise.

        (int) length_units: length units (see: class LENGTH_UNIT)

        (int) rows: number of rows in the MODFLOW model.

        (int) cols: number of cols in the MODFLOW model.

        (int) layers: number of layers in the MODFLOW model.

        (list) bbox: the bounding box to search in.

        (float) angle: this angle (degrees) is used in the adative algorithm.

        (float) dz_min: minimum allowed distance between layers.

        (int) time_units: time units (see: class TIME_UNIT)

        (string) algorithm: grid generation algorithm (see: class ALGORITHM)

        (dict) faults_data: the dictionary keys correpond to the fault names
        and the items are tuples with 4 values: (Kh_x,ani,Kv,i_bound).

        (unicode) faults_method: it is equal to "regular" if the hydraulic conductivity
        values of the fault are aligned with the XYZ axes (princial axes) or
        it is equal to "vectorial" if the hydraulic conductivity values of
        the fault are aligned with the fault plane.

    """

    if (bbox is None):
        bbox = model.bbox

    X_inf = bbox[0]
    X_sup = bbox[3]

    Y_inf = bbox[1]
    Y_sup = bbox[4]

    dX = (X_sup - X_inf)/cols
    dY = (Y_sup - Y_inf)/rows

    #X_vec = np.linspace(X_inf + dX/2., X_sup - dX/2.,cols)
    #Y_vec = np.linspace(Y_inf + dY/2., Y_sup - dY/2.,rows)
 
    Z_top = np.zeros((rows,cols))

    bottom_min = bbox[2] + 1E-5

    if (len(faults_data))>0 and (faults_method!=u'regular'):
        faults_keys = faults_data.keys()
        faults_v = faults_data.values()
        faults_data_aux = { faults_keys[k]: (faults_v[k][0],np.random.rand(),faults_v[k][1],faults_v[k][2]) for k in range(len(faults_data)) }
        faults_data = deepcopy(faults_data_aux)
        del(faults_data_aux)

    # Define Z-top
    for i in np.arange(rows):
        yp = Y_sup - (2*i+1)*dY/2.0
        for j in np.arange(cols):
            xp = X_inf + (2*j+1)*dX/2.0
            Z_top[i,j] = model.height((xp, yp))


    Z_top_min = np.min(Z_top)

    if ((Z_top_min-bottom_min)/layers < dz_min):
        layers = int((Z_top_min-bottom_min)/dz_min)

    if (algorithm == 'regular'):

        Z_bottoms=regular_grid(model, rows,cols,layers,Z_top,X_inf,Y_sup,dX,dY,
            bottom_min,units_data)
    
    elif (algorithm == 'adaptive'):

        Z_bottoms,layers=adaptive_grid(model,rows,cols,layers,Z_top,X_inf,Y_sup,dX,dY,
            bottom_min,units_data,angle,dz_min)

    K_hor, K_anisotropy_hor, K_ver, I_bound,chani_var=set_unit_properties(model,
        units_data,Z_top,Z_bottoms,rows,cols,layers,X_inf,Y_sup,dX,dY,faults_data)

    if len(faults_data)>0:
        faults_intersections(model.faults,faults_data,rows,cols,layers,Z_top,Z_bottoms,X_inf,Y_inf,
            dX,dY,K_hor,K_ver,K_anisotropy_hor,I_bound,faults_method)

    if not(np.isscalar(I_bound)):
        cells_checker(I_bound,rows,cols,layers)

    #  ------- Flowpy Packages ----
    # Grid

    #geo=fp.utils.reference.SpatialReference(delr=dX*np.ones(cols),delc=dY*np.ones(rows),
        #lenuni=length_units, xll=X_inf, yll=Y_inf,units='meters',epsg=3116)

    mf_handle = fp.modflow.mf.Modflow(modelname=name,namefile_ext='nam',
        version='mf2005')
    # Variables for the Dis package
    dis = fp.modflow.ModflowDis(mf_handle,nlay=layers, nrow=rows, ncol=cols,
        top=Z_top, botm=Z_bottoms, delc=dY, delr=dX, xul=X_inf, yul=Y_sup,
        itmuni=time_units, lenuni=length_units, proj4_str='EPSG:3116')

    # Variables for the BAS package
    bas = fp.modflow.ModflowBas(mf_handle,ibound = I_bound)

    # Add LPF package to the MODFLOW model
    lpf = fp.modflow.ModflowLpf(mf_handle, chani=chani_var, hk=K_hor,
        vka=K_ver, hani=K_anisotropy_hor,laytyp=np.ones(layers,dtype=np.int32))#

    mf_handle.write_input()

    output = {'num_layers': layers}
    return(output)

# ===================== AUXILIAR FUNCTIONS ========================

def neighboring_cells(i,j,k,layers,rows,cols,I_bound):

    """
    checks if a cell has only one cell connected at the most.
    
    Args:
        (int) i,j,k : position of the given cell.

        (int) rows: number of rows in the MODFLOW model.

        (int) cols: number of cols in the MODFLOW model.

        (int) layers: number of layers in the MODFLOW model.

        (numpy.array) I_bound: 3D-block array of the ibound data.
    """

    cell_ibound = np.zeros(6)
    if (j>0):
        cell_ibound[0] = I_bound[k,i,j-1]
    if (j<cols-2):
        cell_ibound[1] = I_bound[k,i,j+1]
    if (i>0):
        cell_ibound[2] = I_bound[k,i-1,j]
    if (i<rows-2):
        cell_ibound[3] = I_bound[k,i+1,j]
    if (k>0):
        cell_ibound[4] = I_bound[k-1,i,j]
    if (k<layers-2):
        cell_ibound[5] = I_bound[k+1,i,j]

    return np.sum(cell_ibound)

def cells_checker(I_bound,rows,cols,layers):
    """
    deactivates some cells which are hydraulically isolated or have
    just one connected cell.
    
    Args:
        (numpy.array) I_bound: 3D-block array of the ibound data.

        (int) rows: number of rows in the MODFLOW model.

        (int) cols: number of cols in the MODFLOW model.

        (int) layers: number of layers in the MODFLOW model.
    """
    for i in range(rows):
        for j in range(cols):
            for k in range(layers):
                if neighboring_cells(i,j,k,layers,rows,cols,I_bound) < 2:
                    I_bound[k,i,j]=0


def faults_intersections(faults,faults_data,rows,cols,layers,Z_top,Z_bottoms,X_inf,Y_inf,
    dX,dY,K_hor,K_ver,K_anisotropy_hor,I_bound,faults_method):

    """
    Intersecets the faults of the model with cells of the grid.
    
    Args:
        (dict) faults: dictionary of the faults geometry.

        (dict) faults_data: the dictionary keys correpond to the fault names
        and the items are tuples with 4 values: (Kh_x,ani,Kv,i_bound).

        (int) rows: number of rows in the MODFLOW model.

        (int) cols: number of cols in the MODFLOW model.

        (int) layers: number of layers in the MODFLOW model.

        (numpy.array) Z_top: topogrphy of the grid model.

        (numpy.array) Z_bottoms: bottoms of the layers.

        (float) X_inf: x coordinate of lower left corner of the grid.

        (float) Y_inf: y coordinate of lower left corner of the grid.

        (float) dX: spacings along a col.

        (float) dY: spacings along a row.

        (numpy.array) K_hor: 3D-block array of the horizontal hydraulic conductivity.

        (numpy.array) K_ver: 3D-block array of the vertical hydraulic conductivity.

        (numpy.array) K_anisotropy_hor: 3D-block array of the horizontal anisotropy.

        (numpy.array) I_bound: 3D-block array of the ibound data.

        (unicode) faults_method: it is equal to "regular" if the hydraulic conductivity
        values of the fault are aligned with the XYZ axes (princial axes) or
        it is equal to "vectorial" if the hydraulic conductivity values of
        the fault are aligned with the fault plane.


    """

    corner  = np.array([X_inf,Y_inf,0.0])
    count_tri = -1
    cell =np.array([[-dX/2,-dY/2,0],[dX/2,-dY/2,0],[dX/2,dY/2,0],[-dX/2,dY/2,0],
        [-dX/2,-dY/2,0],[dX/2,-dY/2,0],[dX/2,dY/2,0],[-dX/2,dY/2,0]])
    h_xyz = np.array([dX,dY,0])/2.0

    epsilon = 1e-5
    
    anisotropy_bool = not(np.isscalar(K_anisotropy_hor))
    ibound_bool = not(np.isscalar(I_bound))

    for name,fault in faults.iteritems():
        Data = faults_data[name]
        for plane in fault:
            fplane = np.array(plane)

            #count_tri += 1
            x0 = fplane[0]-corner; B = fplane[1]-corner; C = fplane[2]-corner
            nv = np.cross(B-x0,C-x0); nv /= np.linalg.norm(nv)
            fdata = get_fault_Hydro_data(Data,nv,faults_method)
            nv_abs = np.abs(nv)
            normal_vecs = [B-x0, C-B, x0-C]
            normal_vecs[0] /= np.linalg.norm(normal_vecs[0])
            normal_vecs[1] /= np.linalg.norm(normal_vecs[1])
            normal_vecs[2] /= np.linalg.norm(normal_vecs[2])

            max_vals = np.max(fplane,0) - corner
            min_vals = np.min(fplane,0) - corner

            j_min = int(max(ceil(min_vals[0]/dX),1)); j_max = int(min(ceil(max_vals[0]/dX),cols))            
            i_min = int(max(ceil(min_vals[1]/dY),1)); i_max = int(min(ceil(max_vals[1]/dY),rows))

            for i in range(i_min-1,i_max):
                yc = (i)*dY
                I = rows - (i+1)
                for j in range(j_min-1,j_max):
                    xc = (j)*dX

                    for L in range(layers):
                        if L>0:
                            z_up = Z_bottoms[L-1,I,j]; z_low = Z_bottoms[L,I,j]
                        else:
                            z_up = Z_top[I,j]; z_low = Z_bottoms[0,I,j]

                        if max(z_low,min_vals[2]) < min(z_up,max_vals[2]):

                            center_cell = np.array([xc+dX/2,yc+dY/2,(z_up+z_low)/2.0])
                            D = np.dot(x0-center_cell,nv)
                                                        
                            h_xyz[2] = (z_up-z_low)/2.0
                            r = np.dot(nv_abs,h_xyz)
                            
                            if abs(D)<=r+epsilon:
                                a_cell = x0-center_cell
                                b_cell = B-center_cell
                                c_cell = C-center_cell
                                for p in range(3):
                                    nv_tri = normal_vecs[p]
                                    for q in range(3):                                        
                                        bool_cross = eval_proy(a_cell, b_cell, c_cell,cross_e(q,nv_tri),h_xyz,epsilon)

                                        if bool_cross:
                                            break

                                    if bool_cross:
                                            break

                                if not(bool_cross):
                                    K_ver[L,I,j] = fdata[2]
                                    K_hor[L,I,j] = fdata[0]
                                    if anisotropy_bool:
                                        K_anisotropy_hor[L,I,j] = fdata[1]
                                    if ibound_bool:
                                        I_bound[L,I,j] = fdata[3]
   
def cross_e(k,vec):
    """
    cross product between the canonical vector ek and a given vector, i.e,
    ek x vec.
    
    Args:
        (int) k: integer of the canonical vector, e.g, e1 = (1,0,0)

        (numpy.array) vec: vector.
    """
    if k==0:
        return np.array([0,-vec[2],vec[1]])
    elif k==1:
        return np.array([vec[2],0,-vec[0]])
    else:
        return np.array([-vec[1],vec[0],0])

def eval_proy(a,b,c,vec,h_xyz,epsilon):
    """
    Determines if a triangle and a box have a separeted axis.
    (See Seperating Axis Therom or SAT)
    
    Args:
        (numpy.array) a,b,c: corners of the triangle.

        (numpy.array) vec: vector of the axis.

        (numpy.array) h_xyz: size of the box (hx, hy, hz).

        (float) epsilon: numerical error allowed.

    """
    
    pa = np.dot(a,vec)
    pb = np.dot(b,vec)
    pc = np.dot(c,vec)
    r = np.dot(np.abs(vec),h_xyz)
    b1 = min(pa,pb,pc)>r+epsilon
    b2 = max(pa,pb,pc)< -(r+epsilon)
    if b1 or b2:
        return True
    else:
        return False

def get_fault_Hydro_data(data,nv,f_method):

    """
    Gets the hydraulic parameters of a triangle of a fault depending of f_method.
    
    Args:
        (tuple) data: hydraulic parameters (Kh_x,ani,Kv,i_bound)

        (numpy.array) nv: normal vector of a triangle plane.

        (unicode) f_method: faults method.

    """

    if f_method==u'regular':
        return data
    else:
        eps_angle = 0.98480775301220802 # cos(10)
        Kh = data[0]; Kv = data[2]; ibound = data[3]

        if abs(nv[0])>eps_angle:
            return (Kv,Kh/Kv,Kh,ibound)
        elif abs(nv[1])>eps_angle:
            return (Kh,Kv/Kh,Kh,ibound)
        elif abs(nv[2])>eps_angle:
            return (Kh,1,Kv,ibound)
        else:
            bus_v = np.cross(nv,np.array([nv[1],-nv[0],0.]))
            bus_v /= np.linalg.norm(bus_v)
            data_aprox = Kh*np.abs(bus_v)+Kv*np.abs(nv)
            return (data_aprox[0],data_aprox[1]/data_aprox[0],data_aprox[2],ibound)


def regular_grid(model,rows,cols,layers,Z_top,X_inf,Y_sup,
    dX,dY,bottom_min, units_data):
    """
    Generates a regular finite difference grid to use in MODFLOW.
    
    Args:
        (GeologicalModel) model: geomodlr model to work on.

        (int) rows: number of rows in the MODFLOW model.

        (int) cols: number of cols in the MODFLOW model.

        (int) layers: number of layers in the MODFLOW model.

        (numpy.array) Z_top: topogrphy of the grid model.

        (float) X_inf: x coordinate of upper left corner of the grid.

        (float) Y_sup: y coordinate of upper left corner of the grid.

        (float) dX: spacings along a col.

        (float) dY: spacings along a row.

        (float) bottom_min: model bottom.

        (dic) units_data: units data.
    """

    Z_bottoms = np.zeros((layers,rows,cols))

    I_bound = np.ones((layers, rows, cols), dtype=np.int32)

    K_hor = np.zeros((layers,rows,cols))
    K_anisotropy_hor = np.zeros((layers,rows,cols))
    K_ver = np.zeros((layers,rows,cols))

    for i in np.arange(rows):
        
        yp = Y_sup - (2*i+1)*dY/2.0

        for j in np.arange(cols):

            xp = X_inf + (2*j+1)*dX/2.0
            z_max = Z_top[i,j]

            z_vec = np.linspace(z_max, bottom_min,layers+1)

            Z_bottoms[:,i,j] = z_vec[1:]

    return(Z_bottoms)

def adaptive_grid(model,rows,cols,layers, Z_top,X_inf,Y_sup,
    dX,dY,bottom_min,units_data,angle,dz_min):
    """
    Generates an adaptive finite difference grid to use in MODFLOW.
    
    Args:
        (GeologicalModel) model: geomodlr model to work on.

        (int) rows: number of rows in the MODFLOW model.

        (int) cols: number of cols in the MODFLOW model.

        (int) layers: number of layers in the MODFLOW model.

        (numpy.array) Z_top: topogrphy of the grid model.

        (float) X_inf: x coordinate of upper left corner of the grid.

        (float) Y_sup: y coordinate of upper left corner of the grid.

        (float) dX: spacings along a col.

        (float) dY: spacings along a row.

        (float) bottom_min: model bottom.

        (dic) units_data: units data.

        (float) angle: maximum angle (degrees) between two consecutive grid
        points in the same layer.

        (float) dz_min: minimum allowed distance between layers.
    """

    Z_Bool_Top = np.isfinite(Z_top)

    Z_bottoms = []
    Pos_List = []
    Mat_Order = np.zeros((rows, cols), dtype=np.int16)

    Max_Tan = np.tan(angle*np.pi/180.0)

    Z_Layer_L = np.zeros((rows,cols))

    for i in np.arange(rows):
        yp = Y_sup - (2*i+1)*dY/2.0

        for j in np.arange(cols):

            xp = X_inf + (2*j+1)*dX/2.0
            z_max = Z_top[i,j]
            dz = (z_max - bottom_min)/layers            
            zl = z_max-dz

            #z_mean,change = find_unit_limits(model, xp, yp, z_max, zl, 1E-3)
            z_mean,change = find_unit_limits_CPP(model, xp, yp, z_max, zl, 1E-3)
            Z_Bool_Top[i,j] = change

            if change:
                Pos_List.append((i,j))

            Z_Layer_L[i,j] = min(z_mean - z_max,-dz_min) + z_max


    Z_bottoms.append(Z_Layer_L.copy())

    Layer_Correction(Pos_List,Mat_Order,Z_Bool_Top,Z_top,Z_bottoms,0,
        Max_Tan, rows,cols,dX,dY,dz_min)
    
    Mat_Order*=0

    Count_Lay = 0
    Z_Bool_Bot = np.isfinite(Z_top)     

    for L in np.arange(1,layers-1):

        for i in np.arange(rows):
            yp = Y_sup - (2*i+1)*dY/2.00

            for j in np.arange(cols):

                xp = X_inf + (2*j+1)*dX/2.0                
                zp = Z_bottoms[Count_Lay][i,j]
                #zl = zp - (zp- bottom_min)/(layers-L)
                zl = (Z_top[i,j]*(layers-(L+1)) + bottom_min*(L+1))/layers
                
                #z_mean,change = find_unit_limits(model, xp, yp, zp, zl,1E-3)
                z_mean,change = find_unit_limits_CPP(model, xp, yp, zp, zl, 1E-3)

                if (change) and (not(Z_Bool_Top[i,j])):
                    Pos_List.append((i,j))

                Z_Bool_Bot[i,j] = change
                Z_Layer_L[i,j] = min(z_mean - zp,-dz_min) + zp


        #if (np.sum(Z_Bool_Top & Z_Bool_Bot) == 0):
        if not(np.any(Z_Bool_Top & Z_Bool_Bot)):

            Z_bottoms[Count_Lay][Z_Bool_Bot | np.logical_not(Z_Bool_Top)] = Z_Layer_L[Z_Bool_Bot | np.logical_not(Z_Bool_Top)]
            Z_Bool_Top = Z_Bool_Top | Z_Bool_Bot

            Layer_Correction(Pos_List,Mat_Order,Z_Bool_Top,Z_top,Z_bottoms,
                Count_Lay, Max_Tan, rows,cols,dX,dY,dz_min)

            Mat_Order*=0

        else:

            Count_Lay += 1
            Z_Bool_Top = Z_Bool_Bot.copy()
            Pos_List = np.where(Z_Bool_Top)
            Pos_List = zip(Pos_List[0],Pos_List[1])
            Z_bottoms.append(Z_Layer_L.copy())

        Mat_Order*=0

    Z_bottoms.append(bottom_min*np.ones((rows,cols)))
    layers = len(Z_bottoms)
    Z_bottoms= np.array(Z_bottoms)

    return((Z_bottoms,layers))

def set_unit_properties(model,units_data,Z_top,Z_bottoms,rows,cols,layers,X_inf,Y_sup,
    dX,dY,faults_data):
    """
    Generates the hidraulic conductivy and I_bound (active or not active) matrices.
    
    Args:
        (GeologicalModel) model: geomodlr model to work on.

        (dic) units_data: units data.

        (numpy.array) Z_top: topogrphy of the grid model.

        (numpy.array) Z_bottoms: bottoms of the layers.

        (int) rows: number of rows in the MODFLOW model.

        (int) cols: number of cols in the MODFLOW model.

        (int) layers: number of layers in the MODFLOW model.        

        (float) X_inf: x coordinate of upper left corner of the grid.

        (float) Y_sup: y coordinate of upper left corner of the grid.
        
    """

    aux_data=np.array(units_data.values() + faults_data.values())
    anisotropy_bool = np.max(aux_data[:,1])==np.min(aux_data[:,1])
    ibound_bool = np.max(aux_data[:,3])==np.min(aux_data[:,3])

    if anisotropy_bool:
        K_anisotropy_hor = aux_data[0,1]
        chani_var = aux_data[0,1]
    else:
        K_anisotropy_hor = np.zeros((layers,rows,cols))
        chani_var = -1

    if ibound_bool:
        I_bound = aux_data[0,3]
    else:
        I_bound = np.ones((layers, rows, cols), dtype=np.int32)

    K_hor = np.zeros((layers,rows,cols))
    K_ver = np.zeros((layers,rows,cols))

    anisotropy_bool = not(anisotropy_bool)
    ibound_bool = not(ibound_bool)
    
    for L in np.arange(0,layers):

        if (L>0):
            mid_points = ( Z_bottoms[L-1,:,:] + Z_bottoms[L,:,:])/2.0
        else:
            mid_points = ( Z_top + Z_bottoms[L,:,:])/2.0
        
        for i in np.arange(rows):
            yp = Y_sup - (2*i+1)*dY/2.0

            for j in np.arange(cols):
    
                xp = X_inf + (2*j+1)*dX/2.0

                Unit = model.closest([xp,yp,mid_points[i,j]])[0]
                Data = units_data[Unit]
                K_hor[L,i,j] = Data[0];
                if anisotropy_bool:
                    K_anisotropy_hor[L,i,j] = Data[1]
                K_ver[L,i,j] = Data[2]
                if ibound_bool:
                    I_bound[L,i,j] = Data[3]

    return((K_hor, K_anisotropy_hor, K_ver, I_bound,chani_var))
    
def find_unit_limits_CPP(model, xp, yp, z_max, z_min, eps):
    return model.find_unit_limits(xp, yp, z_max, z_min, eps)    

def find_unit_limits(model, xp, yp, z_max, z_min, eps):

    #return model.find_unit_limits(xp, yp, z_max, z_min, eps)
    """
    Given two points with the same x and y coordinates, and different z
    coordinates, it determines if exist a unit change between these points.
    If this happens, it return the z coordinate of this change. This algorithm
    is based on the Bisection method.
    
    Args:
        (GeologicalModel) model: geomodlr model to work on.

        (float) xp: X coordinate of the two points.

        (float) yp: Y coordinate of the two points.

        (float) z_max: z coordinate of the upper point.

        (float) z_min: z coordinate of the lower point.

        (float) eps: absolute error to converge.
    """

    Unit_max = model.closest([xp,yp,z_max])[0]
    Unit_min = model.closest([xp,yp,z_min])[0]

    z_mean = (z_max+z_min)/2.0
    Unit_mean = model.closest([xp, yp, z_mean])[0]

    if (Unit_max == Unit_mean) and (Unit_min == Unit_mean):
        change = False
    else:
        change = True
        while (z_max-z_min)>eps:

            if (Unit_max == Unit_mean):
                z_max = z_mean
            elif (Unit_min == Unit_mean):
                z_min = z_mean
            else:
                z_min = z_mean
                Unit_min = Unit_mean

            z_mean = (z_max+z_min)/2.0
            Unit_mean = model.closest([xp, yp, z_mean])[0]

    return((z_min, change))

def Layer_Correction(Pos_List,Mat_Order,Z_Bool_Top,Z_top,Z_Bottoms,Layer,
    Max_Tan,Rows,Cols,dX,dY,dz_min):
    """
    Corrects the layers points to guarantee a given maximum angle between
    two consecutive points.
    
    Args:
        (list) Pos_List: list of tuples with the row-col indexes of the points,
        which have been fixed in the layer.

        (numpy.array) Mat_Order: integer matrix that shows in which order
        the points of the layer have been fixed. 

        (numpy.array Z_Bool_Top: Boolean matrix that shows which points has
        been fixed in the above layer.

        (numpy.array) Z_top: topogrphy of the grid model.

        (list) Z_Bottoms: bottoms of the layers.

        (int) Layer: layer

        (float) Max_Tan: numpy.tan(angle)

        (int) Rows: number of rows in the MODFLOW model.

        (int) Cols: number of cols in the MODFLOW model.

        (float) dX: spacings along a col.

        (float) dY: spacings along a row.

        (float) dz_min: minimum allowed distance between layers.
    """
    
    c=0

    while c<len(Pos_List):

        i,j = Pos_List[c]

        # Down
        if i>0:
            I = i-1

            if not(Z_Bool_Top[I,j]) or (Mat_Order[i,j]<Mat_Order[I,j]):
                
                dz=Z_Bottoms[Layer][i,j]-Z_Bottoms[Layer][I,j]

                Change,dz = Angle(0.,dY,dz,Max_Tan)
                if Change:

                    if Layer==0:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[I,j]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer-1][I,j]-dz_min)

                    Z_Bottoms[Layer][I,j] = Val
                    Mat_Order[I,j] = Mat_Order[i,j] + 1
                    Z_Bool_Top[I,j] = Change
                    Pos_List.append((I,j))

        # Up
        if i< (Rows-1):
            I = i +1

            if not(Z_Bool_Top[I,j]) or (Mat_Order[i,j]<Mat_Order[I,j]):

                dz=Z_Bottoms[Layer][i,j]-Z_Bottoms[Layer][I,j]

                Change,dz = Angle(0.,dY,dz,Max_Tan)
                if Change:
                    
                    if Layer==0:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[I,j]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer-1][I,j]-dz_min)
                    
                    Z_Bottoms[Layer][I,j] = Val
                    Mat_Order[I,j] = Mat_Order[i,j] + 1
                    Z_Bool_Top[I,j] = Change
                    Pos_List.append((I,j))

        # Left
        if j>0:
            J = j-1

            if not(Z_Bool_Top[i,J]) or (Mat_Order[i,j]<Mat_Order[i,J]):

                dz=Z_Bottoms[Layer][i,j]-Z_Bottoms[Layer][i,J]

                Change,dz = Angle(dX,0.,dz,Max_Tan)
                if Change:

                    if Layer==0:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[i,J]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer-1][i,J]-dz_min)

                    Z_Bottoms[Layer][i,J] = Val
                    Mat_Order[i,J] = Mat_Order[i,j] + 1
                    Z_Bool_Top[i,J] = Change
                    Pos_List.append((i,J))

        # Right
        if j<Cols-1:
            J = j+1

            if not(Z_Bool_Top[i,J]) or (Mat_Order[i,j]<Mat_Order[i,J]):
                
                dz=Z_Bottoms[Layer][i,j]-Z_Bottoms[Layer][i,J]

                Change,dz = Angle(dX,0.,dz,Max_Tan)                
                if Change:

                    if Layer==0:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[i,J]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer-1][i,J]-dz_min)

                    Z_Bottoms[Layer][i,J] = Val
                    Mat_Order[i,J] = Mat_Order[i,j] + 1
                    Z_Bool_Top[i,J] = Change
                    Pos_List.append((i,J))

        c += 1

def Angle(dx,dy,dz,Max_Tan):
    """
    Determines the angle between two points is lower or greater than the
    given max angle.
    
    Args:
        (float) dx: distance in the X direction between the two points.

        (float) dy: distance in the Y direction between the two points.

        (float) dz: distance in the Z direction between the two points.

        (float) Max_Tan: numpy.tan(angle)
    """

    Tan_O = dz/np.sqrt(dx**2+dy**2)

    if (Tan_O>Max_Tan):
        dz = Max_Tan*np.sqrt(dx**2+dy**2)
        Change = True
    else:
        Change = False

    return((Change,dz))
