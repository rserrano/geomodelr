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

class TIME_UNIT:
    UNDEFINED = 0
    FEET =1
    METERS = 2
    CENTIMETERS = 3

class LENGTH_UNIT:
    UNDEFINED = 0
    SECONDS = 1
    MINUTES = 2
    HOURS = 3
    DAYS = 4
    YEARS = 5

class ALGORITHM:
    REGULAR = 'regular'
    ADAPTIVE = 'adaptive'

def create_modflow_inputs(name, model, units_data,
    length_units=LENGTH_UNIT.SECONDS, rows=100, cols=100, layers=100,
    bbox=None, angle=20, dz_min = 1.0, time_units=TIME_UNIT.METERS, 
    algorithm=ALGORITHM.REGULAR):
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

    # Define Z-top
    for i in np.arange(rows):
        yp = Y_sup - (2*i+1)*dY/2.0
        for j in np.arange(cols):
            xp = X_inf + (2*j+1)*dX/2.0
            Z_top[i,j] = model.height((xp, yp))


    Z_top_min = np.min(Z_top)

    if ((Z_top_min-bottom_min)/layers < dz_min):
        layers = int((Z_top_min-bottom_min)/dz_min)
    
    if (algorithm is 'regular'):

        Z_bottoms=regular_grid(model, rows,cols,layers,Z_top,X_inf,Y_sup,dX,dY,
            bottom_min,units_data)
    
    elif (algorithm is 'adaptive'):

        Z_bottoms,layers=adaptive_grid(model,rows,cols,layers,Z_top,X_inf,Y_sup,dX,dY,
            bottom_min,units_data,angle,dz_min)

    K_hor, K_anisotropy_hor, K_ver, I_bound,chani_var=set_unit_properties(model,
        units_data,Z_top,Z_bottoms,rows,cols,layers,X_inf,Y_sup,dX,dY)

    #  ------- Flowpy Packages ----
    # Grid

    #geo=fp.utils.reference.SpatialReference(delr=dX*np.ones(cols),delc=dY*np.ones(rows),
        #lenuni=length_units, xll=X_inf, yll=Y_inf,units='meters',epsg=3116)

    mf_handle = fp.modflow.mf.Modflow(modelname=name,namefile_ext='nam',
        version='mf2005')
    # Variables for the Dis package
    dis = fp.modflow.ModflowDis(mf_handle,nlay=layers, nrow=rows, ncol=cols,
        top=Z_top, botm=Z_bottoms, delc=dY, delr=dX, xul=X_inf, yul=Y_sup,
        itmuni=time_units, lenuni=length_units)

    #dis.sr=geo

    # Variables for the BAS package
    bas = fp.modflow.ModflowBas(mf_handle,ibound = I_bound)

    # Add LPF package to the MODFLOW model
    lpf = fp.modflow.ModflowLpf(mf_handle, chani=chani_var, hk=K_hor,
        vka=K_ver, hani=K_anisotropy_hor,laytyp=np.ones(layers,dtype=np.int32))#

    mf_handle.write_input()
#
    #return((mf_handle,dis,geo))

# ===================== AUXILIAR FUNCTIONS ========================

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

            z_mean,change = find_unit_limits(model, xp, yp, z_max,
                zl, 1E-3)

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
                

                z_mean,change = find_unit_limits(model, xp, yp, zp, zl,
                    1E-3)

                if (change) & (not(Z_Bool_Top[i,j])):
                    Pos_List.append((i,j))

                Z_Bool_Bot[i,j] = change
                Z_Layer_L[i,j] = min(z_mean - zp,-dz_min) + zp


        if (np.sum(Z_Bool_Top & Z_Bool_Bot) == 0):

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
    dX,dY):
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

    aux_data=np.array(units_data.values())
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
                if not(anisotropy_bool):
                    K_anisotropy_hor[L,i,j] = Data[1]
                K_ver[L,i,j] = Data[2]
                if not(ibound_bool):
                    I_bound[L,i,j] = Data[3]

    return((K_hor, K_anisotropy_hor, K_ver, I_bound,chani_var))

def find_unit_limits(model, xp, yp, z_max, z_min, eps):
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

    if (Unit_max == Unit_mean) & (Unit_min == Unit_mean):
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
                Unit_min == Unit_mean

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

            if not(Z_Bool_Top[I,j]) | (Mat_Order[i,j]<Mat_Order[I,j]):
                
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

            if not(Z_Bool_Top[I,j]) | (Mat_Order[i,j]<Mat_Order[I,j]):

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

            if not(Z_Bool_Top[i,J]) | (Mat_Order[i,j]<Mat_Order[i,J]):

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

            if not(Z_Bool_Top[i,J]) | (Mat_Order[i,j]<Mat_Order[i,J]):
                
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
