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
import os
from sys import exit
from math import ceil, floor
from scipy import interpolate

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

def create_modflow_inputs(model_name, model, Units_data,
    length_units=LENGTH_UNIT.SECONDS, N_row=100, N_col=100, N_layers=100,
    Bbox=None, angle=20, dz_min = 1.0, time_units=TIME_UNIT.METERS, 
    Class=ALGORITHM.REGULAR):
    """
    Generates the DIS, BAS, LPF and NAM files, which are used by classical
    MODFLOW processors. The user has to import the NAM file from his MODFLOW
    program. 
    
    Args:
        (string) model_name: name of generated files.

        (GeologicalModel) model: geomodlr model to work on.

        (dict) Units_data: the dictionary keys correpond to the unit names
        and the items are tuples with 4 values: (Kx,rate,Kz,i_bound).
        Kx and Kz are the hydraulic conductivy in the respective axis and
        rate variable is the rate between Ky and Kx, i.e., rate=Ky/Kx.
        Finally, i_bound is equal to 1 if the respective unit is active and 0
        otherwise.

        (int) length_units: length units (see: class LENGTH_UNIT)

        (int) N_row: number of rows in the MODFLOW model.

        (int) N_col: number of cols in the MODFLOW model.

        (int) N_layers: number of layers in the MODFLOW model.

        (list) Bbox: the bounding box to search in.

        (float) angle: this angle (degrees) is used in the adative algorithm.

        (float) dz_min: minimum allowed distance between layers.

        (int) time_units: time units (see: class TIME_UNIT)

        (string) Class: grid generation algorithm (see: class ALGORITHM)

    """

    try:
        num_unit = len(model.units)
    except:
        exit('You have not defined the Geomodelr model')

    if (Bbox is None):
        Bbox = model.bbox

    X_inf = Bbox[0]
    X_sup = Bbox[3]

    Y_inf = Bbox[1]
    Y_sup = Bbox[4]

    dX = (X_sup - X_inf)/N_col
    dY = (Y_sup - Y_inf)/N_row

    #X_vec = np.linspace(X_inf + dX/2., X_sup - dX/2.,N_col)
    #Y_vec = np.linspace(Y_inf + dY/2., Y_sup - dY/2.,N_row)
 
    Z_top = np.zeros((N_row,N_col))

    Bottom_min = Bbox[2] + 1E-5

    # Define Z-top
    for i in np.arange(N_row):
        for j in np.arange(N_col):

            xp = X_inf + (2*j+1)*dX/2.0
            yp = Y_inf + (2*i+1)*dY/2.0
            Z_top[i,j] = model.height((xp, yp))


    Z_top_min = np.min(Z_top)

    if ((Z_top_min-Bottom_min)/N_layers < dz_min):
        N_layers = int((Z_top_min-Bottom_min)/dz_min)
        print('Maximum number of initial layers: ' + str(N_layers))

    
    # ---------- Define Z-top, Z-bottoms and Hydraulic conductivity
    # Class Fine
    if (Class is 'regular'):

        Z_bottoms,K_hor,K_ratio_hor,K_ver,I_bound=regular_grid(model,
            N_row,N_col,N_layers,Z_top,X_inf,Y_inf,dX,dY,Bottom_min,
            Units_data)
    
    elif (Class is 'adaptive'):

        Z_bottoms,K_hor,K_ratio_hor,K_ver,I_bound,N_layers=adaptive_grid(model,
            N_row,N_col,N_layers,Z_top,X_inf,Y_inf,dX,dY,Bottom_min,
            Units_data,angle,dz_min)


    #  ------- Flowpy Package ----
    # Grid

    chani_var = -np.ones(N_layers,dtype=np.int32)

    mf_handle = fp.modflow.Modflow(model_name, exe_name='mf2005',verbose=False)
    # Variables for the Dis package
    dis = fp.modflow.ModflowDis(mf_handle,nlay=N_layers, nrow=N_row, ncol=N_col,
        top=Z_top, botm=Z_bottoms, delc=dY, delr=dX, xul=X_inf, yul=Y_sup,
        itmuni=time_units, lenuni=length_units)

    # Variables for the BAS package
    bas = fp.modflow.ModflowBas(mf_handle,ibound = I_bound)

    # Add LPF package to the MODFLOW model
    lpf = fp.modflow.ModflowLpf(mf_handle, chani=chani_var, hk=K_hor,
        vka=K_ver, hani=K_ratio_hor,laytyp=np.ones(N_layers,dtype=np.int32))

    mf_handle.write_input()
    
    mv_comand = 'mv ' + model_name + '* /media/sf_CompartidaVB/'
    os.system(mv_comand)

# ===================== AUXILIAR FUNCTIONS ========================

def regular_grid(model,N_row,N_col,N_layers,Z_top,X_inf,Y_inf,
    dX,dY,Bottom_min, Units_data):
    """
    Generates a regular grid to use in MODFLOW.
    
    Args:
        (GeologicalModel) model: geomodlr model to work on.

        (int) N_row: number of rows in the MODFLOW model.

        (int) N_col: number of cols in the MODFLOW model.

        (int) N_layers: number of layers in the MODFLOW model.

        (numpy.array) Z_top: topogrphy of the grid model.

        (float) X_inf: origin X coordinate.

        (float) Y_inf: origin Y coordinate.

        (float) dX: spacings along a col.

        (float) dY: spacings along a row.

        (float) Bottom_min: model bottom.

        (dic) Units_data: units data.
    """


    Z_bottoms = np.zeros((N_layers,N_row,N_col))

    I_bound = np.ones((N_layers, N_row, N_col), dtype=np.int32)

    K_hor = np.zeros((N_layers,N_row,N_col))
    K_ratio_hor = np.zeros((N_layers,N_row,N_col))
    K_ver = np.zeros((N_layers,N_row,N_col))

    for i in np.arange(N_row):
        
        yp = Y_inf + (2*i+1)*dY/2.0

        for j in np.arange(N_col):

            xp = X_inf + (2*j+1)*dX/2.0
            z_max = Z_top[i,j]

            z_vec = np.linspace(z_max, Bottom_min,N_layers+1)

            Z_bottoms[:,i,j] = z_vec[1:]

            for L in np.arange(N_layers):
                zp = (z_vec[L] + z_vec[L+1])/2.0

                Unit = model.closest([xp,yp,zp])[0]
                Data = Units_data[Unit]
                K_hor[L,i,j] = Data[0]; K_ratio_hor[L,i,j] = Data[1]
                K_ver[L,i,j] = Data[2]
                I_bound[L,i,j] = Data[3]

    return((Z_bottoms,K_hor, K_ratio_hor, K_ver, I_bound))


def adaptive_grid(model,N_row,N_col,N_layers, Z_top,X_inf,Y_inf,
    dX,dY,Bottom_min,Units_data,angle,dz_min):
    """
    Generates an adaptive grid to use in MODFLOW.
    
    Args:
        (GeologicalModel) model: geomodlr model to work on.

        (int) N_row: number of rows in the MODFLOW model.

        (int) N_col: number of cols in the MODFLOW model.

        (int) N_layers: number of layers in the MODFLOW model.

        (numpy.array) Z_top: topogrphy of the grid model.

        (float) X_inf: origin X coordinate.

        (float) Y_inf: origin Y coordinate.

        (float) dX: spacings along a col.

        (float) dY: spacings along a row.

        (float) Bottom_min: model bottom.

        (dic) Units_data: units data.

        (float) angle: maximum angle (degrees) between two consecutive grid
        points in the same layer.

        (float) dz_min: minimum allowed distance between layers.
    """

    Z_Bool_Top = np.isfinite(Z_top)

    Z_bottoms = []
    Pos_List = []
    Mat_Order = np.zeros((N_row, N_col), dtype=np.int16)

    Max_Tan = np.tan(angle*np.pi/180.0)

    Z_Layer_L = np.zeros((N_row,N_col))

    for i in np.arange(N_row):
        yp = Y_inf + (2*i+1)*dY/2.0

        for j in np.arange(N_col):

            xp = X_inf + (2*j+1)*dX/2.0
            z_max = Z_top[i,j]

            dz = (z_max - Bottom_min)/N_layers
            
            zl = z_max-dz

            z_mean,change = find_unit_limits(model, xp, yp, z_max,
                zl, 1E-3)

            Z_Bool_Top[i,j] = change

            if change:
                Pos_List.append((i,j))

            Z_Layer_L[i,j] = min(z_mean - z_max,-dz_min) + z_max


    Z_bottoms.append(Z_Layer_L.copy())

    Layer_Correction(Pos_List,Mat_Order,Z_Bool_Top,Z_top,Z_bottoms,0,True,
        Max_Tan, N_row,N_col,dX,dY,dz_min)
    
    Mat_Order*=0

    Count_Lay = 0
    Z_Bool_Bot = np.isfinite(Z_top)     

    for L in np.arange(1,N_layers-1):

        for i in np.arange(N_row):
            yp = Y_inf + (2*i+1)*dY/2.0

            for j in np.arange(N_col):

                xp = X_inf + (2*j+1)*dX/2.0
                
                zp = Z_bottoms[Count_Lay][i,j]
                #zl = zp - (zp- Bottom_min)/(N_layers-L)
                zl = (Z_top[i,j]*(N_layers-(L+1)) + Bottom_min*(L+1))/N_layers
                

                z_mean,change = find_unit_limits(model, xp, yp, zp, zl,
                    1E-3)

                if (change) & (not(Z_Bool_Top[i,j])):
                    Pos_List.append((i,j))

                Z_Bool_Bot[i,j] = change
                #Z_bottoms[L,i,j] = min(z_mean - zp,-1.01) + zp
                Z_Layer_L[i,j] = min(z_mean - zp,-dz_min) + zp


        if (np.sum(Z_Bool_Top & Z_Bool_Bot) == 0):

            #Layers_Bool[L-1] = False
            #Z_bottoms[L,Z_Bool_Top] = Z_bottoms[L-1,Z_Bool_Top]
            #Z_bottoms[Count_Lay][Z_Bool_Bot] = Z_Layer_L[Z_Bool_Bot]
            Z_bottoms[Count_Lay][Z_Bool_Bot | np.logical_not(Z_Bool_Top)] = Z_Layer_L[Z_Bool_Bot | np.logical_not(Z_Bool_Top)]
            Z_Bool_Top = Z_Bool_Top | Z_Bool_Bot
            #Z_bottoms[Count_Lay][np.logical_not(Z_Bool_Top)] = Z_Layer_L[np.logical_not(Z_Bool_Top)]

            Layer_Correction(Pos_List,Mat_Order,Z_Bool_Top,Z_top,Z_bottoms,
                Count_Lay, False, Max_Tan, N_row,N_col,dX,dY,dz_min)

            Mat_Order*=0


        else:
            # Smoother

            #Z_bottoms[L-1,Aux_Bool] = Z_bottoms[Last_Layer,Aux_Bool]

            Count_Lay += 1
            Z_Bool_Top = Z_Bool_Bot.copy()
            Pos_List = np.where(Z_Bool_Top)
            Pos_List = zip(Pos_List[0],Pos_List[1])
            Z_bottoms.append(Z_Layer_L.copy())

        Mat_Order*=0


    Z_bottoms.append(Bottom_min*np.ones((N_row,N_col)))
    N_layers = len(Z_bottoms)

    Z_bottoms= np.array(Z_bottoms)

    K_hor = np.zeros((N_layers,N_row,N_col))
    K_ratio_hor = np.zeros((N_layers,N_row,N_col))
    K_ver = np.zeros((N_layers,N_row,N_col))

    I_bound = np.ones((N_layers, N_row, N_col), dtype=np.int32)
    for L in np.arange(0,N_layers):

        if (L>0):
            mid_points = ( Z_bottoms[L-1,:,:] + Z_bottoms[L,:,:])/2.0
        else:
            mid_points = ( Z_top + Z_bottoms[L,:,:])/2.0
        
        for i in np.arange(N_row):
            yp = Y_inf + (2*i+1)*dY/2.0

            for j in np.arange(N_col):
    
                xp = X_inf + (2*j+1)*dX/2.0

                Unit = model.closest([xp,yp,mid_points[i,j]])[0]
                Data = Units_data[Unit]
                K_hor[L,i,j] = Data[0]; K_ratio_hor[L,i,j] = Data[1]
                K_ver[L,i,j] = Data[2]
                I_bound[L,i,j] = Data[3]

    return((Z_bottoms,K_hor, K_ratio_hor, K_ver, I_bound,
        N_layers))

# find_unit_limits: It finds the limit between two geological units.
def find_unit_limits(model, xp, yp, z_max, z_min, eps):
    """
    Given two points with the same x and y coordinates, and different z
    coordinates, determines if exist a unit change between these points.
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
        #z_mean = z_min
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



def Layer_Correction(Pos_List,Mat_Order,Z_Bool_Top,Z_top,Z_Bottoms,Layer,Code,
    Max_Tan,Rows,Cols,dX,dY,dz_min):
    """
    Corrects the layers points to guarantee a given maximum angle between
    two consecutive points.
    
    Args:
        (list) Pos_List: list of tuples with the row-col indexes of the points,
        which have been fixed in the layer.

        (numpy.array Mat_Order: integer matrix that shows in which order
        the points of the layer have been fixed. 

        (numpy.array Z_Bool_Top: Boolean matrix that shows which points has
        been fixed in the above layer.

        (numpy.array) Z_top: topogrphy of the grid model.

        (list) Z_Bottoms: bottoms of the layers.

        (int) Layer: layer

        (bool) Code: If it is equal True, it is working at the first layer,
        and False otherwise.

        (float) Max_Tan: numpy.tan(angle)

        (float) dY: spacings along a row.

        (float) Bottom_min: model bottom.

        (dic) Units_data: units data.

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

                    if Code:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[I,j]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer][I,j]-dz_min)

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
                    
                    if Code:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[I,j]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer][I,j]-dz_min)
                    
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

                    if Code:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[i,J]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer][i,J]-dz_min)

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

                    if Code:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_top[i,J]-dz_min)
                    else:
                        Val = min(Z_Bottoms[Layer][i,j]-dz,Z_Bottoms[Layer][i,J]-dz_min)

                    Z_Bottoms[Layer][i,J] = Val
                    Mat_Order[i,J] = Mat_Order[i,j] + 1
                    Z_Bool_Top[i,J] = Change
                    Pos_List.append((i,J))

        c += 1


def Angle(dx,dy,dz,Max_Tan):
    """
    Determines if two points has an angle lower or greater than the
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
