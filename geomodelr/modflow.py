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

# Creates the data from Geomodelr to ModFlow
def create_modflow_inputs( model_name='no_name', model=None, N_row=10, N_col=10,
	N_layers=5, properties=None, Class=1):
	
	try:
		num_unit = len(model.units)
	except:
		exit('You have not defined the Geomodelr model')

	try:
		num_properties = (np.shape(properties))[1]
		if (not( num_unit==num_properties )):
			exit('The Geomodelr model has ' + str(num_unit) +
				'units. Your properties array has ' + str(num_propierties) +
				'units.')
	except:
		exit('You have not defined the properties array.')


	del(num_properties)

	#properties[1,:] /= properties[0,:]

	K_data = { (model.units)[k]: properties[:,k] for k in range(num_unit) }

	X_inf = model.bbox[0]
	X_sup = model.bbox[3]

	Y_inf = model.bbox[1]
	Y_sup = model.bbox[4]

	dX = (X_sup - X_inf)/N_col
	dY = (Y_sup - Y_inf)/N_row

	X_vec = np.linspace(X_inf + dX/2., X_sup - dX/2.,N_col)
	Y_vec = np.linspace(Y_inf + dY/2., Y_sup - dY/2.,N_row)
 
 	Z_top = np.zeros((N_row,N_col))
 	Z_bottoms = np.zeros((N_layers,N_row,N_col))

	Bottom_min = model.bbox[2]


 	# ---------- Define Z-top, Z-bottoms and Hydraulic conductivity
 	# Class Fine
	if (Class == 1):

		chani_var = -np.ones(N_layers, dtype=np.int32)

		K_hor = np.zeros((N_layers,N_row,N_col))
		K_ratio_hor = np.zeros((N_layers,N_row,N_col))
		K_ver = np.zeros((N_layers,N_row,N_col))

		for i in np.arange(N_row):
			for j in np.arange(N_col):

				xp = X_vec[j]
				yp = Y_vec[i]
				z_max = model.height([xp, yp])
				Z_top[i,j] = z_max



				#if (j>0):
				#	z_vec[0] = z_max
				#elif (i > 0):
				#	z_vec = np.insert( Z_bottoms[:,i-1,0] ,0,z_max)
				#else:
				#	z_vec = np.linspace(z_max,model.bbox[2],N_layers+1)

				z_vec = np.linspace(z_max, Bottom_min,N_layers+1)
				define_bottoms(model, xp, yp, N_layers,	z_vec)

				Z_bottoms[:,i,j] = z_vec[1:]

				for L in np.arange(N_layers):
					zp = (z_vec[L] + z_vec[L+1])/2.0

					Unit = model.closest([xp,yp,zp])[0]
					Data = K_data[Unit]
					K_hor[L,i,j] = Data[0]; K_ratio_hor[L,i,j] = Data[1]
					K_ver[L,i,j] = Data[2]
	
	else:

	# ---------- Define Z-top, Z-bottoms and Hydraulic conductivity
 	# Class Smart

		Z_Bool_Top = np.isfinite(Z_top)

		for i in np.arange(N_row):
			for j in np.arange(N_col):

				xp = X_vec[j]
				yp = Y_vec[i]
				z_max = model.height([xp, yp])
				Z_top[i,j] = z_max

				dz = (z_max - Bottom_min)/N_layers

				z_mean,change = find_units_limit(model, xp, yp, z_max,
					z_max-dz, 1E-3)

				Z_Bool_Top[i,j] = change

				Z_bottoms[0,i,j] = min(z_mean - z_max,-1.01) + z_max


		Z_Bool_Bot = np.isfinite(Z_top)
		Layers_Bool = np.isfinite(np.arange(N_layers))

		for L in np.arange(1,N_layers):

			for i in np.arange(N_row):
				for j in np.arange(N_col):

					xp = X_vec[j]
					yp = Y_vec[i]
					zp = Z_bottoms[L-1,i,j]
					
					dz = (zp- Bottom_min)/(N_layers-L)

					z_mean,change = find_units_limit(model, xp, yp, zp,
					zp-dz, 1E-3)

					#print L, '\t', i, '\t', j, '\t',  change, '\n'
					Z_Bool_Bot[i,j] = change
					Z_bottoms[L,i,j] = min(z_mean - zp,-1.01) + zp

			if (np.sum(Z_Bool_Top & Z_Bool_Bot) == 0):

				Layers_Bool[L-1] = False
				Z_bottoms[L,Z_Bool_Top] = Z_bottoms[L-1,Z_Bool_Top]
				Z_Bool_Top = Z_Bool_Top | Z_Bool_Bot

			else:

				Layers_Bool[L-1] = True
				Z_Bool_Top = Z_Bool_Bot.copy()


		N_layers = np.sum(Layers_Bool)
		Z_bottoms = Z_bottoms[Layers_Bool,:,:]
		Z_bottoms[-1,:,:] = Bottom_min

		K_hor = np.zeros((N_layers,N_row,N_col))
		K_ratio_hor = np.zeros((N_layers,N_row,N_col))
		K_ver = np.zeros((N_layers,N_row,N_col))

		chani_var = -np.ones(N_layers, dtype=np.int32)

		for L in np.arange(0,N_layers):

			if (L>0):
				mid_points = ( Z_bottoms[L-1,:,:] + Z_bottoms[L,:,:])/2.0
			else:
				mid_points = ( Z_top + Z_bottoms[L,:,:])/2.0
			
			for i in np.arange(N_row):
				for j in np.arange(N_col):
		
					xp = X_vec[j]
					yp = Y_vec[i]

					Unit = model.closest([xp,yp,mid_points[i,j]])[0]
					Data = K_data[Unit]
					K_hor[L,i,j] = Data[0]; K_ratio_hor[L,i,j] = Data[1]
					K_ver[L,i,j] = Data[2]


	#return((N_layers,Z_bottoms,K_hor))

	#  ------- Flopy Package ----
	# Grid

	mf_handle = fp.modflow.Modflow(model_name, exe_name='mf2005',verbose=False)
	dis = fp.modflow.ModflowDis(mf_handle,nlay=N_layers, nrow=N_row, ncol=N_col,
		top=Z_top, botm=Z_bottoms, delc=dY, delr=dX)


	# Variables for the BAS package
	# ibound = np.ones((N_layers, N_row, N_col), dtype=np.int32)
	# ibound[:, :, 0] = -1
	# ibound[:, :, 1] = -2
	# ibound[:, :, 2] = -3
	# strt = np.ones((N_layers, N_row, N_col), dtype=np.float32)
	# strt[:, :, 0] = 10.5
	# strt[:, :, 1] = 5.5
	# strt[:, :, 2] = 1.5
	bas = fp.modflow.ModflowBas(mf_handle)

	# Add LPF package to the MODFLOW model
	lpf = fp.modflow.ModflowLpf(mf_handle, chani=chani_var, hk=K_hor,
		vka=K_ver, hani=K_ratio_hor)
	
	# oc = fp.modflow.ModflowOc(mf_handle)
	# pcg = fp.modflow.ModflowPcg(mf_handle)

	mf_handle.write_input()
	
	mv_comand = 'mv ' + model_name + '* /media/sf_CompartidaVB/'
	os.system(mv_comand)

	return(dis)

# ===================== AUXILIAR FUNCTIONS ========================

# FIND_UNITS_LIMIT: It finds the limit between two geological units.

def find_units_limit(model, xp, yp, Z_top, z_min, eps):

		Unit_max = model.closest([xp,yp,Z_top])[0]
		Unit_min = model.closest([xp,yp,z_min])[0]

		z_mean = (Z_top+z_min)/2.0
		Unit_mean = model.closest([xp, yp, z_mean])[0]

		if (Unit_max == Unit_mean) & (Unit_min == Unit_mean):
			change = False
			#z_mean = z_min
		else:
			change = True
			while (Z_top-z_min)>eps:

				if (Unit_max == Unit_mean):
					Z_top = z_mean
				elif (Unit_min == Unit_mean):
					z_min = z_mean
				else:
					z_min = z_mean
					Unit_min == Unit_mean

				z_mean = (Z_top+z_min)/2.0
				Unit_mean = model.closest([xp, yp, z_mean])[0]

		return((z_min, change))


# DEFINE_BOOTOMS: It finds the layers bottoms
def define_bottoms(model, xp, yp, N_layers,	Z_vec):
	
	K_index = 0
	counter = 0

	#Aux_Bool = np.isfinite(Z_vec)

	for k in np.arange(N_layers-1):

		z_up = Z_vec[k]; z_low = Z_vec[k+1]
		z_mean,Bool = find_units_limit(model, xp, yp, z_up, z_low, 1E-7)

		if Bool:
			aux_val = Z_vec[k]
			Z_vec[k+1] = min(z_mean - aux_val,-1.01) + aux_val
			#Z_vec[k+1] = z_mean
			if (counter > 0):
				Aux_vec = np.linspace(Z_vec[K_index], Z_vec[k+1], counter+2)
				Z_vec[K_index:k+2] = Aux_vec
				counter = 0

			K_index = k+1
		else:
			counter += 1
			#Aux_Bool[k+1] = False

	if (counter > 0):
		Aux_vec = np.linspace(Z_vec[K_index], Z_vec[N_layers], counter+2)
		Z_vec[K_index:N_layers+1] = Aux_vec

# DEFINE_BOOTOMS: It finds the layers bottoms
#def define_bottoms(model, xp, yp, N_layers,	Z_vec):
	#return(Aux_Bool)
		




	
	# model.closest([1,1,1]) -> unidad geologica a la que eprtenece
	# model.height([1,1]) -> altura del punto ingresado
	# model.bbox -> cubo que contienen todo el modelo geologico
