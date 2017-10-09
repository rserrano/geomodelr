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

# Creates the data from Geomodelr to ModFlow
def create_modflow_inputs( model_name='no_name', model=None, N_row=100,
	N_col=100, N_layers=100, properties=None, act_uni=None ,Bbox=None, Class=1):
	
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


	if (Bbox is None):
		Bbox = model.bbox

	del(num_properties)

	#properties[1,:] /= properties[0,:]
	if (act_uni is None):
		Units_data = { (model.units)[k]: np.append(properties[:,k], 1) for k in range(num_unit) }
	else:
		Units_data = { (model.units)[k]: np.append(properties[:,k],act_uni[k]) for k in range(num_unit) }

	X_inf = Bbox[0]
	X_sup = Bbox[3]

	Y_inf = Bbox[1]
	Y_sup = Bbox[4]

	dX = (X_sup - X_inf)/N_col
	dY = (Y_sup - Y_inf)/N_row

	X_vec = np.linspace(X_inf + dX/2., X_sup - dX/2.,N_col)
	Y_vec = np.linspace(Y_inf + dY/2., Y_sup - dY/2.,N_row)
 
 	Z_top = np.zeros((N_row,N_col))
 	Z_bottoms = np.zeros((N_layers,N_row,N_col))

	Bottom_min = Bbox[2] + 1E-5

	# Define Z-top

	for i in np.arange(N_row):
		for j in np.arange(N_col):

			xp = X_vec[j]
			yp = Y_vec[i]
			Z_top[i,j] = model.height([xp, yp])


	Z_top_min = np.min(Z_top)
 	# ---------- Define Z-top, Z-bottoms and Hydraulic conductivity
 	# Class Fine
	if (Class == 1):

		ibound = np.ones((N_layers, N_row, N_col), dtype=np.int32)

		K_hor = np.zeros((N_layers,N_row,N_col))
		K_ratio_hor = np.zeros((N_layers,N_row,N_col))
		K_ver = np.zeros((N_layers,N_row,N_col))

		for i in np.arange(N_row):
			yp = Y_vec[i]

			for j in np.arange(N_col):

				xp = X_vec[j]
				z_max = Z_top[i,j]

				z_vec = np.linspace(z_max, Bottom_min,N_layers+1)
				#define_bottoms(model, xp, yp, N_layers,	z_vec)

				Z_bottoms[:,i,j] = z_vec[1:]

				for L in np.arange(N_layers):
					zp = (z_vec[L] + z_vec[L+1])/2.0

					Unit = model.closest([xp,yp,zp])[0]
					Data = Units_data[Unit]
					K_hor[L,i,j] = Data[0]; K_ratio_hor[L,i,j] = Data[1]
					K_ver[L,i,j] = Data[2]
					I_bound[L,i,j] = Data[3]
	
	elif (Class==2):

	# ---------- Define Z-top, Z-bottoms and Hydraulic conductivity
 	# Class Smart

		Z_Bool_Top = np.isfinite(Z_top)

		Bottom_BOOL = np.isfinite(Z_bottoms)

		for i in np.arange(N_row):
			yp = Y_vec[i]

			for j in np.arange(N_col):

				xp = X_vec[j]
				z_max = Z_top[i,j]

				dz = (z_max - Bottom_min)/N_layers

				z_mean,change = find_units_limit(model, xp, yp, z_max,
					z_max-dz, 1E-3)

				Z_Bool_Top[i,j] = change

				Z_bottoms[0,i,j] = min(z_mean - z_max,-1.01) + z_max


		Z_Bool_Bot = np.isfinite(Z_top)
		Layers_Bool = np.isfinite(np.arange(N_layers))
		Bottom_BOOL[0,:,:] = Z_Bool_Top

		Last_Layer = 0

		for L in np.arange(1,N_layers-1):

			for i in np.arange(N_row):
				for j in np.arange(N_col):

					xp = X_vec[j]
					yp = Y_vec[i]
					
					zp = Z_bottoms[L-1,i,j]					
					#zl = zp - (zp- Bottom_min)/(N_layers-L)
					zl = (Z_top[i,j]*(N_layers-(L+1)) + Bottom_min*(L+1))/N_layers

					z_mean,change = find_units_limit(model, xp, yp, zp, zl,
						1E-3)

					Z_Bool_Bot[i,j] = change
					Z_bottoms[L,i,j] = min(z_mean - zp,-1.01) + zp

			if (np.sum(Z_Bool_Top & Z_Bool_Bot) == 0):

				Layers_Bool[L-1] = False
				Z_bottoms[L,Z_Bool_Top] = Z_bottoms[L-1,Z_Bool_Top]
				Z_Bool_Top = Z_Bool_Top | Z_Bool_Bot

			else:

				# Correction of FALSE values
				Aux_Bool=np.logical_not(Z_Bool_Top)
				Z_bottoms[L-1,Aux_Bool] = Z_bottoms[Last_Layer,Aux_Bool]
				Last_Layer = L

				Z_Bool_Top = Z_Bool_Bot.copy()

			Bottom_BOOL[L,:,:] = Z_Bool_Top

			# Correction of the layer besfore the last layer
			Aux_Bool=np.logical_not(Z_Bool_Top)
			Z_bottoms[L,Aux_Bool] = Z_bottoms[Last_Layer,Aux_Bool]



		N_layers = np.sum(Layers_Bool)
		Z_bottoms = Z_bottoms[Layers_Bool,:,:]
		Z_bottoms[-1,:,:] = Bottom_min

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
				yp = Y_vec[i]
				for j in np.arange(N_col):
		
					xp = X_vec[j]

					Unit = model.closest([xp,yp,mid_points[i,j]])[0]
					Data = Units_data[Unit]
					K_hor[L,i,j] = Data[0]; K_ratio_hor[L,i,j] = Data[1]
					K_ver[L,i,j] = Data[2]
					I_bound[L,i,j] = Data[3]

	else:
		pass

	# ---------- Define Z-top, Z-bottoms and Hydraulic conductivity
 	# Class Smart Top and Bottom



	#return((N_layers,Z_bottoms,K_hor))

	#  ------- Flopy Package ----
	# Grid

	print('Meshing processes: Done')

	chani_var = -np.ones(N_layers, dtype=np.int32)

	mf_handle = fp.modflow.Modflow(model_name, exe_name='mf2005',verbose=False)
	dis = fp.modflow.ModflowDis(mf_handle,nlay=N_layers, nrow=N_row, ncol=N_col,
		top=Z_top, botm=Z_bottoms, delc=dY, delr=dX, xul=X_inf, yul=Y_sup)


	# Variables for the BAS package
	# I_bound = np.ones((N_layers, N_row, N_col), dtype=np.int32)
	# ibound[:, :, 0] = -1
	# ibound[:, :, 1] = -2
	# ibound[:, :, 2] = -3
	# strt = np.ones((N_layers, N_row, N_col), dtype=np.float32)
	# strt[:, :, 0] = 10.5
	# strt[:, :, 1] = 5.5
	# strt[:, :, 2] = 1.5
	bas = fp.modflow.ModflowBas(mf_handle,ibound = I_bound)

	# Add LPF package to the MODFLOW model
	lpf = fp.modflow.ModflowLpf(mf_handle, chani=chani_var, hk=K_hor,
		vka=K_ver, hani=K_ratio_hor)
	
	# oc = fp.modflow.ModflowOc(mf_handle)
	# pcg = fp.modflow.ModflowPcg(mf_handle)

	mf_handle.write_input()
	
	mv_comand = 'mv ' + model_name + '* /media/sf_CompartidaVB/'
	os.system(mv_comand)

	return((dis,Layers_Bool,Bottom_BOOL))
	#return(dis)

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

# UNIIT_ZTOP_BOT: It finds the layers bottoms
def units_Ztop_bot(model, Top, Bottom, X_vec, Y_vec,
	num_unit, N_row, N_col, N_layers, N_div=None):
	
	Ztop_bot = { (model.units)[k]: np.array([np.inf,-np.inf,0.,0.]) for k in range(num_unit) }
	Lay_Range = { (model.units)[k]: np.array([N_layers,1]) for k in range(num_unit) }
	

	for i in np.arange(N_row):
		yp = Y_vec[i]

		for j in np.arange(N_col):
			xp = X_vec[j]
			z_top = Top[i,j]
			dz = (z_top-Bottom)/N_div
			dz_ray = (z_top-Bottom)/N_layers

			for k in np.arange(N_div+1):

				zp = z_top-k*dz
				Unit = model.closest([xp,yp,zp])[0]
				data = Ztop_bot[Unit]
				lay_data = Lay_Range[Unit]

				Bool = False
				if (zp<data[0]):
					data[0] = zp
					Bool = True
					data[2] = z_top

				if (zp>data[1]):
					data[1] = zp
					Bool = True
					#Lay_up = N_layers+1-ceil((zp - Bottom)/dz_ray)
					#lay_data[1] = Lay_up
					data[3] = z_top
				
				if Bool:
					
					Ztop_bot[Unit] = data

					#Lay_up = lay_data[1]
					#Lay_low = Lay_up + ceil((data[1]-data[0])/dz_ray)-1
					#lay_data[0] = Lay_low
	
					#Lay_Range[Unit] = lay_data

	return((Lay_Range,Ztop_bot))

		




	
	# model.closest([1,1,1]) -> unidad geologica a la que eprtenece
	# model.height([1,1]) -> altura del punto ingresado
	# model.bbox -> cubo que contienen todo el modelo geologico
