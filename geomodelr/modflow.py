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
	N_layers=5, properties=None):
	
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

	X_vec = np.linspace(X_inf,X_sup,N_col + 1)
	Y_vec = np.linspace(Y_inf,Y_sup,N_row + 1)

	dX_vec = np.diff(X_vec)
	dY_vec = np.diff(Y_vec)
 
 	Z_top = np.zeros((N_row,N_col))
 	Z_bottoms = np.zeros((N_layers,N_row,N_col))

 	chani_var = -np.ones(N_layers, dtype=np.int32)

	K_hor = np.zeros((N_layers,N_row,N_col))
	K_ratio_hor = np.zeros((N_layers,N_row,N_col))
	K_ver = np.zeros((N_layers,N_row,N_col))

 	# ---------- Define Z-top, Z-bottoms and Hydraulic conductivity

	for i in np.arange(N_row):
		for j in np.arange(N_col):

			xp = X_vec[j] + dX_vec[j]/2.0
			yp = Y_vec[i] + dY_vec[i]/2.0
			z_max = model.height([xp, yp])
			Z_top[i,j] = z_max

			if (j>0):
				z_vec[0] = z_max
			elif (i > 0):
				z_vec = np.insert( Z_bottoms[:,i-1,0] ,0,z_max)
			else:
				z_vec = np.linspace(z_max,model.bbox[2],N_layers+1)

			z_vec = np.linspace(z_max,model.bbox[2],N_layers+1)

			#if (i==10)&(j==36):
			#	print z_vec[0]
			#	print xp,yp
			#	exit(0)

			define_bottoms(model, xp, yp, N_layers,	z_vec)

			#if (i==10)&(j==36):
			#	print z_vec[:3]
			#	print xp,yp
			#	exit(0)

			#Sum = np.sum(np.diff(z_vec)>-1.1)
			#if (Sum>0):
			#	print i, j
			#	exit(0)
			


			Z_bottoms[:,i,j] = z_vec[1:]

			for L in np.arange(N_layers):
				zp = (z_vec[L] + z_vec[L+1])/2.0

				Unit = model.closest([xp,yp,zp])[0]
				Data = K_data[Unit]
				K_hor[L,i,j] = Data[0]; K_ratio_hor[L,i,j] = Data[1]
				K_ver[L,i,j] = Data[2]


	#  ------- Flopy Package ----
	# Grid

	mf_handle = fp.modflow.Modflow(model_name, exe_name='mf2005',verbose=False)
	dis = fp.modflow.ModflowDis(mf_handle,nlay=N_layers, nrow=N_row, ncol=N_col,
		top=Z_top, botm=Z_bottoms, delc=dY_vec, delr=dX_vec)


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
def find_unit_boundary(model, xp, yp, Z_top, z_min, eps):

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


def define_bottoms(model, xp, yp, N_layers,	Z_vec):
	
	K_index = 0
	counter = 0

	#Aux_Bool = np.isfinite(Z_vec)

	for k in np.arange(N_layers-1):

		z_up = Z_vec[k]; z_low = Z_vec[k+1]
		z_mean,Bool = find_unit_boundary(model, xp, yp, z_up, z_low, 1E-7)

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

	#return(Aux_Bool)
		




	
	# model.closest([1,1,1]) -> unidad geologica a la que eprtenece
	# model.height([1,1]) -> altura del punto ingresado
	# model.bbox -> cubo que contienen todo el modelo geologico
