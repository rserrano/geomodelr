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

def create_modflow_inputs( outname, model, grid_divisions, properties):
	
	num_units = 3
	X_inf = model.bbox[0]
	X_sup = model.bbox[3]

	Y_inf = model.bbox[1]
	Y_sup = model.bbox[4]

	N_row = grid_divisions[0] + 1
	N_col = grid_divisions[1] + 1

	X_vec = np.linspace(X_inf,X_sup,N_col)
	Y_vec = np.linspace(Y_inf,Y_sup,N_row)
 
 	Z_max = np.zeros((N_row,N_col))

	for i in np.arange(N_row):
		for j in np.arange(N_col):
			Z_max[i,j] = model.height([X_vec[j],Y_vec[i]])

	dz = 0.5*( (X_sup-X_inf)/(N_col-1) + (Y_sup-Y_inf)/(N_row-1)) 

	Z_bottoms = np.zeros((num_units,N_row,N_col))

	for k in np.arange(num_units):
		Z_bottoms[k,:,:] = Z_max - dz*(k+1)



	#  ------- Flopy Package ----
	# Grid
	
	modelname = '/media/sf_CompartidaVB/Proof_Mesh'
	mf_handle = fp.modflow.Modflow(modelname, exe_name='mf2005',verbose=True)
	dis = fp.modflow.ModflowDis(mf_handle,nlay=num_units, nrow=N_row, ncol=N_col,top=Z_max, botm=Z_bottoms)

	# Variables for the BAS package
	ibound = np.ones((num_units, N_row, N_col), dtype=np.int32)
	ibound[:, :, 0] = -1
	ibound[:, :, -1] = -1
	strt = np.ones((num_units, N_row, N_col), dtype=np.float32)
	strt[:, :, 0] = 10.
	strt[:, :, -1] = 0.
	bas = fp.modflow.ModflowBas(mf_handle, ibound=ibound, strt=strt)

	# Add LPF package to the MODFLOW model
	lpf = fp.modflow.ModflowLpf(mf_handle, hk=10., vka=10.)
	oc = fp.modflow.ModflowOc(mf_handle)
	pcg = fp.modflow.ModflowPcg(mf_handle)

	#lst=fp.modflow.mf.ModflowList(mf_handle,extension='list',unitnumber=2)

	mf_handle.write_input()
	return(dis)

	# model es un model_from_file de geomodelr
	# propierties deberia ser diccionario
    

    # import geomodelr
	# model = geomodelr.model_from_file('./modelo.json')

	
	# model.closest([1,1,1]) -> unidad geologica a la que eprtenece
	# model.height([1,1]) -> altura del punto ingresado
	# model.bbox -> cubo que contienen todo el modelo geologico
