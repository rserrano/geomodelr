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
	
	num_units = 8
	X_inf = model.bbox[0]
	X_sup = model.bbox[3]

	Y_inf = model.bbox[1]
	Y_sup = model.bbox[4]

	X_vec = np.linspace(X_inf,X_sup,grid_divisions[0])
	Y_vec = np.linspace(Y_inf,Y_sup,grid_divisions[1])
 
 	Z_max = np.zeros((grid_divisions[1],grid_divisions[0]))

	for i in np.arange(grid_divisions[1]):
		for j in np.arange(grid_divisions[0]):
			Z_max[i,j] = model.height([X_vec[j],Y_vec[i]])

	return(0)

	# model es un model_from_file de geomodelr
	# propierties deberia ser diccionario
    

    # import geomodelr
	# model = geomodelr.model_from_file('./modelo.json')

	
	# model.closest([1,1,1]) -> unidad geologica a la que eprtenece
	# model.height([1,1]) -> altura del punto ingresado
	# model.bbox -> cubo que contienen todo el modelo geologico
