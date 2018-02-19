/*
	Geomodelr query tool. Tools for using a geomodelr.com geological model.
	Copyright (C) 2016 Geomodelr, Inc.
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Affero General Public License for more details.
	
	You should have received a copy of the GNU Affero General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include "speed_up.hpp"
#include <algorithm>    // std::min_element, std::max_elemen
#include <stdlib.h>
#include <iomanip>

std::pair<double, bool> find_unit_limits_cpp(const Model* model, double xp, double yp,
	double z_max, double z_min, double eps) {
  
	wstring unit_max = g0(model->closest(point3(xp,yp,z_max)));

	if ( unit_max.compare(g0(model->closest(point3(xp,yp,z_min))))==0 ){
		return std::make_pair( z_min, false );
	}    
    do {
    	double z_mean = (z_max+z_min)/2.0;
    	if ( unit_max.compare( g0(model->closest(point3(xp, yp, z_mean)))) == 0 ) {
    		z_max = z_mean;
    	}
		else {
    		z_min = z_mean;
    	}
    } while( (z_max-z_min)>eps );

    return std::make_pair( z_min, true );
}
