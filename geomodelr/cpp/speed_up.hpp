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
#ifndef GEOMODELR_SPEED_UP_HPP
#define GEOMODELR_SPEED_UP_HPP

#include "basic.hpp"
//#include "model.hpp"
#include <math.h>

using std::cout;
using std::endl;

class Model;

std::pair<double, bool> find_unit_limits_cpp(const Model* model,double xp, double yp,
	double z_max, double z_min, double eps);

template <typename Function, typename Derivative>
double Newton_Raphson(const Function& fx, const Derivative& diff_fx, double x0){
    double xf;
    for (size_t iter = 1; iter <= 10000; iter++){
        xf = x0 - fx(x0)/diff_fx(x0);
        if (std::abs(xf-x0) < tolerance){
            return xf;
        }
        x0 = xf;
    }
    return xf;
}

#endif //GEOMODELR_SPEED_UP_HPP
