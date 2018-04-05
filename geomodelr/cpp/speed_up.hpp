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

double distance_poly_fault_pt(const point2& pt, const polygon& poly,const rtree_seg* poly_seg_tree,
    const vector<rtree_seg *>& fault_lines);

double distance_poly_fault_pt2(int idx, const point2& pt, double x_rigth,const rtree_s* poly_seg_tree,
    const vector<rtree_seg *>& fault_lines);

#endif //GEOMODELR_SPEED_UP_HPP
