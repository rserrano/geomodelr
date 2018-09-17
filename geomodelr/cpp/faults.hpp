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
#ifndef GEOMODELR_FAULTS_HPP
#define GEOMODELR_FAULTS_HPP

#include "basic.hpp"
#include <math.h>
using std::cout;
using std::endl;

void find_faults_plane_intersection(const map<wstring, vector<triangle_pt> >& faults_cpp, const line_3d& plane_info,
	map<wstring, vector<line> >& faults_intersection,const int f_index, double& start_x);

pydict map_to_pydict(const map<wstring, vector<line> >& intersections);

map<wstring, vector<line>> find_faults_multiple_planes_intersection(const map<wstring, vector<triangle_pt>>& faults_cpp,
    const vector<line_3d>& planes_cpp);

pydict find_faults_multiple_planes_intersection_python(const pydict& fplanes, const pylist& planes);

map<wstring, vector<line>> find_faults_topography_intersection(const map<wstring, vector<triangle_pt>>& faults_cpp,
    const vector<vector<double>>& topography_array, double z_max, double z_min,double x_inf, double y_inf,
    double dx, double dy, int rows, int cols);

vector<vector<double>> topography_to_vector(const pylist& topography, int rows, int cols, double& z_max,double& z_min);

pydict find_faults_topography_intersection_python(const pydict& fplanes, const pydict& topography_info);

#endif // GEOMODELR_FAULTS_HPP
