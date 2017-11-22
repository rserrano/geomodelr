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
#include "faults.hpp"

vector<vector<point2>> find_fault_plane_intersection( const vector<triangle_pt>& fplane, const point3& x0, const point3& v1, const point3& v2, const point3& nv) {
	return vector<vector<point2>>();
}

pydict find_faults_plane_intersection(const pydict& fplanes, const pylist& plane) {
	std::cout << "executing find_faults_plane_intersection\n";
	return pydict();
}
pydict find_faults_multiple_planes_intersection(const pydict& fplanes, const pylist& planes) {
	std::cout << "executing find_faults_multiple_planes_intersection\n";
	return pydict();
}
 
