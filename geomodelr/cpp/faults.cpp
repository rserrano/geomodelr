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
#include <algorithm>	// std::min_element, std::max_elemen
#include <stdlib.h>
#include <iomanip>

// ====================== AUXILIAR FUNCTIONS =================================
// Converts the python dictionary of faults into a C++ map.
//	  faults_python: python dictionary of faults.

map<wstring, vector<triangle_pt> > pydict_to_map(const pydict& faults_python){

	pylist dict_keys = faults_python.keys();

	map<wstring, vector<triangle_pt> > faults_cpp;

	for (int k = 0; k<python::len(dict_keys); k++){
		wstring fkey = python::extract<wstring>(dict_keys[k]);
		
		vector<triangle_pt> faults_triangles;
		const pylist& dict_values = python::extract<pylist>((faults_python.values())[k]);

		for (int f=0; f< python:: len(dict_values); f++){

			point3 node_a(python::extract<double>(dict_values[f][0][0]),python::extract<double>(dict_values[f][0][1]),
				python::extract<double>(dict_values[f][0][2]));

			point3 node_b(python::extract<double>(dict_values[f][1][0]),python::extract<double>(dict_values[f][1][1]),
				python::extract<double>(dict_values[f][1][2]));

			point3 node_c(python::extract<double>(dict_values[f][2][0]),python::extract<double>(dict_values[f][2][1]),
				python::extract<double>(dict_values[f][2][2]));

			faults_triangles.push_back(triangle_pt(node_a,node_b,node_c));
		}

		faults_cpp.insert(map<wstring, vector<triangle_pt>>::value_type(fkey, faults_triangles));
	}

	return faults_cpp;
}

// Converts a python list of planes into a C++ vector.
//	  planes: pylist of planes.
vector<line_3d> pylist_to_vector(const pylist& planes){

	vector<line_3d> planes_cpp;
	for (int i=0; i < python::len(planes); i++) {
		int N = python::len(python::extract<pylist>(planes[i]));
		line_3d plane;
		for (int j = 0; j < N; j++ ){
			plane.push_back(point3(	python::extract<double>(planes[i][j][0]),
						python::extract<double>(planes[i][j][1]),
						python::extract<double>(planes[i][j][2])) );
		}
		planes_cpp.push_back(plane);
	}
	return planes_cpp;
}



// Converts the C++ map of intersections into a python dictionary.
//	  intersections: C++ map of intersections.
pydict map_to_pydict(const map<wstring, vector<line> >& intersections){
	pydict output;
	for ( auto it = intersections.begin(); it != intersections.end(); it++ ){
		output[it->first] = vector_to_pylist(it->second);
	}
	return output;
}

// Joints the straight segments to create unique lines.
//	  faults: vector of straight segments of given fault.
//	  start_x: norm of v1 vector of the previous plane.
vector<line> join_lines(vector<line_segment>& faults, const double start_x){

	vector<line> output;
	size_t C;
	double dist;
	point2 p0, pf, pa, pb;
	bool check;
	point2 aux_point(start_x,0.0);

	while (!faults.empty()){

		// Takes the first straight segment.
		C=0; 
		p0 = faults[0].first; 
		pf = faults[0].second; 
		line line_p;
		line_p.push_back(point2(gx(p0)+start_x,gy(p0)));
		line_p.push_back(point2(gx(pf)+start_x,gy(pf)));
		faults.erase(faults.begin());

		while (C<faults.size()){

			for (int i=0; i<2; i++){

				check = true;
				if (i==0){
					pa = faults[C].first;
					pb = faults[C].second;
				} else {
					pb = faults[C].first;
					pa = faults[C].second;
				}
				// Computes the distace between pa and the initial point of the line.
				dist = geometry::distance(pa,p0);
				
				if (dist<epsilon){

					//pb = (faults[C])[1-i];
					p0 = pb;
					geometry::add_point(pb,aux_point);
					line_p.insert(line_p.begin(),pb);
					faults.erase(faults.begin() + C);
					C=0; check = false;
					break;
				}

				// Computes the distace between pa and the final point of the line.
				dist = geometry::distance(pa,pf);
				if (dist<epsilon){

					//pb = (faults[C])[1-i];
					pf=pb;
					geometry::add_point(pb,aux_point);
					line_p.push_back(pb);
					faults.erase(faults.begin() + C);
					C=0; check = false;
					break;
				}
			}

			if (check){
				C++;
			}
		}

		p0 = line_p[0];
		geometry::subtract_point(p0,line_p.back());
		if (std::sqrt(geometry::dot_product(p0,p0))<epsilon){
			line_p.pop_back();
		}
		output.push_back(line_p);
	}	
	return output;
	
}

// Finds the lines that intersect a plane with the fault planes.
//	  faults_cpp:	the set of fault planes to intersect with the plane.
//	  plane_info:	plane to intersect. It's the four corners of the plane.
void find_faults_plane_intersection(const map<wstring, vector<triangle_pt> >& faults_cpp, const line_3d& plane_info,
	map<wstring, vector<line> >& output, const int f_index, double& start_x) {

	// Find the vectors v1 and v2 in the plane that generate the space.
	point3 x0 = plane_info[0];
	point3 v1 = plane_info[1]; 
	geometry::subtract_point(v1, x0);
	point3 v2 = plane_info[3]; 
	geometry::subtract_point(v2, x0);
	double norm_v1 = std::sqrt( geometry::dot_product( v1, v1 ));
	double dot_val = geometry::dot_product(v1,v2);

	//Gram-Schmidt method is used to make them perpendicular.
	if (std::abs(dot_val)>tolerance){
		geometry::multiply_value( v1, dot_val/geometry::dot_product(v1,v1) );
		geometry::subtract_point( v2, v1 );
	}
	// both vectors are normalized.
	geometry::divide_value( v1, std::sqrt( geometry::dot_product( v1, v1 )));
	geometry::divide_value( v2, std::sqrt( geometry::dot_product( v2, v2 )));
	
	polygon plane_poly; ring& outer = plane_poly.outer();
	outer.push_back(point2(0.0,0.0));
	// creates the polygon class using the four points of the plane.
	for (size_t k=1; k<plane_info.size(); k++){
		point3 aux_p = plane_info[k]; geometry::subtract_point(aux_p,x0);
		outer.push_back(point2(geometry::dot_product(aux_p,v1),geometry::dot_product(aux_p,v2)));
	}

	// normal vector to the plane.
	point3 nv = cross_product(v1,v2);

	for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){

		// finds the lines that intersect a plane with the fault planes.
		vector<line_segment> fault_lines = find_mesh_plane_intersection(iter->second, x0, v1, v2, nv, plane_poly);
	
		vector<line> aux_vector = join_lines_tree(fault_lines,start_x,true);
		
	if (f_index==0){
			output.insert(map<wstring, vector<line>>::value_type(iter->first,aux_vector));
		}
		else{
			output[iter->first].reserve(output[iter->first].size() + distance(aux_vector.begin(),aux_vector.end()));
			output[iter->first].insert(output[iter->first].end(),aux_vector.begin(),aux_vector.end());
		}
	}

	start_x += norm_v1;

}

// Finds the lines that intersect a set of planes with the fault planes.
// fplanes: the set of fault planes to intersect with the plane.
// planes: set of planes to intersect.
map<wstring, vector<line>> find_faults_multiple_planes_intersection(const map<wstring, vector<triangle_pt>>& faults_cpp, const vector<line_3d>& planes_cpp) {
	
	map<wstring, vector<line>> faults_intersection;
	double start_x = 0.0;
	for (size_t k=0; k<planes_cpp.size(); k++){
		find_faults_plane_intersection(faults_cpp,planes_cpp[k],faults_intersection, k, start_x);
	}
	
	return faults_intersection;

}

// Finds the lines that intersect a set of planes with the fault planes.
//	  fplanes: the set of fault planes to intersect with the plane.
//	  planes: set of planes to intersect.
pydict find_faults_multiple_planes_intersection_python(const pydict& fplanes, const pylist& planes) {

	map<wstring, vector<triangle_pt> > faults_cpp = pydict_to_map(fplanes);
	vector<line_3d> planes_cpp = pylist_to_vector(planes);	
	return map_to_pydict(find_faults_multiple_planes_intersection(faults_cpp, planes_cpp));

}
 

// converts from pylist with the values of the topography to a C++ vector.
//	  topography: pylist with the values of the topography.
//	  rows: rows of the topography grid.
//	  cols: cols of the topography grid.
//	  z_max: reference variable to calculate the maximum value of the topography.
//	  z_min: reference variable to calculate the minimum value of the topography.
vector<vector<double>> topography_to_vector(const pylist& topography, int rows, int cols, double& z_max,double& z_min){

	vector<vector<double>> topography_array(rows, vector<double>(cols));

	z_min = python::extract<double>(topography[0][0]);
	z_max = z_min;
	double val;
	for (int i=0; i<rows;i++){
		for (int j=0; j<cols;j++){
			val = python::extract<double>(topography[j][i]);
			topography_array[i][j] = val;
			if (val<z_min){z_min=val;}
			if (val>z_max){z_max=val;}
		}
	}

	return topography_array;

}

// Finds the lines that intersect the topography with a set of faults.
map<wstring, vector<line>> find_faults_topography_intersection(const map<wstring, vector<triangle_pt>>& faults_cpp,
	const vector<vector<double>>& topography_array, double z_max, double z_min,double x_inf, double y_inf,
	double dx, double dy, int rows, int cols){
	map<wstring, vector<line>> output;

	for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){
		
		vector<line_segment> faults_intersection, tmp_intersection;
		
		for (const triangle_pt& tri_fault: iter->second){ //for each triangle in the fault plane.
			std::vector<segment> tmp_intersection = find_triangle_topography_intersection( tri_fault, topography_array, z_max, z_min, x_inf, y_inf, dx, dy, rows, cols );
			faults_intersection.insert(faults_intersection.end(), tmp_intersection.begin(), tmp_intersection.end());
		}
		
		output[iter->first] = join_lines_tree(faults_intersection,0.0,true);
		//output[iter->first] = join_lines(faults_intersection,0.0);
	}
	return output;
}

// Finds the lines that intersect the topography with a set of faults (Python).
pydict find_faults_topography_intersection_python(const pydict& fplanes, const pydict& topography_info){

	double x_inf = python::extract<double>(topography_info["point"][0]);
	double y_inf = python::extract<double>(topography_info["point"][1]);

	double dx = python::extract<double>(topography_info["sample"][0]);
	double dy = python::extract<double>(topography_info["sample"][1]);

	int cols = python::extract<int>(topography_info["dims"][0]);
	int rows = python::extract<int>(topography_info["dims"][1]);

	map<wstring, vector<triangle_pt> > faults_cpp = pydict_to_map(fplanes);

	double z_max, z_min;
	vector<vector<double>> topography_array = topography_to_vector(python::extract<pylist>(topography_info["heights"]),
		rows, cols, z_max, z_min);

	return map_to_pydict(find_faults_topography_intersection(faults_cpp,topography_array, z_max, z_min, x_inf, y_inf, dx, dy,
		rows, cols));
}
