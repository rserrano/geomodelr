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

using namespace std;

/*vector<vector<point2>> find_fault_plane_intersection( const vector<triangle_pt>& fplane, const point3& x0, const point3& v1, const point3& v2, const point3& nv) {
	return vector<vector<point2>>();
}

pylist find_fault_plane_intersection(const pylist& fplanes, const point3& x0, const point3& vx, 
	const point3& vy, const point3& nv){
}*/

pydict find_faults_plane_intersection(const pydict& fplanes, const pylist& plane) {

	pylist dict_keys = fplanes.keys();
	//pylist dict_values = fplanes.values();

	map<string, vector<triangle_pt> > faults_cpp;

	point3 node_a, node_b, node_c;

	for (int k = 0; k<len(dict_keys); k++){
		//std::cout << python::extract<const char*>(python::str(dict_keys[k])) << std::endl;
		string fkey = string(python::extract<char*>(python::str(dict_keys[k])));
		wstring name    = python::extract<wstring>(dict_keys[k]);
		
		vector<triangle_pt> faults_triangles;
		pylist dict_values = python::extract<pylist>((fplanes.values())[k]);

		for (int f=0; f< len(dict_values); f++){

			point3 node_a(python::extract<double>(dict_values[f][0][0]),python::extract<double>(dict_values[f][0][1]),
				python::extract<double>(dict_values[f][0][2]));

			point3 node_b(python::extract<double>(dict_values[f][1][0]),python::extract<double>(dict_values[f][1][1]),
				python::extract<double>(dict_values[f][1][2]));

			point3 node_c(python::extract<double>(dict_values[f][2][0]),python::extract<double>(dict_values[f][2][1]),
				python::extract<double>(dict_values[f][2][2]));

			faults_triangles.push_back(triangle_pt(node_a,node_b,node_c));
		}

		//faults_cpp[fkey] = faults_triangles;
		faults_cpp.insert(map<string, vector<triangle_pt>>::value_type(fkey, faults_triangles));
	}

	cout << endl << "Size of the fault map:  " << faults_cpp.size() << endl << endl;
	for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){
		cout << endl << iter->first << endl;

		auto& aux_vec = iter->second;
		for (size_t k=0; k<aux_vec.size(); k++){

			triangle_pt& aux_tri = aux_vec[k];

			cout << gx(g0(aux_tri)) << " - " << gy(g0(aux_tri)) << " - " << gz(g0(aux_tri)) << "\t";
			cout << gx(g1(aux_tri)) << " - " << gy(g1(aux_tri)) << " - " << gz(g1(aux_tri)) << "\t";
			cout << gx(g2(aux_tri)) << " - " << gy(g2(aux_tri)) << " - " << gz(g2(aux_tri)) << endl;

		}
	}

	return pydict();
}
pydict find_faults_multiple_planes_intersection(const pydict& fplanes, const pylist& planes) {
	std::cout << "Calvoooo\n";
	return pydict();
}
 
