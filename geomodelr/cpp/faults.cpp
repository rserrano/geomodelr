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

pylist find_fault_plane_intersection(const vector<triangle_pt>& fplane, const point3& x0, const point3& v1, const point3& v2,
	const point3& nv) {

	point3 node_a, node_b, node_c, dir, ini;
	double D = geometry::dot_product(x0,nv);
	double eval_a, eval_b, eval_c, T;

	pylist output;

	for (const triangle_pt& it: fplane) {

		node_a = g0(it); node_b = g1(it); node_c = g2(it);
		eval_a = geometry::dot_product(node_a,nv) - D;
		eval_b = geometry::dot_product(node_b,nv) - D;
		eval_c = geometry::dot_product(node_c,nv) - D;
		pylist tuple_line;

		if (eval_a*eval_b < 1e-9){

			dir = node_b; ini = node_a;
			geometry::subtract_point(dir,ini);
			geometry::subtract_point(ini,x0);
			T = -geometry::dot_product(nv,ini)/geometry::dot_product(nv,dir);
			assert(T<1+1e-9); assert(T>-1e-9);

			geometry::multiply_value(dir,T); geometry::add_point(dir,node_a);
			tuple_line.append( python::make_tuple(gx(dir),gy(dir),gz(dir)) );
		}

		if (eval_b*eval_c < 1e-9){

			dir = node_c; ini = node_b;
			geometry::subtract_point(dir,ini);
			geometry::subtract_point(ini,x0);
			T = -geometry::dot_product(nv,ini)/geometry::dot_product(nv,dir);
			assert(T<1+1e-9); assert(T>-1e-9);

			geometry::multiply_value(dir,T); geometry::add_point(dir,node_b);
			tuple_line.append( python::make_tuple(gx(dir),gy(dir),gz(dir)) );
		}

		if (eval_c*eval_a < 0.0){

			dir = node_a; ini = node_c;
			geometry::subtract_point(dir,ini);
			geometry::subtract_point(ini,x0);
			T = -geometry::dot_product(nv,ini)/geometry::dot_product(nv,dir);
			assert(T<1+1e-9); assert(T>-1e-9);

			geometry::multiply_value(dir,T); geometry::add_point(dir,node_c);
			tuple_line.append( python::make_tuple(gx(dir),gy(dir),gz(dir)) );
		}

		if (python::len(tuple_line)==2){
			output.append(tuple_line);
		}

	}

	return output;
}

/*
pylist find_fault_plane_intersection(const pylist& fplanes, const point3& x0, const point3& vx, 
	const point3& vy, const point3& nv){
}*/

pydict find_faults_plane_intersection(const pydict& fplanes, const pylist& plane_info) {

	pylist dict_keys = fplanes.keys();
	//pylist dict_values = fplanes.values();

	map<wstring, vector<triangle_pt> > faults_cpp;

	//point3 node_a, node_b, node_c, x0, v1, v2;

	for (int k = 0; k<len(dict_keys); k++){
		//string fkey = string(python::extract<char*>(python::str(dict_keys[k])));
		wstring fkey    = python::extract<wstring>(dict_keys[k]);
		
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

		faults_cpp.insert(map<wstring, vector<triangle_pt>>::value_type(fkey, faults_triangles));
	}

	// Plane information
	point3 x0, v1, v2;
	
	x0.set<0>(python::extract<double>(plane_info[0][0])); x0.set<1>(python::extract<double>(plane_info[0][1]));
	x0.set<2>(python::extract<double>(plane_info[0][2]));

	v1.set<0>(python::extract<double>(plane_info[1][0]) - gx(x0)); v1.set<1>(python::extract<double>(plane_info[1][1]) - gy(x0));
	v1.set<2>(python::extract<double>(plane_info[1][2]) - gz(x0));	

	v2.set<0>(python::extract<double>(plane_info[3][0]) - gx(x0)); v2.set<1>(python::extract<double>(plane_info[3][1]) - gy(x0));
	v2.set<2>(python::extract<double>(plane_info[3][2]) - gz(x0));

	double dot_val = geometry::dot_product(v1,v2);
	if (std::abs(dot_val)>1e-5){
		geometry::multiply_value( v1,dot_val/geometry::dot_product(v1,v1));
		geometry::subtract_point(v2,v1);	
	}

	geometry::divide_value( v1, std::sqrt( geometry::dot_product( v1, v1 )));
	geometry::divide_value( v2, std::sqrt( geometry::dot_product( v2, v2 )));

	point3 nv;
	nv.set<0>(gy(v1)*gz(v2) - gy(v2)*gz(v1));
	nv.set<1>(gz(v1)*gx(v2) - gz(v2)*gx(v1));
	nv.set<2>(gx(v1)*gy(v2) - gx(v2)*gy(v1));

	//nv = geometry::

	pydict faults_intersection;

	for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){
		faults_intersection[iter->first] = find_fault_plane_intersection(iter->second, x0, v1, v2, nv);
	}

	return faults_intersection;

}
pydict find_faults_multiple_planes_intersection(const pydict& fplanes, const pylist& planes) {
	std::cout << "Calvoooo\n";
	return pydict();
}
 
