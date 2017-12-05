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

// Joints the straight segments to create unique lines.
//      faults: vector of straight segments of given fault.
//      start_x: norm of v1 vector of the previous plane.
pylist joint_lines(vector<line>& faults, const double start_x){

	pylist output;
	size_t C;
	double dist;
	point2 p0, pf, pa, pb;
	bool check; 

    while (faults.size()>0){

        // Takes the first straight segment.
    	C=0; 
    	p0 = (faults[0])[0];
    	pf = (faults[0])[1];
    	pylist line_p;
    	line_p.append(python::make_tuple(gx(p0) + start_x,gy(p0)));
    	line_p.append(python::make_tuple(gx(pf) + start_x,gy(pf)));
    	faults.erase(faults.begin());

    	while (C<faults.size()){

    		for (int i=0; i<2; i++){

    			check = true;
    			pa = (faults[C])[i];
                // Computes the distace between pa and the initial point of the line.
    			dist = geometry::distance(pa,p0);
    			
    			if (dist<1E-5){

                    pb = (faults[C])[1-i];
    				line_p.insert(0,python::make_tuple(gx(pb) + start_x,gy(pb)));
    				p0 = pb;
					faults.erase(faults.begin() + C);
					C=0; check = false;
					break;
    			}

                // Computes the distace between pa and the final point of the line.
                dist = geometry::distance(pa,pf);
                if (dist<1E-5){

                    pb = (faults[C])[1-i];
                    line_p.append(python::make_tuple(gx(pb) + start_x,gy(pb)));
                    pf = pb;
                    faults.erase(faults.begin() + C);
                    C=0; check = false;
                    break;
                }
    		}

            if (check){
                C++;
            }
    	}
        output.append(line_p);
    }    
    return output;
    
}

// Finds the lines that intersect a plane with a given fault plane.
vector<line> find_fault_plane_intersection(const vector<triangle_pt>& fplane, const point3& x0, const point3& v1, const point3& v2,
	const point3& nv, const polygon& plane_poly) {

	point3 node_a, node_b, node_c, dir;
	double D = geometry::dot_product(x0,nv);
	double eval_a, eval_b, eval_c, T, dot_val;

	vector<line> output;

	for (const triangle_pt& it: fplane) { //for each triangle in the fault plane.

		node_a = g0(it); node_b = g1(it); node_c = g2(it);
		eval_a = geometry::dot_product(node_a,nv) - D; // evals the first point of the triangle in the plane equation.
		eval_b = geometry::dot_product(node_b,nv) - D; // evals the second point of the triangle in the plane equation.
		eval_c = geometry::dot_product(node_c,nv) - D; // evals the third point of the triangle in the plane equation.
        line straight_segment;

        // determines if the straight segment given by the first and second point intersects the fault plane.
		if (eval_a*eval_b <= 0.0){

			dir = node_b;
			geometry::subtract_point(dir,node_a);
            dot_val = geometry::dot_product(nv,dir);
            if (std::abs(dot_val)>1e-20){
                T = (D - geometry::dot_product(nv,node_a)) /dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_a);

                geometry::subtract_point(dir,x0);
                geometry::append(straight_segment,point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2)));
            }
		}

        // determines if the straight segment given by the second and third point intersects the fault plane.
		if (eval_b*eval_c <= 0.0){

			dir = node_c;
			geometry::subtract_point(dir,node_b);
			dot_val = geometry::dot_product(nv,dir);
            if (std::abs(dot_val)>1e-20){
                T = (D - geometry::dot_product(nv,node_b)) /dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_b);

                geometry::subtract_point(dir,x0);
                geometry::append(straight_segment,point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2)));
            }
		}

        // determines if the straight segment given by the third and first point intersects the fault plane.
		if (eval_c*eval_a <= 0.0){

			dir = node_a;
			geometry::subtract_point(dir,node_c);
			dot_val = geometry::dot_product(nv,dir);
            if (std::abs(dot_val)>1e-20){
                T = (D - geometry::dot_product(nv,node_c)) /dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_c);

                geometry::subtract_point(dir,x0);
                geometry::append(straight_segment,point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2)));
            }
		}

		if (straight_segment.size()==2){

            vector<line> intersect_line;
            geometry::intersection(plane_poly,straight_segment,intersect_line);
            if ((intersect_line.size()==1) && (geometry::length(intersect_line[0])>=5E-6) ){
                output.push_back(intersect_line[0]);
            }
		}

	}

	return output;
}


void find_faults_plane_intersection(const map<wstring, vector<triangle_pt> >& faults_cpp, const pylist& plane_info,
    pydict& output, const int f_index, double& start_x) {

    // Plane information
    point3 x0, x_aux, v1, v2;
    
    x0.set<0>(python::extract<double>(plane_info[0][0])); x0.set<1>(python::extract<double>(plane_info[0][1]));
    x0.set<2>(python::extract<double>(plane_info[0][2]));

    v1.set<0>(python::extract<double>(plane_info[1][0]) - gx(x0)); v1.set<1>(python::extract<double>(plane_info[1][1]) - gy(x0));
    v1.set<2>(python::extract<double>(plane_info[1][2]) - gz(x0));
    double norm_v1 = std::sqrt( geometry::dot_product( v1, v1 ));

    v2.set<0>(python::extract<double>(plane_info[3][0]) - gx(x0)); v2.set<1>(python::extract<double>(plane_info[3][1]) - gy(x0));
    v2.set<2>(python::extract<double>(plane_info[3][2]) - gz(x0));

	double dot_val = geometry::dot_product(v1,v2);
	if (std::abs(dot_val)>1e-9){
		geometry::multiply_value( v1,dot_val/geometry::dot_product(v1,v1));
		geometry::subtract_point(v2,v1);	
	}
	geometry::divide_value( v1, std::sqrt( geometry::dot_product( v1, v1 )));
	geometry::divide_value( v2, std::sqrt( geometry::dot_product( v2, v2 )));

    polygon plane_poly; ring& outer = plane_poly.outer();
    outer.push_back(point2(0.0,0.0));

    for (int k=1; k<python::len(plane_info); k++){
        x_aux.set<0>(python::extract<double>(plane_info[k][0]) - gx(x0)); x_aux.set<1>(python::extract<double>(plane_info[k][1]) - gy(x0));
        x_aux.set<2>(python::extract<double>(plane_info[k][2]) - gz(x0));
        outer.push_back(point2(geometry::dot_product(x_aux,v1),geometry::dot_product(x_aux,v2)));
    }

	point3 nv;
	nv.set<0>(gy(v1)*gz(v2) - gy(v2)*gz(v1));
	nv.set<1>(gz(v1)*gx(v2) - gz(v2)*gx(v1));
	nv.set<2>(gx(v1)*gy(v2) - gx(v2)*gy(v1));

	//nv = geometry::
    vector<line> faults_lines;

	for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){
        faults_lines = find_fault_plane_intersection(iter->second, x0, v1, v2, nv, plane_poly);
        if (f_index==0){
            output[iter->first] = pylist();
        }
        pylist aux_list = python::extract<pylist>(output[iter->first]);
        aux_list.extend(joint_lines( faults_lines, start_x ));
        output[iter->first] = aux_list;
	}

    start_x += norm_v1;

}

pydict find_faults_multiple_planes_intersection(const pydict& fplanes, const pylist& planes) {

    pylist dict_keys = fplanes.keys();

    map<wstring, vector<triangle_pt> > faults_cpp;

    for (int k = 0; k<python::len(dict_keys); k++){
        //string fkey = string(python::extract<char*>(python::str(dict_keys[k])));
        wstring fkey = python::extract<wstring>(dict_keys[k]);
        
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

    pydict faults_intersection;
    

    double start_x = 0.0;
    for (int k=0; k<python::len(planes); k++){
        pylist plane_info = python::extract<pylist>(planes[k]);
        find_faults_plane_intersection(faults_cpp,plane_info,faults_intersection, k, start_x);
    }

    return faults_intersection;

}
 
