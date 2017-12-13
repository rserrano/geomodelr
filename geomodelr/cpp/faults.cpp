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


point3 cross_product(const point3& v1,const point3& v2){

    point3 output;
    output.set<0>(gy(v1)*gz(v2) - gy(v2)*gz(v1));
    output.set<1>(gz(v1)*gx(v2) - gz(v2)*gx(v1));
    output.set<2>(gx(v1)*gy(v2) - gx(v2)*gy(v1));

    return output;
}


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

	for (const triangle_pt& it: fplane){ //for each triangle in the fault plane.

		node_a = g0(it); node_b = g1(it); node_c = g2(it);
		eval_a = geometry::dot_product(node_a,nv) - D; // evals the first point of the triangle in the plane equation.
		eval_b = geometry::dot_product(node_b,nv) - D; // evals the second point of the triangle in the plane equation.
		eval_c = geometry::dot_product(node_c,nv) - D; // evals the third point of the triangle in the plane equation.
        line straight_segment;
        point2 aux_point;

        // determines if the straight segment given by the first and second point intersects the fault plane.
		if (eval_a*eval_b <= 0.0){

			dir = node_b;
			geometry::subtract_point(dir,node_a);
            dot_val = geometry::dot_product(nv,dir);
            if (std::abs(dot_val)>1e-50){
                T = -eval_a/dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_a);
                
                geometry::subtract_point(dir,x0);
                // saves the point coordinates in the coordinate system of the plane.
                geometry::append(straight_segment,point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2)));
            }
		}

        // determines if the straight segment given by the second and third point intersects the fault plane.
		if (eval_b*eval_c <= 0.0){

			dir = node_c;
			geometry::subtract_point(dir,node_b);
			dot_val = geometry::dot_product(nv,dir);
            if (std::abs(dot_val)>1e-50){

                T = -eval_b/dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_b);
                geometry::subtract_point(dir,x0);
                if (straight_segment.size()==1){

                    // Find if the point is already in the list.
                    aux_point = point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2));
                    if (geometry::distance(aux_point,straight_segment[0])>=2E-5){
                        geometry::append(straight_segment,aux_point); // saves the point coordinates in the coordinate system of the plane.
                    }
                }
                else{
                    // saves the point coordinates in the coordinate system of the plane.
                    geometry::append(straight_segment,point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2)));    
                }
                
            }
		}

        /* determines if the straight segment given by the third and first point intersects the fault plane as long as the
           straight segment has not yet been defined.*/
		if ((eval_c*eval_a <= 0.0) && (straight_segment.size()==1)) {

			dir = node_a;
			geometry::subtract_point(dir,node_c);
			dot_val = geometry::dot_product(nv,dir);

            if (std::abs(dot_val)>1e-50){

                T = -eval_c/dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_c);
                geometry::subtract_point(dir,x0);
                aux_point = point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2));
                
                // Find if the point is already in the list.
                if (geometry::distance(aux_point,straight_segment[0])>=2E-5){
                    geometry::append(straight_segment,aux_point); // saves the point coordinates in the coordinate system of the plane.
                }
            }
		}

		if (straight_segment.size()==2){
            vector<line> intersect_line;
            geometry::intersection(plane_poly,straight_segment,intersect_line);

            // ensures that the straight segment is greater than 2E-5 and, in addition, it must intersects the polygon of the plane.
            if ((intersect_line.size()==1) && (geometry::length(intersect_line[0])>=2E-5) ){
                output.push_back(intersect_line[0]);
            }
		}

	}
	return output;
}

// Finds the lines that intersect a plane with the fault planes.
// faults_cpp:    the set of fault planes to intersect with the plane.
// plane_info:      plane to intersect. It's the four corners of the plane.
void find_faults_plane_intersection(const map<wstring, vector<triangle_pt> >& faults_cpp, const pylist& plane_info,
    pydict& output, const int f_index, double& start_x) {

    // Plane information
    point3 x0, x_aux, v1, v2;
    
    x0.set<0>(python::extract<double>(plane_info[0][0])); x0.set<1>(python::extract<double>(plane_info[0][1]));
    x0.set<2>(python::extract<double>(plane_info[0][2]));

    // Find the vectors v1 and v2 in the plane that generate the space.
    v1.set<0>(python::extract<double>(plane_info[1][0]) - gx(x0)); v1.set<1>(python::extract<double>(plane_info[1][1]) - gy(x0));
    v1.set<2>(python::extract<double>(plane_info[1][2]) - gz(x0));
    double norm_v1 = std::sqrt( geometry::dot_product( v1, v1 ));

    v2.set<0>(python::extract<double>(plane_info[3][0]) - gx(x0)); v2.set<1>(python::extract<double>(plane_info[3][1]) - gy(x0));
    v2.set<2>(python::extract<double>(plane_info[3][2]) - gz(x0));

	double dot_val = geometry::dot_product(v1,v2);
    //Gram-Schmidt method is used to make them perpendicular.
	if (std::abs(dot_val)>1e-50){
		geometry::multiply_value( v1,dot_val/geometry::dot_product(v1,v1));
		geometry::subtract_point(v2,v1);	
	}
    // both vectors are normalized.
	geometry::divide_value( v1, std::sqrt( geometry::dot_product( v1, v1 )));
	geometry::divide_value( v2, std::sqrt( geometry::dot_product( v2, v2 )));

    polygon plane_poly; ring& outer = plane_poly.outer();
    outer.push_back(point2(0.0,0.0));

    // creates the polygon class using the four points of the plane.
    for (int k=1; k<python::len(plane_info); k++){
        x_aux.set<0>(python::extract<double>(plane_info[k][0]) - gx(x0)); x_aux.set<1>(python::extract<double>(plane_info[k][1]) - gy(x0));
        x_aux.set<2>(python::extract<double>(plane_info[k][2]) - gz(x0));
        outer.push_back(point2(geometry::dot_product(x_aux,v1),geometry::dot_product(x_aux,v2)));
    }

    // normal vector to the plane.
	point3 nv = cross_product(v1,v2);

    vector<line> faults_lines;

	for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){

        // finds the lines that intersect a plane with the fault planes.
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

// Finds the lines that intersect a set of planes with the fault planes.
// fplanes: the set of fault planes to intersect with the plane.
// planes: set of planes to intersect.
pydict find_faults_multiple_planes_intersection(const pydict& fplanes, const pylist& planes) {

    pylist dict_keys = fplanes.keys();

    map<wstring, vector<triangle_pt> > faults_cpp;

    // converts from pydict to map of C++.
    for (int k = 0; k<python::len(dict_keys); k++){
        wstring fkey = python::extract<wstring>(dict_keys[k]);
        
        vector<triangle_pt> faults_triangles;
        const pylist& dict_values = python::extract<pylist>((fplanes.values())[k]);

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

    pydict faults_intersection;
    double start_x = 0.0;
    for (int k=0; k<python::len(planes); k++){
        pylist plane_info = python::extract<pylist>(planes[k]);
        find_faults_plane_intersection(faults_cpp,plane_info,faults_intersection, k, start_x);
    }

    return faults_intersection;

}
 
// =========================================================================================

void get_triangle_normal(const triangle_pt& tri, const point3& X0, point3& nv){
    
    point3 v1 = g1(tri); geometry::subtract_point(v1,X0);
    point3 v2 = g2(tri); geometry::subtract_point(v2,X0);
    nv = cross_product(v1,v2);
}

bool intersect_triangle_plane(const triangle_pt& tri, const point3& x0,const point3& nv){

    double D = geometry::dot_product(x0,nv);
    double eval_a = geometry::dot_product(g0(tri),nv) - D; // evals the first point of the triangle in the plane equation.
    double eval_b = geometry::dot_product(g1(tri),nv) - D; // evals the second point of the triangle in the plane equation.
    double eval_c = geometry::dot_product(g2(tri),nv) - D; // evals the third point of the triangle in the plane equation.

    return ((eval_a*eval_b<=0) || (eval_b*eval_c<=0) || (eval_c*eval_a<=0));
}

vector<size_t> sort_indexes(const vector<double> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

pylist joint_lines_3d(vector<line_3d>& faults){

    pylist output;
    size_t C;
    double dist;
    point3 p0, pf, pa, pb;
    bool check; 

    while (faults.size()>0){

        // Takes the first straight segment.
        C=0; 
        p0 = (faults[0])[0];
        pf = (faults[0])[1];
        pylist line_p;
        line_p.append(python::make_tuple(gx(p0) ,gy(p0), gz(p0)));
        line_p.append(python::make_tuple(gx(pf) ,gy(pf), gz(pf)));
        faults.erase(faults.begin());

        while (C<faults.size()){

            for (int i=0; i<2; i++){

                check = true;
                pa = (faults[C])[i];
                // Computes the distace between pa and the initial point of the line.
                dist = geometry::distance(pa,p0);
                
                if (dist<1E-5){

                    pb = (faults[C])[1-i];
                    line_p.insert(0,python::make_tuple(gx(pb) ,gy(pb), gz(pb)));
                    p0 = pb;
                    faults.erase(faults.begin() + C);
                    C=0; check = false;
                    break;
                }

                // Computes the distace between pa and the final point of the line.
                dist = geometry::distance(pa,pf);
                if (dist<1E-5){

                    pb = (faults[C])[1-i];
                    line_p.append(python::make_tuple(gx(pb) ,gy(pb), gz(pb)));
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

// Finds the lines that intersect a fault plane with a given triangle.
void find_triangle_plane_intersection(const triangle_pt& tri, const point3& x0, const point3& nv, const point3& nv_line,
    vector<point3>& intersection_points, vector<double>& point_proyections){

    point3 node_a, node_b, node_c, dir;
    double D = geometry::dot_product(x0,nv);
    double eval_a, eval_b, eval_c, T, dot_val;

    node_a = g0(tri); node_b = g1(tri); node_c = g2(tri);
    eval_a = geometry::dot_product(node_a,nv) - D; // evals the first point of the triangle in the plane equation.
    eval_b = geometry::dot_product(node_b,nv) - D; // evals the second point of the triangle in the plane equation.
    eval_c = geometry::dot_product(node_c,nv) - D; // evals the third point of the triangle in the plane equation.

    // determines if the straight segment given by the first and second point intersects the fault plane.
    if (eval_a*eval_b <= 0.0){

        dir = node_b;
        geometry::subtract_point(dir,node_a);
        dot_val = geometry::dot_product(nv,dir);
        if (std::abs(dot_val)>1e-50){

            T = -eval_a/dot_val;
            geometry::multiply_value(dir,T); geometry::add_point(dir,node_a);            
            intersection_points.push_back(dir);
            point_proyections.push_back(geometry::dot_product(dir,nv_line));            
        }
    }

    // determines if the straight segment given by the second and third point intersects the fault plane.
    if (eval_b*eval_c <= 0.0){

        dir = node_c;
        geometry::subtract_point(dir,node_b);
        dot_val = geometry::dot_product(nv,dir);
        if (std::abs(dot_val)>1e-50){

            T = -eval_b/dot_val;
            geometry::multiply_value(dir,T); geometry::add_point(dir,node_b);
            intersection_points.push_back(dir);
            point_proyections.push_back(geometry::dot_product(dir,nv_line));            
        }
    }

    /* determines if the straight segment given by the third and first point intersects the fault plane as long as the
       straight segment has not yet been defined.*/
    if (eval_c*eval_a <= 0.0){

        dir = node_a;
        geometry::subtract_point(dir,node_c);
        dot_val = geometry::dot_product(nv,dir);

        if (std::abs(dot_val)>1e-50){

            T = -eval_c/dot_val;
            geometry::multiply_value(dir,T); geometry::add_point(dir,node_c);
            intersection_points.push_back(dir);
            point_proyections.push_back(geometry::dot_product(dir,nv_line));
        }
    }

}

pydict find_faults_topography_intersection(const pydict& fplanes, const pydict& topography_info, double up_faults){

    pylist dict_keys = fplanes.keys();

    map<wstring, vector<triangle_pt> > faults_cpp;

    double x_inf = python::extract<double>(topography_info["point"][0]);
    double y_inf = python::extract<double>(topography_info["point"][1]);
    point3 xy_0(x_inf,y_inf,-up_faults);

    double dx = python::extract<double>(topography_info["sample"][0]);
    double dy = python::extract<double>(topography_info["sample"][1]);

    int rows = python::extract<int>(topography_info["dims"][1]);
    int cols = python::extract<int>(topography_info["dims"][0]);

    // converts from pydict to map of C++.
    for (int k = 0; k<python::len(dict_keys); k++){
        wstring fkey = python::extract<wstring>(dict_keys[k]);
        
        vector<triangle_pt> faults_triangles;
        const pylist& dict_values = python::extract<pylist>((fplanes.values())[k]);

        for (int f=0; f< python::len(dict_values); f++){

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

    double topography_array[rows][cols];
    for (int i=0; i<rows;i++){
        for (int j=0; j<cols;j++){
            topography_array[i][j] = python::extract<double>(topography_info["heights"][j][i]);
        }
    }
    double max_x, max_y, min_x, min_y;
    int i_max, i_min, j_max, j_min;
    pydict output;

    for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){

        vector<line_3d> faults_intersection;
        for (const triangle_pt& tri_fault: iter->second){ //for each triangle in the fault plane.
            point3 x0_f = g0(tri_fault);
            point3 nv_f; get_triangle_normal(tri_fault,x0_f,nv_f);
            //geometry::subtract_point(x0_f,xy_0);

            max_x = std::max(gx(g0(tri_fault)),std::max(gx(g1(tri_fault)),gx(g2(tri_fault))))-x_inf;
            min_x = std::min(gx(g0(tri_fault)),std::min(gx(g1(tri_fault)),gx(g2(tri_fault))))-x_inf;
            max_y = std::max(gy(g0(tri_fault)),std::max(gy(g1(tri_fault)),gy(g2(tri_fault))))-y_inf;
            min_y = std::min(gy(g0(tri_fault)),std::min(gy(g1(tri_fault)),gy(g2(tri_fault))))-y_inf;

            // defines the bouding box of the fault triangle.
            j_min = std::max(int(ceil(min_x/dx)),1); j_max = std::min(int(ceil(max_x/dx)),cols-1);
            i_min = std::max(int(ceil(min_y/dy)),1); i_max = std::min(int(ceil(max_y/dy)),rows-1);

            for (int i=i_min; i<=i_max;i++){
                for (int j=j_min; j<=j_max;j++){

                    point3 B(dx*j,dy*(i-1),topography_array[i-1][j]); geometry::add_point(B,xy_0);
                    point3 C(dx*(j-1),dy*i,topography_array[i][j-1]); geometry::add_point(C,xy_0);

                    for (int t=0; t<=1;t++){

                        point3 A(dx*(j-1+t),dy*(i-1+t),topography_array[i-1+t][j-1+t]); geometry::add_point(A,xy_0);
                        triangle_pt tri_topo(A,B,C);
                        point3 nv_t; get_triangle_normal(tri_topo,A,nv_t);

                        point3 nv_line = cross_product(nv_t,nv_f);
                        double norm_nt_line = std::sqrt(geometry::dot_product(nv_line,nv_line));

                        if ((norm_nt_line>1E-20) && intersect_triangle_plane(tri_fault,A,nv_t) && intersect_triangle_plane(tri_topo,x0_f,nv_f)){
                            
                            vector<point3> intersection_points; vector<double> point_proyections;

                            find_triangle_plane_intersection(tri_fault, A, nv_t, nv_line, intersection_points, point_proyections);
                            find_triangle_plane_intersection(tri_topo, x0_f, nv_f, nv_line, intersection_points, point_proyections);

                            if (point_proyections.size()==4){
                                
                                vector<size_t> vec_pos = sort_indexes(point_proyections);
                                
                                if ((vec_pos[0]+vec_pos[1]!=1) && (vec_pos[0]+vec_pos[1]!=5)){

                                    line_3d segment;
                                    geometry::append(segment,intersection_points[vec_pos[1]]);
                                    geometry::append(segment,intersection_points[vec_pos[2]]);

                                    if (geometry::length(segment)>2E-5){
                                        faults_intersection.push_back(segment);
                                    }
                                }
                            }
                        }

                    }

                }
            }

        }

        output[iter->first] = joint_lines_3d(faults_intersection);

    }

    return output;
}