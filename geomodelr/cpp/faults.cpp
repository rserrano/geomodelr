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
#include <algorithm>    // std::min_element, std::max_elemen
#include <stdlib.h>
#include <iomanip>

// ====================== AUXILIAR FUNCTIONS =================================
// Converts the python dictionary of faults into a C++ map.
//      faults_python: python dictionary of faults.

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
//      planes: pylist of planes.
vector<line_3d> pylist_to_vector(const pylist& planes){

    vector<line_3d> planes_cpp;
    for (int i=0; i<python::len(planes); i++){
        int N = python::len(python::extract<pylist>(planes[i]));
        line_3d plane;
        for (int j=0; j<N; j++){
            plane.push_back(point3(python::extract<double>(planes[i][j][0]),python::extract<double>(planes[i][j][1]),
                python::extract<double>(planes[i][j][2])) );
        }
        planes_cpp.push_back(plane);
    }

    return planes_cpp;
}

// Converts a C++ vector of intersections into a python list.
//      input: C++ vector of intersections.
pylist vector_to_pylist(const vector<line>& input){
    pylist output;
    for (auto& it_line: input){
        pylist line_p;
        for (auto& it_point: it_line){
            pylist point; point.append(gx(it_point)); point.append(gy(it_point));
            //line_p.append(python::make_tuple(gx(it_point),gy(it_point)));
            line_p.append(point);
        }
        output.append(line_p);
    }
    return output;
}


pylist join_lines_tree_test(const pylist& segments){
    
    vector<line_segment> lines_cpp;
    for (int k=0; k<python::len(segments);k++){
        point2 a(python::extract<double>(segments[k][0][0]),python::extract<double>(segments[k][0][1]));
        point2 b(python::extract<double>(segments[k][1][0]),python::extract<double>(segments[k][1][1]));
        lines_cpp.push_back(line_segment(a,b));
    }
    return (vector_to_pylist(join_lines_tree(lines_cpp,0.0)));
}
// Converts the C++ map of intersections into a python dictionary.
//      intersections: C++ map of intersections.
pydict map_to_pydict(const map<wstring, vector<line> >& intersections){
    pydict output;
    for ( auto it = intersections.begin(); it != intersections.end(); it++ ){
        output[it->first] = vector_to_pylist(it->second);
    }
    return output;
}

// Computes the cross product between two point3.
//      v1: first vector.
//      v2: second vector.
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

vector<line> join_lines_tree(const vector<line_segment>& lines,const double start_x){

    vector<value_l> intersections(lines.size());
    int C=0;
    for (auto& iter_segment: lines){
        intersections.push_back(std::make_tuple(iter_segment,C));
        C++;
    }

    rtree_l segments_tree(intersections);
    vector<bool> check_lines(lines.size(), true);
    //vector<size_t> positions(lines.size()); iota(positions.begin(), positions.end(), 0);

    vector<line> output;
    point2 p0, pf, pa, pb, pt;
    point2 aux_pt(start_x,0.0);
    int pos;

    auto bool_it = std::find (check_lines.begin(), check_lines.end(), true);
    while (bool_it != check_lines.end()){

        size_t bool_pos = bool_it - check_lines.begin();
        // Takes the first straight segment.
        p0 = lines[bool_pos].first;
        pf = lines[bool_pos].second; 
        line line_p;
        line_p.push_back(point2(gx(p0)+start_x,gy(p0)));
        line_p.push_back(point2(gx(pf)+start_x,gy(pf)));
        check_lines[bool_pos] = false;

        for (int k=0; k<2; k++){
            if (k==0){
                pt=pf;
            } else {
                pt=p0;
            }
            bool check = true;
            while (check){

                box bx(point2(gx(pt) - epsilon, gy(pt) - epsilon), point2(gx(pt) + epsilon, gy(pt) + epsilon));
                check = false;
                for ( auto it = segments_tree.qbegin(geometry::index::intersects(bx));it != segments_tree.qend(); it++ ) {
                    pos = g1(*it);
                    if (check_lines[pos]){
                        pa = g0(*it).first;
                        pb = g0(*it).second;

                        if (geometry::distance(pa,pt)<epsilon){
                            pt = pb; geometry::add_point(pb,aux_pt);
                            if (k==0) {
                                line_p.push_back(pb);
                            }else {
                                line_p.insert(line_p.begin(),pb);
                            }
                            check_lines[pos] = false; check = true; break;

                        } else if (geometry::distance(pb,pt)<epsilon){
                            pt = pa; geometry::add_point(pa,aux_pt);
                            if (k==0){
                                line_p.push_back(pa);
                            }else {
                                line_p.insert(line_p.begin(),pa);
                            }

                            check_lines[pos] = false; check = true; break;
                        }
                    }
                }
            }

        }

        p0 = line_p[0];
        geometry::subtract_point(p0,line_p.back());
        if (std::sqrt(geometry::dot_product(p0,p0))<epsilon){
            line_p.pop_back();
        }
        output.push_back(line_p);
        bool_it = std::find (check_lines.begin(), check_lines.end(), true);
    }
    return output;
}
// Finds the lines that intersect a plane with a given fault.
//      fplane: triangles of the faulrt.
//      x0: point of the plane.
//      v1: first direction vector of the plane.
//      v2: second direction vector of the plane.
//      nv: normal vector of the plane.
//      plane_poly: polygonal representation of the plane using the four points of its corners.
vector<line_segment> find_fault_plane_intersection(const vector<triangle_pt>& fplane, const point3& x0, const point3& v1, const point3& v2,
	const point3& nv, const polygon& plane_poly) {

	point3 node_a, node_b, node_c, dir;
	double D = geometry::dot_product(x0,nv);
	double eval_a, eval_b, eval_c, T, dot_val;

	vector<line_segment> output;
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
            if (std::abs(dot_val)>tolerance){
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
            if (std::abs(dot_val)>tolerance){

                T = -eval_b/dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_b);
                geometry::subtract_point(dir,x0);
                if (straight_segment.size()==1){

                    // Find if the point is already in the list.
                    aux_point = point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2));
                    if (geometry::distance(aux_point,straight_segment[0])>=2*epsilon){
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

            if (std::abs(dot_val)>tolerance){

                T = -eval_c/dot_val;
                geometry::multiply_value(dir,T); geometry::add_point(dir,node_c);
                geometry::subtract_point(dir,x0);
                aux_point = point2(geometry::dot_product(dir,v1),geometry::dot_product(dir,v2));
                
                // Find if the point is already in the list.
                if (geometry::distance(aux_point,straight_segment[0])>=2*epsilon){
                    geometry::append(straight_segment,aux_point); // saves the point coordinates in the coordinate system of the plane.
                }
            }
		}

		if (straight_segment.size()==2){
            vector<line> intersect_line;
            geometry::intersection(plane_poly,straight_segment,intersect_line);

            // ensures that the straight segment is greater than 2*epsilon and, in addition, it must intersects the polygon of the plane.
            if ((intersect_line.size()==1) && (geometry::length(intersect_line[0])>=2*epsilon) ){
                output.push_back(line_segment(intersect_line[0][0],intersect_line[0][1]));
            }
		}

	}
	return output;
}

// Finds the lines that intersect a plane with the fault planes.
//      faults_cpp:    the set of fault planes to intersect with the plane.
//      plane_info:    plane to intersect. It's the four corners of the plane.
void find_faults_plane_intersection(const map<wstring, vector<triangle_pt> >& faults_cpp, const line_3d& plane_info,
    map<wstring, vector<line> >& output, const int f_index, double& start_x) {

    // Find the vectors v1 and v2 in the plane that generate the space.
    point3 x0 = plane_info[0];
    point3 v1 = plane_info[1]; geometry::subtract_point(v1,x0);
    point3 v2 = plane_info[3]; geometry::subtract_point(v2,x0);
    double norm_v1 = std::sqrt( geometry::dot_product( v1, v1 ));
	double dot_val = geometry::dot_product(v1,v2);

    //Gram-Schmidt method is used to make them perpendicular.
	if (std::abs(dot_val)>tolerance){
		geometry::multiply_value( v1,dot_val/geometry::dot_product(v1,v1));
		geometry::subtract_point(v2,v1);	
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
        vector<line_segment> fault_lines = find_fault_plane_intersection(iter->second, x0, v1, v2, nv, plane_poly);

        //vector<line> aux_vector = join_lines(fault_lines,start_x);
        vector<line> aux_vector = join_lines_tree(fault_lines,start_x);
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
map<wstring, vector<line>> find_faults_multiple_planes_intersection(const map<wstring, vector<triangle_pt>>& faults_cpp,
    const vector<line_3d>& planes_cpp) {

    map<wstring, vector<line>> faults_intersection;
    double start_x = 0.0;
    for (size_t k=0; k<planes_cpp.size(); k++){
        find_faults_plane_intersection(faults_cpp,planes_cpp[k],faults_intersection, k, start_x);
    }

    return faults_intersection;

}

// Finds the lines that intersect a set of planes with the fault planes.
//      fplanes: the set of fault planes to intersect with the plane.
//      planes: set of planes to intersect.
pydict find_faults_multiple_planes_intersection_python(const pydict& fplanes, const pylist& planes) {

    map<wstring, vector<triangle_pt> > faults_cpp = pydict_to_map(fplanes);
    vector<line_3d> planes_cpp = pylist_to_vector(planes);    
    return map_to_pydict(find_faults_multiple_planes_intersection(faults_cpp, planes_cpp));

}
 
// =========================================================================================

// Finds the normal vector of a triangle.
//      tri: triangle.
//      x0: point of the triangle.
point3 get_triangle_normal(const triangle_pt& tri, const point3& X0){
    
    point3 v1 = g1(tri); geometry::subtract_point(v1,X0);
    point3 v2 = g2(tri); geometry::subtract_point(v2,X0);
    return cross_product(v1,v2);

}

// Determines if a triangle can intersects with a plane.
//      tri: triangle.
//      x0: point of the plane.
//      nv: normal vector of the plane.
bool intersect_triangle_plane(const triangle_pt& tri, const point3& x0,const point3& nv){

    double D = geometry::dot_product(x0,nv);
    double eval_a = geometry::dot_product(g0(tri),nv) - D; // evals the first point of the triangle in the plane equation.
    double eval_b = geometry::dot_product(g1(tri),nv) - D; // evals the second point of the triangle in the plane equation.
    double eval_c = geometry::dot_product(g2(tri),nv) - D; // evals the third point of the triangle in the plane equation.

    return ((eval_a*eval_b<=0) || (eval_b*eval_c<=0) || (eval_c*eval_a<=0));
}

// Sorts the indexes of a vector<double> based on its values.
//      v: vector.
vector<size_t> sort_indexes(const vector<double> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// converts from pylist with the values of the topography to a C++ vector.
//      topography: pylist with the values of the topography.
//      rows: rows of the topography grid.
//      cols: cols of the topography grid.
//      z_max: reference variable to calculate the maximum value of the topography.
//      z_min: reference variable to calculate the minimum value of the topography.
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

vector<line> join_lines_3d(vector<line_3d>& faults){

    vector<line> output;
    size_t C;
    double dist;
    point3 p0, pf, pa, pb;
    bool check; 

    while (faults.size()>0){

        // Takes the first straight segment.
        C=0; 
        p0 = (faults[0])[0];
        pf = (faults[0])[1];
        line line_p;
        line_p.push_back(point2(gx(p0),gy(p0)));
        line_p.push_back(point2(gx(pf),gy(pf)));
        faults.erase(faults.begin());

        while (C<faults.size()){

            for (int i=0; i<2; i++){

                check = true;
                pa = (faults[C])[i];
                // Computes the distace between pa and the initial point of the line.
                dist = geometry::distance(pa,p0);
                
                if (dist<epsilon){

                    pb = (faults[C])[1-i];
                    line_p.insert(line_p.begin(),point2(gx(pb) ,gy(pb)));
                    p0 = pb;
                    faults.erase(faults.begin() + C);
                    C=0; check = false;
                    break;
                }

                // Computes the distace between pa and the final point of the line.
                dist = geometry::distance(pa,pf);
                if (dist<epsilon){

                    pb = (faults[C])[1-i];
                    line_p.push_back(point2(gx(pb) ,gy(pb)));
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
        output.push_back(line_p);
    }    
    return output;    
}

// Finds the lines that intersect a fault plane with a given triangle.
void find_triangle_plane_intersection(const triangle_pt& tri, const point3& x0, const point3& nv, const point3& nv_line,
    vector<point2>& intersection_points, vector<double>& point_proyections){

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
        if (std::abs(dot_val)>tolerance){

            T = -eval_a/dot_val;
            geometry::multiply_value(dir,T); geometry::add_point(dir,node_a);            
            intersection_points.push_back(point2(gx(dir),gy(dir)));
            point_proyections.push_back(geometry::dot_product(dir,nv_line));            
        }
    }

    // determines if the straight segment given by the second and third point intersects the fault plane.
    if (eval_b*eval_c <= 0.0){

        dir = node_c;
        geometry::subtract_point(dir,node_b);
        dot_val = geometry::dot_product(nv,dir);
        if (std::abs(dot_val)>tolerance){

            T = -eval_b/dot_val;
            geometry::multiply_value(dir,T); geometry::add_point(dir,node_b);
            intersection_points.push_back(point2(gx(dir),gy(dir)));
            point_proyections.push_back(geometry::dot_product(dir,nv_line));            
        }
    }

    /* determines if the straight segment given by the third and first point intersects the fault plane as long as the
       straight segment has not yet been defined.*/
    if ((eval_c*eval_a <= 0.0) && (intersection_points.size()==1 || intersection_points.size()==3)){

        dir = node_a;
        geometry::subtract_point(dir,node_c);
        dot_val = geometry::dot_product(nv,dir);

        if (std::abs(dot_val)>tolerance){

            T = -eval_c/dot_val;
            geometry::multiply_value(dir,T); geometry::add_point(dir,node_c);
            intersection_points.push_back(point2(gx(dir),gy(dir)));
            point_proyections.push_back(geometry::dot_product(dir,nv_line));
        }
    }

}

// Finds the lines that intersect the topography with a set of faults.
map<wstring, vector<line>> find_faults_topography_intersection(const map<wstring, vector<triangle_pt>>& faults_cpp,
    const vector<vector<double>>& topography_array, double z_max, double z_min,double x_inf, double y_inf,
    double dx, double dy, int rows, int cols){

    point3 xy_0(x_inf,y_inf,0.0);
    double max_x, max_y, min_x, min_y, max_z, min_z;
    int i_max, i_min, j_max, j_min;
    map<wstring, vector<line>> output;

    for (auto iter = faults_cpp.begin(); iter != faults_cpp.end(); iter++){

        vector<line_segment> faults_intersection;
        vector<value_l> intersections;

        for (const triangle_pt& tri_fault: iter->second){ //for each triangle in the fault plane.
            point3 x0_f = g0(tri_fault);
            point3 nv_f = get_triangle_normal(tri_fault,x0_f);

            max_z = std::max(gz(g0(tri_fault)),std::max(gz(g1(tri_fault)),gz(g2(tri_fault))));
            min_z = std::min(gz(g0(tri_fault)),std::min(gz(g1(tri_fault)),gz(g2(tri_fault))));

            // determines if the triangle intersects the range defined by z_min and z_max.
            if (std::min(max_z,z_max)>=std::max(min_z,z_min)){

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
                            point3 nv_t = get_triangle_normal(tri_topo,A);

                            point3 nv_line = cross_product(nv_t,nv_f);
                            double norm_nt_line = std::sqrt(geometry::dot_product(nv_line,nv_line));

                            // guarantees that the intersection can exist.
                            if ((norm_nt_line>tolerance) && intersect_triangle_plane(tri_fault,A,nv_t) && intersect_triangle_plane(tri_topo,x0_f,nv_f)){
                                
                                vector<point2> intersection_points; vector<double> point_proyections;

                                find_triangle_plane_intersection(tri_fault, A, nv_t, nv_line, intersection_points, point_proyections);
                                find_triangle_plane_intersection(tri_topo, x0_f, nv_f, nv_line, intersection_points, point_proyections);

                                if (point_proyections.size()==4){
                                    
                                    vector<size_t> vec_pos = sort_indexes(point_proyections);
                                    
                                    if ((vec_pos[0]+vec_pos[1]!=1) && (vec_pos[0]+vec_pos[1]!=5)){

                                        line_segment segment(intersection_points[vec_pos[1]],intersection_points[vec_pos[2]]);

                                        if (geometry::length(segment)>2*epsilon){
                                            faults_intersection.push_back(segment);
                                        }
                                    }
                                }
                            }

                        }

                    }
                }
            }

        }

        output[iter->first] = join_lines_tree(faults_intersection,0.0);
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
