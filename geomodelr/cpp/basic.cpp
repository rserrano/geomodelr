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
#include "basic.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

// Is geomodelr verbose.
bool geomodelr_verbose = false;

std::tuple<point2, double> point_segment_projection( const point2& pt, const point2& ps, const point2& pe ) 
{
	double x  = gx(pt);
	double y  = gy(pt);
	double x1 = gx(ps);
	double y1 = gy(ps);
	double x2 = gx(pe);
	double y2 = gy(pe);
	
	double A = x - x1;
	double B = y - y1;
	double C = x2 - x1;
	double D = y2 - y1;
	
	double dot = A * C + B * D;
	double len_sq = C * C + D * D;
	double param = dot / len_sq;
	
	double xx, yy;
	
	if (param < 0 || (x1 == x2 && y1 == y2)) {
	  xx = x1;
	  yy = y1;
	}
	else if (param > 1) {
	  xx = x2;
	  yy = y2;
	}
	else {
	  xx = x1 + param * C;
	  yy = y1 + param * D;
	}
	
	double dx = x - xx;
	double dy = y - yy;
	double dst = dx * dx + dy * dy;
	return std::make_tuple( point2( xx, yy ), dst );
}

std::tuple<point2, double> point_line_projection( const point2& pt, const line& l ) 
{
	std::tuple<point2, double> closest = std::make_tuple(l[0], std::numeric_limits<double>::infinity());
	for ( size_t i = 0; i < l.size()-1; i++ ) {
		std::tuple<point2, double> cl_segment = point_segment_projection( pt, l[i], l[i+1] );
		if ( g1(cl_segment) < g1(closest) ) {
			closest = cl_segment;
		}
	}
	double d = std::sqrt( g1(closest) );
	double p = (d + epsilon)/d;
	double xp = gx(g0(closest));
	double yp = gy(g0(closest));
	double xo = gx(pt);
	double yo = gy(pt);
	double xe = xo + p*(xp-xo);
	double ye = yo + p*(yp-yo);
	return std::make_tuple( point2(xe, ye), d+epsilon );
}

bool always_true( const value_f& v )
{
	return true;
}

// Sorts the indexes of a vector<double> based on its values.
//	  v: vector.
vector<size_t> sort_indexes(const vector<double> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
	   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// Finds the lines that intersect a plane with a set of triangles.
//	  triangles: triangles list.
//	  x0: point of the plane.
//	  v1: first direction vector of the plane.
//	  v2: second direction vector of the plane.
//	  nv: normal vector of the plane.
//	  plane_poly: polygonal representation of the plane using the four points of its corners.
vector<line_segment> find_mesh_plane_intersection(const vector<triangle_pt>& fplane, const point3& x0, const point3& v1, const point3& v2,
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

pylist join_lines_tree_test(const pylist& segments){
	
	vector<line_segment> lines_cpp;
	for (int k=0; k<python::len(segments);k++){
		point2 a(python::extract<double>(segments[k][0][0]),python::extract<double>(segments[k][0][1]));
		point2 b(python::extract<double>(segments[k][1][0]),python::extract<double>(segments[k][1][1]));
		lines_cpp.push_back(line_segment(a,b));
	}
	return (vector_to_pylist(join_lines_tree(lines_cpp,0.0,true)));
}

vector<line> join_lines_tree(const vector<line_segment>& lines, double start_x, bool check_loop){

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
	if ( check_loop ) {
			p0 = line_p[0];
			geometry::subtract_point(p0,line_p.back());
	
			if (std::sqrt(geometry::dot_product(p0,p0))<epsilon){
			line_p.pop_back();
			}
	}
		output.push_back(line_p);
		bool_it = std::find (check_lines.begin(), check_lines.end(), true);
	}
	return output;
}

pylist find_mesh_plane_intersection_python(const pylist& triangles, const pylist& plane_poly) {
	vector<triangle_pt> mesh;
	int n = python::len(triangles);
	for ( int i = 0; i < n; i++ ) {
		point3 node_a(python::extract<double>(triangles[i][0][0]),
				  python::extract<double>(triangles[i][0][1]), 
				  python::extract<double>(triangles[i][0][2]));
		
		point3 node_b(python::extract<double>(triangles[i][1][0]),
				  python::extract<double>(triangles[i][1][1]),
				  python::extract<double>(triangles[i][1][2]));
		
		point3 node_c(python::extract<double>(triangles[i][2][0]),
				  python::extract<double>(triangles[i][2][1]),
				  python::extract<double>(triangles[i][2][2]));
		
		mesh.push_back(triangle_pt(node_a, node_b, node_c));
	}
	n = python::len(plane_poly);
		line_3d plane;
		for (int j=0; j<n; j++){
		plane.push_back(point3(python::extract<double>(plane_poly[j][0]),
					   python::extract<double>(plane_poly[j][1]), 
					   python::extract<double>(plane_poly[j][2])) );
		}
	// Find the vectors v1 and v2 in the plane that generate the space.
	point3 x0 = plane[0];
	point3 v1 = plane[1]; 
	geometry::subtract_point(v1, x0);
	point3 v2 = plane[3]; 
	geometry::subtract_point(v2, x0);
	
	double dot_val = geometry::dot_product(v1,v2);

	//Gram-Schmidt method is used to make them perpendicular.
	if (std::abs(dot_val)>tolerance){
		geometry::multiply_value( v1,dot_val/geometry::dot_product(v1,v1));
		geometry::subtract_point(v2,v1);
	}
	// both vectors are normalized.
	geometry::divide_value( v1, std::sqrt( geometry::dot_product( v1, v1 )));
	geometry::divide_value( v2, std::sqrt( geometry::dot_product( v2, v2 )));
	
	polygon projected_poly; 
	ring& outer = projected_poly.outer();
	outer.push_back(point2(0.0,0.0));
	// creates the polygon class using the four points of the plane.
	for (size_t k=1; k < plane.size(); k++){
		point3 aux_p = plane[k]; 
		geometry::subtract_point(aux_p,x0);
		outer.push_back(point2(geometry::dot_product(aux_p,v1),geometry::dot_product(aux_p,v2)));
	}
	
	// normal vector to the plane.
	point3 nv = cross_product(v1,v2);
	
	vector<line_segment> out_lines = find_mesh_plane_intersection(mesh, x0, v1, v2, nv, projected_poly);
	
	vector<line> aux_vector = join_lines_tree(out_lines, 0.0, false);
	pylist out;
	for ( const line& l: aux_vector ) {
		pylist ln;
		for ( const point2& p: l ) {
			ln.append( python::make_tuple( gx( p ), gy( p ) ) );
		}
		out.append( ln );
	}
	return out;
}

// Computes the cross product between two point3.
//	  v1: first vector.
//	  v2: second vector.
point3 cross_product(const point3& v1,const point3& v2){
	point3 output;
	output.set<0>(gy(v1)*gz(v2) - gy(v2)*gz(v1));
	output.set<1>(gz(v1)*gx(v2) - gz(v2)*gx(v1));
	output.set<2>(gx(v1)*gy(v2) - gx(v2)*gy(v1));

	return output;
}

// Converts a C++ vector of intersections into a python list.
//	  input: C++ vector of intersections.
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

// =========================================================================================

// Finds the normal vector of a triangle.
//	  tri: triangle.
//	  x0: point of the triangle.
point3 get_triangle_normal(const triangle_pt& tri, const point3& X0){
	
	point3 v1 = g1(tri); geometry::subtract_point(v1,X0);
	point3 v2 = g2(tri); geometry::subtract_point(v2,X0);
	return cross_product(v1,v2);

}

// Determines if a triangle can intersects with a plane.
//	  tri: triangle.
//	  x0: point of the plane.
//	  nv: normal vector of the plane.
bool intersect_triangle_plane(const triangle_pt& tri, const point3& x0,const point3& nv){

	double D = geometry::dot_product(x0,nv);
	double eval_a = geometry::dot_product(g0(tri),nv) - D; // evals the first point of the triangle in the plane equation.
	double eval_b = geometry::dot_product(g1(tri),nv) - D; // evals the second point of the triangle in the plane equation.
	double eval_c = geometry::dot_product(g2(tri),nv) - D; // evals the third point of the triangle in the plane equation.

	return ((eval_a*eval_b<=0) || (eval_b*eval_c<=0) || (eval_c*eval_a<=0));
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



vector<line_segment> find_triangle_topography_intersection( const triangle_pt& tri_fault, 
								const vector<vector<double>>& topography_array, 
								double z_max, double z_min,double x_inf, double y_inf,
								double dx, double dy, int rows, int cols ) {
	
	
	point3 xy_0(x_inf,y_inf,0.0);
	double max_x, max_y, min_x, min_y, max_z, min_z;
	int i_max, i_min, j_max, j_min;
	
	point3 x0_f = g0(tri_fault);
	point3 nv_f = get_triangle_normal(tri_fault,x0_f);
	
	max_z = std::max(gz(g0(tri_fault)),std::max(gz(g1(tri_fault)),gz(g2(tri_fault))));
	min_z = std::min(gz(g0(tri_fault)),std::min(gz(g1(tri_fault)),gz(g2(tri_fault))));
	
	// determines if the triangle intersects the range defined by z_min and z_max.
	if ( std::min(max_z,z_max) < std::max(min_z,z_min) ) {
		return vector<line_segment>();
	}
	max_x = std::max(gx(g0(tri_fault)), std::max(gx(g1(tri_fault)), gx(g2(tri_fault))))-x_inf;
	min_x = std::min(gx(g0(tri_fault)), std::min(gx(g1(tri_fault)), gx(g2(tri_fault))))-x_inf;
	max_y = std::max(gy(g0(tri_fault)), std::max(gy(g1(tri_fault)), gy(g2(tri_fault))))-y_inf;
	min_y = std::min(gy(g0(tri_fault)), std::min(gy(g1(tri_fault)), gy(g2(tri_fault))))-y_inf;
	
	// defines the bouding box of the fault triangle.
	j_min = std::max(int(ceil(min_x/dx)),1); j_max = std::min(int(ceil(max_x/dx)),cols-1);
	i_min = std::max(int(ceil(min_y/dy)),1); i_max = std::min(int(ceil(max_y/dy)),rows-1);
	
	vector<line_segment> faults_intersection;
		
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
	return faults_intersection;
}


