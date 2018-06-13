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

#include "polygon.hpp"
#include "section.hpp"
#include <ctime>
#include <cmath>

Polygon::~Polygon()
{
	if ( this->poly_lines != nullptr ) {
		delete this->poly_lines;
	}
}

Polygon::Polygon(): x_corner(std::numeric_limits<double>::infinity()), poly_lines(nullptr)
{
	this->distance_point = std::bind(&Polygon::distance_point_basic_faults, this, std::placeholders::_1);
}

double Polygon::distance_point_basic_faults(const point2& pt) const {
	return geometry::distance(pt, this->boost_poly);
}

double Polygon::distance_point_cover_faults(const point2& pt) const {
	/* Checks if pt is inside the bounding square of the polygon.
	   After that, checks if pt is inside the polygon.*/
	if (geometry::within(pt,this->bbox) && geometry::covered_by(pt,this->boost_poly)){
		return 0.0;
	} else {
		std::pair<line_segment,double> ray_dist_pair = this->ray_distance(pt);
		for ( auto it = this->section->fault_lines->qbegin( geometry::index::intersects(ray_dist_pair.first));
			   it != this->section->fault_lines->qend(); it++ ) {
			const point2& ps = this->section->line_ends[g1(*it)].first;
			const point2& pe = this->section->line_ends[g1(*it)].second;
			return std::min( this->ray_crossing( pt, ps ), this->ray_crossing(pt, pe) );
		}
		return ray_dist_pair.second;
	}
}
	
vector<line_segment> poly_to_vecsegs(const polygon& poly){
	// Exterior ring
	int Nnodes = poly.outer().size();
	vector<line_segment> poly_segments;
	for (int k=0; k<Nnodes; k++){
		poly_segments.push_back(line_segment( poly.outer()[k],poly.outer()[(k+1)%Nnodes]) );
	}

	// Interior rings
	for (auto& hole: poly.inners()){
		Nnodes = hole.size();
		for (int k=0; k<Nnodes; k++){
			poly_segments.push_back(line_segment( hole[k],hole[(k+1)%Nnodes]) );
		}
	}

	return poly_segments;
}

Polygon::Polygon( const polygon& poly, const box& bbox, const Section * section ):boost_poly(poly), bbox(bbox), section(section), poly_lines(nullptr) {
	this->x_corner = bbox.max_corner().get<0>() + epsilon;
	this->poly_lines = new rtree_seg( poly_to_vecsegs(poly) );
	this->boost_poly = poly;
	this->distance_point = std::bind(&Polygon::distance_point_basic_faults, this, std::placeholders::_1);
}

std::pair<line_segment,double> cross_segment(const point2& pt, const line_segment& poly_edge){
	
	double x0 = gx(poly_edge.first); double y0 = gy(poly_edge.first);
	double xf = gx(poly_edge.second); double yf = gy(poly_edge.second);
	double xp = gx(pt); double yp = gy(pt);
	
	double dx = xf-x0; double dy = yf-y0;
	double t = ((xp-x0)*dx + (yp-y0)*dy)/(dx*dx + dy*dy);
	
	point2 output;
	if (t>=1.0){
		//output = point2(xf - xp, yf - yp);
		xf -= xp; yf -= yp;
	} else if(t<=0.0){
		//output = point2(x0 - xp, y0 - yp);
		xf = x0 - xp; yf = y0 - yp;
	} else{
		//output = point2(x0 + t*dx - xp, y0 + t*dy - yp);
		xf = x0 + t*dx - xp; yf = y0 + t*dy -yp;
	}
	double norm = std::sqrt(xf*xf + yf*yf);
	double Re = (norm + epsilon)/norm;

	return std::make_pair(line_segment(pt,point2(xf*Re + xp, yf*Re + yp)),norm);

}

std::pair<line_segment,double> Polygon::ray_distance(const point2& pt) const{

	vector<line_segment> results;   
	this->poly_lines->query(geometry::index::nearest(pt,1), std::back_inserter(results));
	return cross_segment(pt, results.front());
}

double Polygon::ray_crossing ( const point2& pt, const point2& nd ) const {
	// Check if the bounding box can hit the ray.
	point2 v = nd;
	geometry::subtract_point(v, pt);
	double xl, yl;
	if ( gx(v) > 0 ) {
		xl = this->bbox.max_corner().get<0>();
		if ( gx(nd) > xl ) {
			return std::numeric_limits<double>::infinity();
		}
	} else {
		xl = this->bbox.min_corner().get<0>();
		if ( gx(nd) <= xl ) {
			return std::numeric_limits<double>::infinity();
		}
	}
	if ( gy(v) > 0 ) {
		yl = this->bbox.max_corner().get<1>();
		if ( gy(nd) > yl ) {
			return std::numeric_limits<double>::infinity();
		}
	} else {
		yl = this->bbox.min_corner().get<1>();
		if ( gy(nd) <= yl ) {
			return std::numeric_limits<double>::infinity();
		}
	}
	
	double minm = std::numeric_limits<double>::infinity();
	
	minm = std::min((xl-gx(nd))/gx(v), minm);
	minm = std::min((yl-gy(nd))/gy(v), minm);
	
	if ( std::isinf( minm ) ) {
		return std::numeric_limits<double>::infinity();
	}
	geometry::multiply_value(v, minm);
	geometry::add_point(v, nd);
	line_segment s(nd, v);
	double mind = std::numeric_limits<double>::infinity();
	for ( auto it = this->poly_lines->qbegin( geometry::index::intersects(s) );
	      it != this->poly_lines->qend(); it++ ) {
		vector<point2> intp;
		geometry::intersection( *it, s, intp );
		mind = std::min( mind, geometry::distance(pt, intp[0]) );
	}
	return mind;
}

double Polygon::set_distance_function( const wstring& s ) {
	if ( s == L"cover" ) {
		this->distance_point = std::bind(&Polygon::distance_point_cover_faults, this, std::placeholders::_1);
	} else {
		this->distance_point = std::bind(&Polygon::distance_point_basic_faults, this, std::placeholders::_1);
	}
}

