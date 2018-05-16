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
#include <ctime>

Polygon::~Polygon()
{
	if ( this->poly_lines != nullptr ) {
		delete this->poly_lines;
	}
}

Polygon::Polygon(): x_corner(std::numeric_limits<double>::infinity()), poly_lines(nullptr) 
{	
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

Polygon::Polygon(const polygon& poly){
	box env;
	geometry::envelope(poly, env);
	this->b_box = env;
	this->x_corner = env.max_corner().get<0>() + epsilon;
	this->poly_lines = new rtree_seg( poly_to_vecsegs(poly) );
	this->boost_poly = poly;
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
    // double norm = std::sqrt(geometry::dot_product(output,output));
    // geometry::multiply_value(output,(norm + epsilon)/norm);
    // geometry::add_point(output,pt);
    double norm = std::sqrt(xf*xf + yf*yf);
    double Re = (norm + epsilon)/norm;

    return std::make_pair(line_segment(pt,point2(xf*Re + xp, yf*Re + yp)),norm);

}

std::pair<line_segment,double> Polygon::ray_distance(const point2& pt) const{

    vector<line_segment> results;   
    this->poly_lines->query(geometry::index::nearest(pt,1), std::back_inserter(results));
    return cross_segment(pt, results.front());
}

double Polygon::distance_point(const point2& pt, const rtree_l* fault_lines) const {

	/* Checks if pt is inside the bounding square of the polygon.
	   After that, checks if pt is inside the polygon.*/
    if (geometry::within(pt,this->b_box) && geometry::covered_by(pt,this->boost_poly)){
        return 0.0;
    } else{

    	std::pair<line_segment,double> ray_dist_pair = this->ray_distance(pt);
        for ( auto it = fault_lines->qbegin( geometry::index::intersects(ray_dist_pair.first));
        	it != fault_lines->qend(); it++ ) {
            return std::numeric_limits<double>::infinity();
        }
        return ray_dist_pair.second;
    }
}


double Polygon::distance_point_new(const point2& pt, const rtree_l* fault_lines) const {

	if (geometry::covered_by(pt,this->boost_poly)){
		return 0.0;
	}

	// Start analyze the first ring.
	const ring& outer = this->boost_poly.outer();
	size_t nnodes = outer.size();
	bool not_intersection;
	double min_dist = std::numeric_limits<double>::infinity();

	for ( size_t k = 0; k < nnodes; k++ ) {
		
		line_segment edge = line_segment(outer[k],outer[(k+1)%nnodes]);
		auto ray_pair = cross_segment(pt,edge);
		not_intersection = true;
		std::cerr << "Segment: " << k << ". Distance: " ray_pair.second << std::endl;;
		if (ray_pair.second<min_dist){
			for ( auto it = fault_lines->qbegin( geometry::index::intersects(ray_pair.first));
	    		it != fault_lines->qend(); it++ ) {
	    	 	not_intersection = false;
	    	 	std::cerr << "intersects: " << g1(*it) <<std::endl;
	        	break;
	        }
	        if (not_intersection){
	        	min_dist = ray_pair.second;
	        }
	    }
	}

	for (auto& inner: geometry::interior_rings(this->boost_poly)){		
		
		nnodes = inner.size();
		for ( size_t k = 0; k < nnodes; k++ ) {
			
			line_segment edge = line_segment(inner[k],inner[(k+1)%nnodes]);
			auto ray_pair = cross_segment(pt,edge);
			not_intersection = true;
			if (ray_pair.second<min_dist){
				for ( auto it = fault_lines->qbegin( geometry::index::intersects(ray_pair.first));
		    		it != fault_lines->qend(); it++ ) {
		    	 	not_intersection = false;
		        	break;
		        }
		        if (not_intersection){
		        	min_dist = ray_pair.second;
		        }
		    }
		}
	}
	return min_dist;
}


PolygonPython::PolygonPython(const pylist& points,const pylist& polygons){

	polygon pol;
	ring& outer = pol.outer();
	
	// Start filling the first ring.
	size_t nnodes = python::len(polygons[0]);

	for ( size_t k = 0; k < nnodes; k++ ) {
		pylist pypt = pylist(points[polygons[0][k]]);
		outer.push_back(point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])));
	}

	//this->num_segs = nnodes;
	// Then fill the rest of the rings.
	size_t nrings = python::len(polygons);
	if ( nrings > 1 ) { 
		pol.inners().resize(nrings-1);
		for ( size_t j = 1; j < nrings; j++ ) {
			ring& inner = pol.inners()[j-1];// jth ring.
			size_t nnodes = python::len(polygons[j]);
			//this->num_segs += nnodes;
			for ( size_t k = 0; k < nnodes; k++ ) {
				pylist pypt = pylist(points[polygons[j][k]]);
				inner.push_back(point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])));
			}
		}
	}
	
	geometry::validity_failure_type failure;
	if ( not geometry::is_valid(pol, failure) ) {
		geometry::correct(pol);
		string reason;
		if ( not geometry::is_valid(pol, reason) ) {
			if ( geomodelr_verbose ) {
				std::wcerr << L"non valid polygon in section" << "\n";
			}
			// continue; not avoiding non valid polygons, as they have been validated by shapely. Somehow these polygons get wronget.
		}
	}

	if (nrings>0){
		box env;
		geometry::envelope(pol, env);
		this->b_box = env;
		this->x_corner = env.max_corner().get<0>() + epsilon;
		this->poly_lines = new rtree_seg( poly_to_vecsegs(pol) );
		this->boost_poly = pol;
	}
}

double PolygonPython::distance_poly_test(const pylist& pt) const{

	double x = python::extract<double>(pt[0]);
	double y = python::extract<double>(pt[1]);

	rtree_l * fault_lines = new rtree_l(vector<value_l>(0));
	return this->distance_point_new(point2(x,y),fault_lines);
}

pytuple PolygonPython::time_poly_test(const pylist& pt,int N) const{

	double x = python::extract<double>(pt[0]);
	double y = python::extract<double>(pt[1]);

	rtree_l * fault_lines = new rtree_l(vector<value_l>(0));

	point2 PT = point2(x,y);

	// Boost distance
	clock_t  begin = clock();
	for (int k=0; k<N; k++){
		double aux = geometry::distance(PT,this->boost_poly);
	}
	clock_t end = clock();
	double elapsed_2 = double(end - begin) / CLOCKS_PER_SEC;

		// Own distance
	begin = clock();
	for (int k=0; k<N; k++){
		double aux = this->distance_point(PT,fault_lines);
	}
	end = clock();
	double elapsed_1 = double(end - begin) / CLOCKS_PER_SEC;

	return python::make_tuple(elapsed_1,elapsed_2);
}