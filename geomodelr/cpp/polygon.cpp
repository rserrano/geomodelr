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

Polygon::~Polygon()
{
	if ( this->poly_lines != nullptr ) {
		delete this->poly_lines;
	}
}

Polygon::Polygon(): x_corner(std::numeric_limits<double>::infinity()), poly_lines(nullptr) 
{	
}

Polygon::Polygon(double x_corner, const vector<line_segment>& poly_segments):x_corner(x_corner) {
	this->poly_lines = new rtree_seg(poly_segments);
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
	this->x_corner = env.max_corner().get<0>() + epsilon;
	this->poly_lines = new rtree_seg( poly_to_vecsegs(poly) );
}

std::pair<line_segment,double> cross_segment2(const point2& pt, const line_segment& poly_edge){

    double x0 = gx(poly_edge.first); double y0 = gy(poly_edge.first);
    double xf = gx(poly_edge.second); double yf = gy(poly_edge.second);
    double xp = gx(pt); double yp = gy(pt);

    double dx = xf-x0; double dy = yf-y0;
    double t = ((xp-x0)*dx + (yp-y0)*dy)/(dx*dx + dy*dy);

    point2 output;
    if (t>=1.0){
        output = point2(xf - xp, yf - yp);
    } else if(t<=0.0){
        output = point2(x0 - xp, y0 - yp);
    } else{
        output = point2(x0 + t*dx - xp, y0 + t*dy - yp);
    }
    double norm = std::sqrt(geometry::dot_product(output,output));
    geometry::multiply_value(output,(norm + epsilon)/norm);
    geometry::add_point(output,pt);

    return std::make_pair(line_segment(pt,output),norm);

}

void clear_vector(double y0, vector<line_segment>& results){
    int c=0;
    while (c<results.size()){
        //line_segment edge = g0(results[c]);
        double y_max = std::max(gy(results[c].first),gy(results[c].second));        
        if (std::abs(y0-y_max)<tolerance){
            results.erase(results.begin()+c);
            c--;            
        }
        c++;
    }
}

std::pair<line_segment,double> Polygon::ray_distance(const point2& pt) const{

	line_segment rigth_ray = line_segment(pt,point2(this->x_corner,gy(pt)));
    vector<line_segment> results;
    this->poly_lines->query(geometry::index::intersects(rigth_ray) , std::back_inserter(results));

    if (results.size()>0){
        clear_vector(gy(pt),results);
    }
    if (results.size()%2>0){
        return std::make_pair(line_segment(),0.0);
    }
    else{
        results.clear();
        this->poly_lines->query(geometry::index::nearest(pt,1), std::back_inserter(results));
        return cross_segment2(pt, results[0]);
    }
}

double Polygon::distance_point( const point2& pt ,const vector<rtree_seg *>& fault_lines) const {

    std::pair<line_segment,double> ray_dist_pair = this->ray_distance(pt);
    if (ray_dist_pair.second==0.0){
        return 0.0;

    } else{
        line_segment ray = ray_dist_pair.first;
        for (auto fault_tree: fault_lines){
            for ( auto it = fault_tree->qbegin( geometry::index::intersects(ray)); it != fault_tree->qend(); it++ ) {
                return std::numeric_limits<double>::infinity();
                //return ray_dist_pair.second;
            }
        }
    return ray_dist_pair.second;
    }
}

void Polygon::info( const point2& pt){
	
	line_segment rigth_ray = line_segment(pt,point2(this->x_corner,gy(pt)));
    vector<line_segment> results;
    this->poly_lines->query(geometry::index::intersects(rigth_ray) , std::back_inserter(results));

    std::cerr << std::setprecision(15) <<geometry::wkt(rigth_ray) << std::endl;
    std::cerr << results.size() << std::endl;
    clear_vector(gy(pt),results);
    std::cerr << results.size() << std::endl;
}


PolygonPython::PolygonPython(const pylist& points,const pylist& polygons){

	polygon pol;
	ring& outer = pol.outer();
	
	// Start filling the first ring.
	vector<line_segment> poly_segments;
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
		this->x_corner = env.max_corner().get<0>() + epsilon;
		this->poly_lines = new rtree_seg( poly_to_vecsegs(pol) );
	}
}

double PolygonPython::distance_poly_test(const pylist& pt) const{

	double x = python::extract<double>(pt[0]);
	double y = python::extract<double>(pt[1]);

	vector<rtree_seg *> fault_lines(0);
	return this->distance_point(point2(x,y),fault_lines);
}