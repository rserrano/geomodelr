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
#include "speed_up.hpp"
#include "model.hpp"
#include <algorithm>    // std::min_element, std::max_elemen
#include <stdlib.h>
#include <iomanip>

std::pair<double, bool> find_unit_limits_cpp(const Model* model, double xp, double yp,
    double z_max, double z_min, double eps) {
  
    wstring unit_max = g0(model->closest(point3(xp,yp,z_max)));

    if ( unit_max.compare(g0(model->closest(point3(xp,yp,z_min))))==0 ){
        return std::make_pair( z_min, false );
    }    
    do {
        double z_mean = (z_max+z_min)/2.0;
        if ( unit_max.compare( g0(model->closest(point3(xp, yp, z_mean)))) == 0 ) {
            z_max = z_mean;
        }
        else {
            z_min = z_mean;
        }
    } while( (z_max-z_min)>eps );

    return std::make_pair( z_min, true );
}

line ray_intersection(const point2& pt, double R, const line_segment& poly_s){

    double x1 = gx(poly_s.first) - gx(pt); double y1 = gy(poly_s.first) - gy(pt);
    double x2 = gx(poly_s.second) - gx(pt); double y2 = gy(poly_s.second) - gy(pt);
    double dx = x2 - x1; double dy = y2 - y1;
    
    double A = dx*dx + dy*dy;
    double B = 2.0*(x1*dx + y1*dy);
    double C = x1*x1 + y1*y1 - R*R;
    double D = B*B - 4*A*C;

    //std::cerr << "D: " << std::setprecision(10) << D << std::endl;
    line output;
    output.push_back(pt);
    point2 int_pt;

    if (D<epsilon){
        double T = -0.5*B/A;
        //std::cerr << "T: " << std::setprecision(10) << T << std::endl;
        int_pt = point2(x1 + T*dx, y1 + T*dy);
    } else{

        double T1 = 0.5*(-B + std::sqrt(D))/A;
        double T2 = -(B/A + T1);
        //std::cerr << "T1 y T2: " <<std::setprecision(10) << T1 << "\t" << T2 << std::endl;
        if ((std::abs(T1)<epsilon) || (std::abs(T2)<epsilon)){
            //int_pt = poly_s.first; geometry::subtract_point(int_pt,pt);
            int_pt = point2(x1,y1);
        }
        else if ((std::abs(T2-1.0)<epsilon) || (std::abs(T1-1.0)<epsilon)){
            //int_pt = poly_s.second; geometry::subtract_point(int_pt,pt);
            int_pt = point2(x2,y2);
        }
        //else if (std::abs(T1-T2)<epsilon){
        else{
            double T = 0.5*(T1+T2);
            //output.push_back(point2(x1 + T*dx + gx(pt), y1 + T*dy+gy(pt)));
            int_pt = point2(x1 + T*dx, y1 + T*dy);
        }
    }
    geometry::multiply_value(int_pt,(R + epsilon)/R);
    geometry::add_point(int_pt,pt);
    output.push_back(int_pt);
    return output;
}


double distance_poly_fault_pt(const point2& pt, const polygon& poly,const rtree_seg* poly_seg_tree,
    const vector<rtree_seg *>& fault_lines){

    double poly_dist = geometry::distance(poly, pt);

    if ((poly_dist==0) or (fault_lines.size()==0)){
        return poly_dist;
    }
    std::vector<line_segment> returned_values;
    poly_seg_tree->query(geometry::index::nearest(pt,1), std::back_inserter(returned_values));
    line ray = ray_intersection(pt, poly_dist, returned_values[0]);

    for (auto fault_tree: fault_lines){
        for ( auto it = fault_tree->qbegin( geometry::index::intersects(ray)); it != fault_tree->qend(); it++ ) {
            return std::numeric_limits<double>::infinity();
        }
    }
    return poly_dist;
}

void clear_vector(double y0,vector<value_s>& results){

    int c=0;
    while (c<results.size()){
        line_segment edge = g0(results[c]);
        double y_max = std::max(gy(edge.first),gy(edge.second));        
        if (std::abs(y0-y_max)<epsilon){
            results.erase(results.begin()+c);
            c--;            
        }
        c++;
    }
}

std::pair<line_segment,double> cross_segment(const point2& pt, const line_segment& poly_edge){

    std::cerr << "Punto 3.1" << std::endl;
    std::cerr << geometry::wkt(pt) << std::endl;
    std::cerr << geometry::wkt(poly_edge) << std::endl;
    double x0 = gx(poly_edge.first); double y0 = gy(poly_edge.first);
    double xf = gx(poly_edge.second); double yf = gy(poly_edge.second);
    double xp = gx(pt); double yp = gy(pt);
    std::cerr << "Punto 3.2" << std::endl;
    double dx = xf-x0; double dy = yf-y0;
    double t = ((xp-x0)*dx + (yp-y0)*dy)/(dx*dx + dy*dy);
    std::cerr << "Punto 3.3" << std::endl;
    point2 output;
    if (t>=1.0){
        output = point2(xf - xp, yf - yp);
    } else if(t<=0.0){
        output = point2(x0 - xp, y0 - yp);
    } else{
        output = point2(x0 + t*dx - xp, y0 +t *dy - yp);
    }
    std::cerr << "Punto 3.4" << std::endl;
    double norm = std::sqrt(geometry::dot_product(output,output));
    std::cerr << "Punto 3.5" << std::endl;
    geometry::multiply_value(output,(norm + epsilon)/norm);
    geometry::add_point(output,pt);
    std::cerr << "Punto 3.6" << std::endl;
    return std::make_pair(line_segment(pt,output),norm);

}

template<typename Predicates>
std::pair<line_segment,double> distance_poly_pt(int idx, const line_segment& rigth_ray,
    const rtree_s* poly_seg_tree, const Predicates& predicates){

    vector<value_s> results;
    poly_seg_tree->query(geometry::index::intersects(rigth_ray) and geometry::index::satisfies(predicates),
        std::back_inserter(results));

    std::cerr << "Right ray: " << geometry::wkt(rigth_ray) << std::endl;
    std::cerr << "size segment vec NO: " << results.size() << std::endl;
    if (results.size()>0){
        clear_vector(gy(rigth_ray.first),results);
    }
    std::cerr << "size segment vec YES: " << results.size() << std::endl;

    if (results.size()%2>0){
        return std::make_pair(line_segment(),0.0);
    }
    else{
        
        results.clear();
        poly_seg_tree->query(geometry::index::satisfies(predicates),std::back_inserter(results));
        if (results.size()>0){
            std::cerr << idx << ", " << g1(results[0]) << ", " << g2(results[0]) << std::endl;
        }

        results.clear();
        poly_seg_tree->query(geometry::index::nearest(rigth_ray.first,1) and geometry::index::satisfies(predicates),
            std::back_inserter(results));
        std::cerr << results.size() << std::endl;
        return cross_segment(rigth_ray.first, g0(results[0]));
    }
}

double distance_poly_fault_pt2(int idx, const point2& pt, double x_rigth,const rtree_s* poly_seg_tree,
    const vector<rtree_seg *>& fault_lines){

    auto just_polidx = [idx](const value_s& seg) {
        return (g1(seg)==idx) || (g2(seg)==idx);
    };

    std::pair<line_segment,double> ray_distance = distance_poly_pt(idx, line_segment(pt,point2(x_rigth,gy(pt))),
        poly_seg_tree, just_polidx);

    if (ray_distance.second==0.0){
        return 0.0;
    } else{

        line_segment ray = ray_distance.first;
        for (auto fault_tree: fault_lines){
            for ( auto it = fault_tree->qbegin( geometry::index::intersects(ray)); it != fault_tree->qend(); it++ ) {
                return std::numeric_limits<double>::infinity();
            }
        }
    return ray_distance.second;
    }

}
