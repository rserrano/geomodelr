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
