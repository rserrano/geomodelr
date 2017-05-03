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
	return closest;
}

bool always_true( const value_f& v )
{
	return true;
}
