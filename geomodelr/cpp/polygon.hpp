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
#ifndef GEOMODELR_Polygon_HPP
#define GEOMODELR_Polygon_HPP
#include "basic.hpp"

class Section;
class SectionPython;


class Polygon {

	friend Section;
	friend SectionPython;

protected:
	double x_corner;
	rtree_seg * poly_lines;
public:
	Polygon();
	Polygon(double x_corner, const vector<line_segment>& poly_segments);
	Polygon(const polygon& poly);
	virtual ~Polygon();

	std::pair<line_segment,double> ray_distance(const point2& pt) const;
	double distance_point( const point2& pt ,const vector<rtree_seg *>& fault_lines) const;
	void info( const point2& pt);
};

class PolygonPython : public Polygon {

public:
	PolygonPython(const pylist& points,const pylist& polygon);
 	double distance_poly_test(const pylist& pt) const;
};
#endif