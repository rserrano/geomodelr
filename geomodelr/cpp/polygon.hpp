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
#include <functional>

class Match;
class Model;
class SectionPython;
class Section;

class Polygon {

	friend Section;
	friend Match;
	friend Model;
	friend SectionPython;

protected:
	double x_corner;
	polygon boost_poly;
	box bbox;
	const Section * section;
	rtree_seg * poly_lines;
	std::function<double(const point2&)> distance_point;
	
public:
	Polygon();
	Polygon(const polygon& poly, const box& bbox, const Section * section);
	virtual ~Polygon();
	
	std::pair<vector<line_segment>,double> ray_distance(const point2& pt) const;
	double ray_crossing ( const point2& pt, const point2& nd ) const;
	double distance_point_basic_faults(const point2& pt) const;
	double distance_point_cover_faults(const point2& pt) const;
	void set_distance_function( const wstring& s );
};

#endif
