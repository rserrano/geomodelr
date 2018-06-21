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
#ifndef GEOMODELR_MATCH_HPP
#define GEOMODELR_MATCH_HPP

#include "basic.hpp"
#include "section.hpp"

// Included for friend of match.
class Model;

class AlignedTriangle {
	friend class Match;
	point3 normal;
	point3 point;
	polygon triangle;
public:
	AlignedTriangle(const std::tuple<point3, point3, point3>& triangle);
	int crosses_triangle(const point2& point, double cut) const;
};

class Match {
	friend class Model;
protected:
	const Section * a;
	const Section * b;
	map<int, vector<int>> a_to_b;
	map<int, vector<int>> b_to_a;
	map<wstring, vector<AlignedTriangle>> faults; // The aligned triangles to test.
	map<wstring, std::tuple<int, int>> rel_faults; // Which fault is related to which other fault.
	map<wstring, multi_polygon> excluded_area;
	rtree_f * faultidx;
	const map<wstring, wstring> * params;
	bool faults_disabled;
public:
	Match( const Section * a, const Section * b );
	virtual ~Match();
	void match_polygons();
	std::tuple<map<wstring, vector<triangle_pt>>, map<wstring, vector<size_t>>> match_lines( const map<wstring, wstring>& feature_types );
	std::tuple<int, int, int> crosses_triangles(const point2& pt, double cut) const;
	void set( const vector<std::pair<int, int>>& match );
	void set_params( const map<wstring, wstring> * params );
};

class MatchPython : public Match {
public:
	pylist get() const;
	void set( const pylist& match );
};

pylist test_faultplane_for_lines(const pylist& pyla, const pylist& pylb);
#endif
