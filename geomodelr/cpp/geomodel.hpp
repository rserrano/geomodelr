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

#include <boost/python.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <iostream>
#include <vector>
#include <tuple>
#include <utility>
#include <string>
#include <limits>
#include <map>

using namespace boost;
using std::vector;
using std::string;
using std::map;

typedef geometry::model::point<double, 2, geometry::cs::cartesian> point2;
typedef geometry::model::box<point2> box;
typedef geometry::model::polygon<point2, false, false> polygon;
typedef polygon::ring_type ring;
typedef geometry::model::linestring<point2> line;
typedef std::pair<box, int> value;
typedef geometry::index::rtree<value, geometry::index::quadratic<16>> rtree;

typedef python::list pylist;
typedef python::dict pydict;
typedef python::object pyobject;
typedef python::tuple pytuple;

class Section;
class Model;

/* Checks if a unit is not empty or "NONE" */
struct ValidUnit {

	const Section * section;
	ValidUnit( const Section * section ):section(section) 
	{
		
	}
	bool operator() (const value& b) const;
};

class Match {
	map<int, vector<int>> a_to_b;
	map<int, vector<int>> b_to_a;
	vector<bool> a_free;
	vector<bool> b_free;
public:
	Match(const vector<std::pair<int, int>>& match, size_t sa, size_t sb);
};

class Model {
	double base_point[2];
	double base_line[2];
	vector<Section> sections;
	vector<Match> match;
	Match make_match( const Section& a, const Section& b );
	Match load_match( const pylist& match, size_t sa, size_t sb );
public:
	Model(const pylist& basepoint, const pylist& direction, 
	      const pylist& sections);
	void make_matches();
	void load_matches(const pylist& matching);
};

/* C++ section that queries points to polygons so much faster. */
class Section {
	double cut;
	vector<polygon> polygons;
	vector<string> units;
	vector<line> lines;
	vector<string> lnames;
	rtree * polidx; // To be initialized after polygons and lines.
	
	ValidUnit valid;
	friend class ValidUnit;
	friend Model;
public:
	Section(double cut, const pylist& points, 
		const pylist& polygons, const pylist& units, 
		const pylist& lines, const pylist& lnames );	
	pydict info() const;
	pytuple closest( const pyobject& pypt ) const;
};

