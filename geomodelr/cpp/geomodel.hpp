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
#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <vector>
#include <tuple>
#include <utility>
#include <string>
#include <limits>
#include <map>
#include <memory>

using namespace boost;
using std::vector;
using std::string;
using std::map;

typedef geometry::model::point<double, 2, geometry::cs::cartesian> point2;
typedef geometry::model::point<double, 3, geometry::cs::cartesian> point3;
typedef geometry::model::box<point2> box;
typedef geometry::model::polygon<point2, false, false> polygon;
typedef polygon::ring_type ring;
typedef geometry::model::linestring<point2> line;
typedef std::pair<box, int> value;
typedef geometry::index::rtree<value, geometry::index::quadratic<16>> rtree;

template<class Point>
inline typename geometry::coordinate_type<Point>::type gx(const Point& p){
	return geometry::get<0, Point>(p);
}

template<class Point>
inline typename geometry::coordinate_type<Point>::type gy(const Point& p){
	return geometry::get<1, Point>(p);
}

template<class Point>
inline typename geometry::coordinate_type<Point>::type gz(const Point& p){
	return geometry::get<2, Point>(p);
}

typedef python::list pylist;
typedef python::dict pydict;
typedef python::object pyobject;
typedef python::tuple pytuple;

class Section;
class Model;

/* Checks if a unit is not empty or "NONE" */
struct ValidUnit {
	const Section * section;
	ValidUnit( const Section * section ); 
	bool operator() (const value& b) const;
};



class Match {
	map<int, vector<int>> a_to_b;
	map<int, vector<int>> b_to_a;
	vector<bool> a_free;
	vector<bool> b_free;
	friend class Model;
public:
	Match(const vector<std::pair<int, int>>& match, size_t sa, size_t sb);
	pylist get() const;
};

class Model {
	point2 base_point;
	point2 direction;
	vector<Section *> sections;
	vector<Match> match;
	vector<double> cuts;
	Match make_match( const Section& a, const Section& b );
	Match load_match( const pylist& match, size_t sa, size_t sb );
	
	std::pair<int, double> closest_match(bool a, int a_idx, int pol_idx, const point2& pt) const;
	struct Possible {
		Possible(int a_match, int b_match, double a_dist, double b_dist);
		int a_match;
		int b_match;
		double a_dist;
		double b_dist;
		bool operator<( const Possible& other ) const;
		double distance( double c ) const;
	};
	// Returns all the possible matches of this 2d point, given the distance is unknown.
	std::pair<point2, double> to_model_point(const point3& pt) const;
	point3 to_inverse_point(const point2& p, double cut) const;
	vector<Possible> get_candidates(size_t a_idx, const point2& pt) const;
	vector<Possible> all_closest(size_t a_idx, const point2& pt) const;
	std::tuple<int, int, double> closest_to( size_t a_idx, const point2& pt, double cut ) const;
public:
	Model(const pyobject& basepoint, const pyobject& direction, 
	      const pylist& sections);
	virtual ~Model();
	// Methods to create matches or load them from files.
	void make_matches();
	void set_matches(const pylist& matching);
	pylist get_matches() const;
	// Methods to query matches.
	pylist possible_closest(const pyobject& pt) const;
	pytuple model_point(const pyobject& pt) const;
	pytuple inverse_point(const pyobject& pt) const;
	pytuple closest(const pyobject& pt) const;
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
	
	template<typename Predicates>
	vector<std::pair<int, double>> closer_than( const point2& pt, double distance, const Predicates& predicates ) const {
		point2 mx(gx(pt) + distance, gy(pt) + distance);
		point2 mn(gx(pt) - distance, gy(pt) - distance);
		box bx(mn, mx);
		vector<std::pair<int, double>> ret;
		for ( auto it = this->polidx->qbegin( geometry::index::intersects(bx) and geometry::index::satisfies(this->valid) and predicates );
			it != this->polidx->qend(); it++ ) {
			// Check the actual distance to a polygon.
			int idx = it->second;
			double poldist = geometry::distance(this->polygons[idx], pt);
			if ( poldist <= distance ) {
				ret.push_back(std::make_pair(idx, poldist));
			}
		}
		return ret;
	}
	template<typename Predicates>
	std::pair<int, double> closest_to( const point2& p, const Predicates& predicates ) const {
		if ( this->polidx == nullptr )
			return std::make_pair(-1, std::numeric_limits<double>::infinity());
		
		double maxboxdist = 0.0;
		double mindist = std::numeric_limits<double>::infinity();
		
		int minidx = -1;
		int knear = 1;
		
		bool new_to_check;
		do {
			int n = 0; 
			new_to_check = false;
			for ( auto it = this->polidx->qbegin( geometry::index::nearest(p, knear) and
							      predicates );
				it != this->polidx->qend(); it++ ) {
				// Skip already checked.
				if ( n < knear/2 ) 
				{
					n++;
					continue;
				}
				// Check if new polygons where checked.
				
				new_to_check = true;
				
				// Check the maximum distance from the box to the point.
				// That distance is always lower than the distance to the polygon.
				double boxdist = geometry::distance(p, it->first);
				maxboxdist = std::max(boxdist, maxboxdist);
				
				// Then check the minimum actual distance to a polygon.
				int idx = it->second;
				double dist = geometry::distance(p, this->polygons[idx]);
				
				if ( dist < mindist ) {
					mindist = dist;
					minidx = idx;
				}
			}
			// Increase the number of knear.
			knear *= 2;
			// Do it until none was checked or we have checked boxes beyond the closest polygon.
		} while ( new_to_check && maxboxdist < mindist );
		
		return std::make_pair(minidx, mindist); 
	}
public:
	Section(double cut, const pylist& points, 
		const pylist& polygons, const pylist& units, 
		const pylist& lines, const pylist& lnames );
	virtual ~Section();
	pydict info() const;
	pytuple closest( const pyobject& pypt ) const;
};

