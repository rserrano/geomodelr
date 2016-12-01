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
using std::wstring;
using std::map;

typedef geometry::model::point<double, 2, geometry::cs::cartesian> point2;
typedef geometry::model::point<double, 3, geometry::cs::cartesian> point3;
typedef geometry::model::box<point2> box;
typedef geometry::model::polygon<point2, false, false> polygon;
typedef polygon::ring_type ring;
typedef geometry::model::linestring<point2> line;
typedef std::pair<box, int> value;
typedef std::tuple<box, wstring, int> value_f;
typedef geometry::index::rtree<value, geometry::index::quadratic<16>> rtree;
typedef geometry::index::rtree<value_f, geometry::index::quadratic<16>> rtree_f;
typedef std::tuple<int, int, int> triangle;
typedef std::tuple<int, int> edge;

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

typedef std::tuple<point3, point3, point3> triangle_pt;
template<class T>
inline const typename std::tuple_element<0, T>::type& g0(const T& t){
	return std::get<0>(t);
}
template<class T>
inline const typename std::tuple_element<1, T>::type& g1(const T& t){
	return std::get<1>(t);
}
template<class T>
inline const typename std::tuple_element<2, T>::type& g2(const T& t){
	return std::get<2>(t);
}

template<class T>
inline typename std::tuple_element<0, T>::type& g0(T& t){
	return std::get<0>(t);
}
template<class T>
inline typename std::tuple_element<1, T>::type& g1(T& t){
	return std::get<1>(t);
}
template<class T>
inline typename std::tuple_element<2, T>::type& g2(T& t){
	return std::get<2>(t);
}

class Section;
class Model;
class Topography;

struct GeomodelrException : std::runtime_error
{
	GeomodelrException(const string&);
};

bool always_true( const value& v ); 

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
	static pyobject pytriangulate;
	const Section * a;
	const Section * b;
	map<int, vector<int>> a_to_b;
	map<int, vector<int>> b_to_a;
	map<wstring, vector<AlignedTriangle>> faults;
	rtree_f * faultidx;
	void set( const vector<std::pair<int, int>>& match );
public:
	Match( const Section * a, const Section * b );
	virtual ~Match();
	void match_polygons();
	map<wstring, vector<triangle_pt>> match_lines();
	void load_polygons_match( const pylist& match );
	pylist get() const;
	static void load_triangulate();
	static vector<triangle> triangulate(const vector<point3>& pa, const vector<point3>& pb);
	int crosses_triangles(const point2& pt, double cut) const;
};

class Model {
	point2 base_point;
	point2 direction;
	vector<Section *> sections;
	vector<Match *> match;
	vector<double> cuts;
	Topography * topography;
	
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
	void clear_matches();
	// Returns all the possible matches of this 2d point, given the distance is unknown.
	std::pair<point2, double> to_model_point(const point3& pt) const;
	point3 to_inverse_point(const point2& p, double cut) const;
	vector<Possible> get_candidates(size_t a_idx, const point2& pt) const;
	vector<Possible> all_closest(size_t a_idx, const point2& pt) const;
	std::tuple<int, int, double> closest_to( size_t a_idx, const point2& pt, double cut ) const;
public:
	Model(const pyobject& basepoint, const pyobject& direction, 
	      const pyobject& map, const pyobject& topography,
	      const pylist& sections);
	virtual ~Model();
	// Methods to create matches or load them from files.
	pydict make_matches(); // Returns the faults in global coordinates, (at least until moving plane-fault intersection to C++).
	void set_matches(const pylist& matching);
	pylist get_matches() const;
	// Methods to query matches.
	pylist possible_closest(const pyobject& pt) const;
	pytuple model_point(const pyobject& pt) const;
	pytuple inverse_point(const pyobject& pt) const;
	pytuple closest(const pyobject& pt) const;
	pytuple closest_topo(const pyobject& pt) const;
	pydict info() const;
	double height(const pyobject& pt) const;
};

class Topography {
	point2 point;
	point2 sample;
	size_t dims[2];
	vector<double> heights;
public:
	Topography( const pyobject& point, const pyobject& sample, const pyobject& dims, const pylist& heights );
	double height(const point2&) const;
};

/* C++ section that queries points to polygons so much faster. */
class Section {
	friend Match;
	friend Model;
	
	wstring name;
	double cut;
	vector<polygon> polygons;
	vector<wstring> units;
	vector<line> lines;
	vector<wstring> lnames;
	rtree * polidx; // To be initialized after polygons and lines.

	template<typename Predicates>
	vector<std::pair<int, double>> closer_than( const point2& pt, double distance, const Predicates& predicates ) const {
		point2 mx(gx(pt) + distance, gy(pt) + distance);
		point2 mn(gx(pt) - distance, gy(pt) - distance);
		box bx(mn, mx);
		vector<std::pair<int, double>> ret;
		for ( auto it = this->polidx->qbegin( geometry::index::intersects(bx) and predicates );
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
	Section(const wstring& name, double cut, const pylist& points, 
		const pylist& polygons, const pylist& units, 
		const pylist& lines, const pylist& lnames );
	virtual ~Section();
	pydict info() const;
	pytuple closest( const pyobject& pypt ) const;
	
};

pylist test_faultplane_for_lines(const pylist& pyla, const pylist& pylb);

