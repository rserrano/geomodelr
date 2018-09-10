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

#ifndef GEOMODELR_BASIC_HPP
#define GEOMODELR_BASIC_HPP

#define CGAL_MESH_2_OPTIMIZER_VERBOSE

#include <boost/python.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/lambda/lambda.hpp>
#include <openvdb/openvdb.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <utility>
#include <string>
#include <limits>
#include <map>
#include <memory>
#include <limits>
#include <assert.h>
#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>

using namespace boost;
using std::vector;
using std::string;
using std::wstring;
using std::map;

wstring human_failure_type( const geometry::validity_failure_type& fail );

static const double tolerance = 1e-15;
static const double epsilon = 1e-5;
static const double boost_tol = std::numeric_limits<double>::epsilon();

typedef CGAL::Exact_predicates_inexact_constructions_kernel 	 Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, Kernel>    Vb;
typedef CGAL::Delaunay_mesh_vertex_base_2<Kernel, Vb>           	 Vbi;
typedef CGAL::Delaunay_mesh_face_base_2<Kernel>             	 Fb;
typedef CGAL::Triangulation_data_structure_2<Vbi, Fb>        	 Tds;
typedef CGAL::Exact_predicates_tag 							 Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds, Itag>  CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            	 Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              	 Mesher;
typedef CGAL::Polygon_2<Kernel> polygonCGAL;

typedef std::pair<vector<openvdb::Vec3f>,vector<openvdb::Vec3I>>   triangMesh;
typedef std::pair<vector<CDT::Point>,vector<openvdb::Vec3I>> 	   triangMesh2D;
typedef std::vector< std::vector< double > >					   vectorLayers;
typedef std::map<wstring, vector<std::pair<size_t, size_t>>>	   pairConstraints;
typedef std::tuple<triangMesh2D, vectorLayers, pairConstraints >   feflowInfo;

typedef geometry::model::point<double, 2, geometry::cs::cartesian> point2;
typedef geometry::model::point<double, 3, geometry::cs::cartesian> point3;
typedef geometry::model::segment<point2> line_segment;
typedef geometry::model::box<point2> box;
typedef geometry::model::polygon<point2, false, false> polygon;
typedef geometry::model::multi_polygon<polygon> multi_polygon;
typedef geometry::model::segment<point2> segment;

typedef polygon::ring_type ring;
typedef geometry::model::linestring<point2> line;
typedef geometry::model::linestring<point3> line_3d;
typedef geometry::model::multi_linestring<line> multi_line;

typedef std::tuple<std::tuple<double, double, double>, 
		   std::tuple<double, double, double>> bbox3;

typedef std::tuple<std::tuple<double, double>, 
		   std::tuple<double, double>> bbox2;

typedef std::pair<int, bool> line_anchor;

//typedef std::pair<box, int> value;
typedef std::tuple<box, wstring, int> value_f;
typedef std::tuple<line_segment, int> value_l; // Value to search for surface lines, fault intersection.
typedef std::tuple<segment, wstring> value_s;

typedef geometry::index::rtree<value_f, geometry::index::quadratic<16>> rtree_f;
typedef geometry::index::rtree<value_l, geometry::index::quadratic<16>> rtree_l; // Tree to search for surface line.
typedef geometry::index::rtree<value_s, geometry::index::quadratic<16>> rtree_s;
typedef geometry::index::rtree<line_segment, geometry::index::quadratic<16>> rtree_seg;

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

struct GeomodelrException : std::runtime_error
{
	GeomodelrException(const string&);
};

bool always_true( const value_f& v );

vector<size_t> sort_indexes(const vector<double> &v);

std::tuple<point2, double> point_line_projection( const point2& p, const line& l );

// Is geomodelr verbose.
extern bool geomodelr_verbose;

template<typename FA, typename FB>
class add_funcs {
public:
	FA fa;
	FB fb;
	add_funcs( const FA& fa, const FB& fb ) : fa(fa), fb(fb) {
	}
	
	bool operator() ( const value_f& v ) const {
		return (this->fa(v) and this->fb(v));
	}

};

#endif

