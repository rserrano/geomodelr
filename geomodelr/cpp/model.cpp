//  Copyright Joel de Guzman 2002-2004. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//  Hello World Example from the tutorial
//  [Joel de Guzman 10/9/2002]

#include <boost/python.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <iostream>
#include <vector>
#include <tuple>
#include <string>
using namespace boost;
using std::vector;
using std::string;

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

class Section {
	double cut;
	vector<polygon> polygons;
	vector<string> units;
	vector<line> lines;
	vector<string> lnames;
	rtree * polidx; // To be initialized after polygons and lines.
	
public:
	typedef std::tuple<int, string, rtree::const_query_iterator> close_t;
	Section(double cut, const pylist& points, 
		const pylist& polygons, const pylist& units, 
		const pylist& lines, const pylist& lnames ) {
		int npols = python::len(polygons);
		vector<box> envelopes;
		for ( int i = 0; i < npols; i++ ) {
			polygon pol;
			int nrings = python::len(polygons[i]);
			ring& outer = pol.outer();
			// Start filling the first ring.
			
			int nnodes = python::len(polygons[i][0]);
			for ( int k = 0; k < nnodes; k++ ) {
				pylist pypt = pylist(points[polygons[i][0][k]]);
				outer.push_back(point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])));
			}
			// Then fill the rest of the rings.
			if ( nrings > 1 ) { 
				pol.inners().resize(nrings-1);
				for ( int j = 1; j < nrings; j++ ) {
					ring& inner = pol.inners()[j-1];// jth ring.
					int nnodes = python::len(polygons[i][j]);
					for ( int k = 0; k < nnodes; k++ ) {
						pylist pypt = pylist(points[polygons[i][j][k]]);
						inner.push_back(point2(	python::extract<double>(pypt[0]), 
									python::extract<double>(pypt[1])));
					}
				}
			}
			if ( not geometry::is_valid(pol) or not geometry::is_simple(pol) ) {
				continue;
			}
			envelopes.push_back(pol);
			this->polygons.push_back(pol);
			this->units.push_back(python::extract<string>( units[i] ));
		}
		if ( envelopes.size() > 0 ) {
			polidx = new rtree( envelopes );
		} else {
			polidx = new rtree;
		}
		int nlines = python::len(lines);
		for ( int i = 0; i < nlines; i++ ) {
			line lin;
			int nnodes = python::len(lines[i]);
			for ( int j = 0; j < nnodes; j++ ) {
				pylist pypt = pylist(points[lines[i][j]]);
				lin.push_back(point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])));
			}
			if ( not geometry::is_valid(lin) or not geometry::is_simple(lin) ) {
				std::cerr << i << std::endl;
				continue;
			}
			this->lines.push_back(lin);
			this->lnames.push_back(python::extract<string>( units[i] ));
		}
	}
	pydict info() const {
		pydict res;
		res["polygons"] = this->polygons.size();
		res["lines"] = this->lines.size();
		return res;
	}
	close_t closest( const point2& p ) const {
		
	}
};

BOOST_PYTHON_MODULE(cpp)
{
	python::class_<Section>("Section", python::init<double, const pylist&, 
							const pylist&, const pylist&, 
							const pylist&, const pylist&>())
							.def("info", &Section::info);
}

