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

struct Section;

struct ValidUnit {
	const Section * section;
	ValidUnit( const Section * section ):section(section) 
	{
		
	}
	bool operator() (const value& b) const;
	
};

pylist match( const Section& a, const Section& b );

class Section {
	double cut;
	vector<polygon> polygons;
	vector<string> units;
	vector<line> lines;
	vector<string> lnames;
	rtree * polidx; // To be initialized after polygons and lines.
	friend class ValidUnit;
	friend bool match;
	ValidUnit valid;
public:
	Section(double cut, const pylist& points, 
		const pylist& polygons, const pylist& units, 
		const pylist& lines, const pylist& lnames ): valid(this) {
		int npols = python::len(polygons);
		vector<value> envelopes;
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
			// Calculate the envelope and add it to build the rtree layer.
			box env;
			geometry::envelope(pol, env);
			envelopes.push_back(std::make_pair(env, this->polygons.size()));
			// Now add the actual polygon and its unit.
			this->polygons.push_back(pol);
			this->units.push_back(python::extract<string>( units[i] ));
		}
		// Build the rtree.
		if ( envelopes.size() > 0 ) {
			this->polidx = new rtree( envelopes );
		} else {
			this->polidx = nullptr;
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
			this->lnames.push_back(python::extract<string>( lnames[i] ));
		}
	}
	pydict info() const {
		pydict res;
		res["polygons"] = this->polygons.size();
		res["lines"] = this->lines.size();
		return res;
	}
	pytuple closest( const pyobject& pypt ) const {
		double x = python::extract<double>(pypt[0]);
		double y = python::extract<double>(pypt[1]);
		point2 p(x, y);
		
		if ( this->polidx == nullptr )
			return python::make_tuple(-1, "NONE");
		
		double maxboxdist = 0.0;
		double mindist = std::numeric_limits<double>::infinity();
		
		int minidx = -1;
		int knear = 1; 
		
		bool new_to_check;
		do {
			int n = 0; 
			new_to_check = false;
			for ( 	auto it = this->polidx->qbegin( geometry::index::nearest(p, knear) and geometry::index::satisfies(this->valid) );
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
		if ( minidx == -1 ) {
			return python::make_tuple(-1, "NONE");
		}
		return python::make_tuple(minidx, this->units[minidx]);
	}
};

bool ValidUnit::operator()(const value& b) const{
	const string& unit = this->section->units[b.second];
	return unit != "NONE" and unit != ""; 
}

std::tuple<vector<std::pair<int, int>>, vector<int>, vector<int>>  match( const Section& a, const Section& b ) {
	map<string, vector<int>> units_a;
	map<string, vector<int>> units_b;
	vector<bool> sel_a(a.polygons.size(), false);
	vector<bool> sel_b(b.polygons.size(), false);
	for ( int i = 0; i < a.polygons.size(); i++ ) {
		units_a[a.units[i]].push_back(i);
	}
	for ( int i = 0; i < b.polygons.size(); i++ ) {
		units_b[b.units[i]].push_back(i);
	}
	vector<std::pair<int, int>> m;
	for ( auto it = units_a.begin(); it != units_a.end(); it++ ) {
		if ( units_b.find(it->first) != units_b.end() ) {
			vector<int>& pols_a = it->second;
			vector<int>& pols_b = units_b[it->first];
			for ( int i = 0; i < pols_a.size(); i++ )
			{
				for ( int j = 0; j < pols_b.size(); j++ ) {
					if ( geometry::intersects(a.polygons[pols_a[i]], b.polygons[pols_b[j]]) ) {
						sel_a[pols_a[i]] = true;
						sel_b[pols_b[j]] = true;
						m.push_back(make_pair(pols_a[i], pols_b[j]));
					}
				}
			}
		}
	}
	vector<int> sa;
	for ( int i = 0; i < sel_a.size(); i++ ) {
		if ( not sel_a[i] ) {
			sa.push_back(i);
		}
	}
	vector<int> sb;
	for ( int i = 0; i < sel_b.size(); i++ ) {
		if ( not sel_b[i] ) {
			sb.push_back(i);
		}
	}
	return std::make_tuple( m, sa, sb );
}

class Model {
	Model(const pyobject& geojson) {
	}
}

BOOST_PYTHON_MODULE(cpp)
{
	python::class_<Section>("Section", python::init<double, const pylist&, 
							const pylist&, const pylist&, 
							const pylist&, const pylist&>())
							.def("info", &Section::info)
							.def("closest", &Section::closest);
}

