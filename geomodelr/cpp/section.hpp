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
#ifndef GEOMODELR_SECTION_HPP
#define GEOMODELR_SECTION_HPP

#include "basic.hpp"
#include "polygon.hpp"
#include <stdlib.h>
#include <math.h> 
#include <set>

class Model;
class Match;
class Polygon;
class ModelPython;
class MatchPython;


/* C++ section that queries points to polygons so much faster. */
class Section {
	friend Match;
	friend Model;
	friend ModelPython;
	friend MatchPython;
	friend Polygon;

protected:
	wstring name;
	double cut;
	bbox2 bbox;
	vector<wstring> units;
	vector<line> lines;
	vector<std::pair<point2, point2>> line_ends;
	std::set<line_anchor> anchored_lines;
	vector<wstring> lnames;
	rtree_f * polidx; // To be initialized after polygons and lines.
	rtree_l * fault_lines;
	vector<Polygon *> poly_trees;
	const map<wstring, wstring> * params;
	
	template<typename Predicates>
	vector<std::pair<int, double>> closer_than( const point2& pt, double distance, const Predicates& predicates ) const {
		point2 mx(gx(pt) + distance, gy(pt) + distance);
		point2 mn(gx(pt) - distance, gy(pt) - distance);
		box bx(mn, mx);
		vector<std::pair<int, double>> ret;
		if ( this->polidx != nullptr ) {
			for ( auto it = this->polidx->qbegin( geometry::index::intersects(bx) and geometry::index::satisfies(predicates) );
				it != this->polidx->qend(); it++ ) {
				// Check the actual distance to a polygon.
				int idx = g2(*it);
				// double poldist = geometry::distance(poly_trees[idx]->boost_poly, pt);
				double poldist = poly_trees[idx]->distance_point(pt);
				if ( poldist <= distance ) {
					ret.push_back(std::make_pair(idx, poldist));
				}
			}
		}
		return ret;
	}
public:

	template<typename Predicates>
	std::pair<int, double> closest( const point2& p, const Predicates& predicates ) const {
		if ( this->polidx == nullptr )
			return std::make_pair(-1, std::numeric_limits<double>::infinity());
		
		double maxboxdist = 0.0;
		double mindist = std::numeric_limits<double>::infinity();
		
		int minidx = -1;
		int knear = 1;
		
		bool new_to_check;
		std::set<int> checked;

		
		do {
			new_to_check = false;
			for ( auto it = this->polidx->qbegin( geometry::index::nearest(p, knear) and geometry::index::satisfies(predicates) );
				it != this->polidx->qend(); it++ ) {
				// Skip already checked.
				int idx = g2(*it);
				if ( checked.find( idx ) != checked.end() )
				{
					continue;
				}
				checked.insert(idx);
				// Check if new polygons where checked.
				new_to_check = true;
				
				// Check the maximum distance from the box to the point.
				// That distance is always lower than the distance to the polygon.
				double boxdist = geometry::distance(p, g0(*it));
				maxboxdist = std::max(boxdist, maxboxdist);
				
				// Then check the minimum actual distance to a polygon.
				// double poldist = geometry::distance(poly_trees[idx]->boost_poly, pt);
				double poldist = poly_trees[idx]->distance_point(p);
				if ( poldist < mindist ) {
					mindist = poldist;
					minidx = idx;
				}
			}
			// Increase the number of knear.
			knear *= 2;
			// Do it until none was checked or we have checked boxes beyond the closest polygon.
		} while ( new_to_check && maxboxdist < mindist );
		
		//if (minidx==-1){std::cerr << "No hay poligono cerca" << std::endl;}

		return std::make_pair(minidx, mindist); 
	}
	std::pair<int, double> closest( const point2& ) const;
	void set_params( const map<wstring, wstring> * params );
	std::tuple<map<wstring, vector<triangle_pt>>, map<wstring, vector<size_t>>> last_lines(bool is_back, double end);
	Section( const wstring& name, double cut, const bbox2& bbox );
	virtual ~Section();
};

class SectionPython : public Section {
	map<wstring, wstring> local_params; // For single section use only.
public:
	SectionPython( const wstring& name, double cut, const pyobject& bbox, const pylist& points, 
		       const pylist& polygons, const pylist& units, const pylist& lines,
		       const pylist& lnames, const pylist& anchored_lines );
	pydict info() const;
	pytuple closest( const pyobject& pypt ) const;
	void set_params( const pydict& params );
	pydict get_params( ) const;
	double distance_poly(const pylist& pypt, int poly_idx) const;
};

class GeologicalMapPython: public SectionPython {
public:
	GeologicalMapPython( const pyobject& bbox, const pylist& points, 
			     const pylist& polygons, const pylist& units, const pylist& lines,
			     const pylist& lnames, const pylist& anchored_lines );
};

void extend_line( bool beg, const bbox2& bbox, line& l );
pylist test_extend_line( bool beg, const pyobject& bbox, const pylist& line );

#endif

