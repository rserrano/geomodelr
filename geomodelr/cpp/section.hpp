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

class Model;
class Match;

/* C++ section that queries points to polygons so much faster. */
class Section {
	friend Match;
	friend Model;
protected:
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
		if ( this->polidx != nullptr ) {
			for ( auto it = this->polidx->qbegin( geometry::index::intersects(bx) and predicates );
				it != this->polidx->qend(); it++ ) {
				// Check the actual distance to a polygon.
				int idx = it->second;
				double poldist = geometry::distance(this->polygons[idx], pt);
				if ( poldist <= distance ) {
					ret.push_back(std::make_pair(idx, poldist));
				}
			}
		}
		return ret;
	}
	template<typename Predicates>
	std::pair<int, double> closest( const point2& p, const Predicates& predicates ) const {
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
	map<wstring, vector<triangle_pt>> last_lines(bool is_back, double end);
public:
	Section(const wstring& name, double cut );
	virtual ~Section();
};

class SectionPython : Section {
	SectionPython(const wstring& name, double cut, const pylist& points, 
		      const pylist& polygons, const pylist& units, 
		      const pylist& lines, const pylist& lnames );
	
	pydict info() const;
	pytuple closest( const pyobject& pypt ) const;

};

#endif

