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
#include "section.hpp"


Section::Section(double cut, const pylist& points, 
	const pylist& polygons, const pylist& units, 
	const pylist& lines, const pylist& lnames ): valid(this) {
	size_t npols = python::len(polygons);
	vector<value> envelopes;
	for ( size_t i = 0; i < npols; i++ ) {
		polygon pol;
		size_t nrings = python::len(polygons[i]);
		ring& outer = pol.outer();
		// Start filling the first ring.
		
		size_t nnodes = python::len(polygons[i][0]);
		for ( size_t k = 0; k < nnodes; k++ ) {
			pylist pypt = pylist(points[polygons[i][0][k]]);
			outer.push_back(point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])));
		}
		// Then fill the rest of the rings.
		if ( nrings > 1 ) { 
			pol.inners().resize(nrings-1);
			for ( size_t j = 1; j < nrings; j++ ) {
				ring& inner = pol.inners()[j-1];// jth ring.
				size_t nnodes = python::len(polygons[i][j]);
				for ( size_t k = 0; k < nnodes; k++ ) {
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

	size_t nlines = python::len(lines);
	for ( size_t i = 0; i < nlines; i++ ) {
		line lin;
		size_t nnodes = python::len(lines[i]);
		for ( size_t j = 0; j < nnodes; j++ ) {
			pylist pypt = pylist(points[lines[i][j]]);
			lin.push_back(point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])));
		}
		if ( not geometry::is_valid(lin) or not geometry::is_simple(lin) ) {
			continue;
		}
		this->lines.push_back(lin);
		this->lnames.push_back(python::extract<string>( lnames[i] ));
	}
}

pydict Section::info() const {
	pydict res;
	res["polygons"] = this->polygons.size();
	res["lines"] = this->lines.size();
	return res;
}

pytuple Section::closest( const pyobject& pypt ) const {
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

bool ValidUnit::operator()(const value& b) const {
	const string& unit = this->section->units[b.second];
	return unit != "NONE" and unit != ""; 
}
