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

Section::~Section()
{
	if ( this->polidx != nullptr ) {
		delete this->polidx;
	}
}

Section::Section( const wstring& name, double cut ): name(name), cut(cut), polidx(nullptr) 
{
}

std::pair<int, double> Section::closest( const point2& pt ) const {
	return this->closest(pt, always_true);
}
	

SectionPython::SectionPython(const wstring& name, double cut, const pylist& points, 
	const pylist& polygons, const pylist& units, 
	const pylist& lines, const pylist& lnames ): Section( name, cut )
{
	size_t npols = python::len(polygons);
	vector<value_f> envelopes;
	for ( size_t i = 0; i < npols; i++ ) {
		polygon pol;
		wstring unit = python::extract<wstring>( units[i] );
		if ( unit == L"NONE" or unit == L"" ) {
			continue;
		}
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
		
		geometry::validity_failure_type failure;
		if ( not geometry::is_valid(pol, failure) ) {
			geometry::correct(pol);
			string reason;
			if ( not geometry::is_valid(pol, reason) ) {
				if ( geomodelr_verbose ) {
					std::wcerr << L"non valid polygon in section " << name << L" from with unit " << unit << L" valid: " << (geometry::is_valid(pol)?L"true":L"false")
						   << L" simple: " << (geometry::is_simple(pol)?L"true":L"false") << "\n";
				}
				// continue; not avoiding non valid polygons, as they have been validated by shapely. Somehow these polygons get wronget.
			}
		}
		if ( not geometry::is_simple(pol) ) {
			continue;
		}
		// Calculate the envelope and add it to build the rtree layer.
		box env;
		geometry::envelope(pol, env);
		envelopes.push_back(std::make_tuple(env, unit, this->polygons.size()));
		
		// Now add the actual polygon and its unit.
		this->polygons.push_back(pol);
		this->units.push_back(unit);
	}
	// Build the rtree.
	if ( envelopes.size() > 0 ) {
		this->polidx = new rtree_f( envelopes );
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
		this->lnames.push_back(python::extract<wstring>( lnames[i] ));
	}
}

pydict SectionPython::info() 
const {
	pydict res;
	res["polygons"] = this->polygons.size();
	res["lines"] = this->lines.size();
	return res;
}

pytuple SectionPython::closest( const pyobject& pypt ) 
const {
	double x = python::extract<double>(pypt[0]);
	double y = python::extract<double>(pypt[1]);
	point2 p(x, y);
	
	std::pair<int, int> cls = Section::closest(p);
	
	if ( cls.first == -1 ) {
		return python::make_tuple(-1, wstring(L"NONE"));
	}
	
	return python::make_tuple(cls.first, this->units[cls.first]);

}



map<wstring, vector<triangle_pt>> Section::last_lines( bool is_back, double end ) {
	map<wstring, vector<triangle_pt>> ret;
	for ( size_t i = 0; i < this->lines.size(); i++ )
	{
		const wstring& name = this->lnames[i];
		if ( name == L"" ) 
		{
			continue;
		}
		for ( size_t j = 1; j < this->lines[i].size(); j++ ) {
			point3 p0(gx(this->lines[i][j-1]), gy(this->lines[i][j-1]), end);
			point3 p1(gx(this->lines[i][j]),   gy(this->lines[i][j]),   end);
			point3 p2(gx(this->lines[i][j-1]), gy(this->lines[i][j-1]), this->cut);
			point3 p3(gx(this->lines[i][j]),   gy(this->lines[i][j]),   this->cut);
			
			if ( ! is_back ) {
				std::swap(p0, p2);
				std::swap(p1, p3);
			}
			
			triangle_pt t1(p0, p1, p2);
			triangle_pt t2(p1, p3, p2);
			
			ret[name].push_back(t1);
			ret[name].push_back(t2);
		}
	}
	return ret;
}

