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
#include<cmath>

Section::~Section()
{
	if ( this->polidx != nullptr ) {
		delete this->polidx;
	}
}

Section::Section( const wstring& name, double cut, const bbox2& bbox ): name(name), cut(cut), bbox( bbox ), polidx(nullptr) 
{
	
}

std::pair<int, double> Section::closest( const point2& pt ) const {
	return this->closest(pt, always_true);
}
	

SectionPython::SectionPython(const wstring& name, double cut, 
	const pyobject& bbox, const pylist& points, 
	const pylist& polygons, const pylist& units, 
	const pylist& lines, const pylist& lnames,
	const pylist& anchored_lines ): Section( name, cut, std::make_tuple( std::tuple<double, double>(python::extract<double>(bbox[0]), python::extract<double>(bbox[1])),
									     std::tuple<double, double>(python::extract<double>(bbox[2]), python::extract<double>(bbox[3])) ) )
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
	// Add the lines.
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
	
	// Add which are the anchors, (for signaling later), but also extend the previously added lines.
	size_t nanchs = python::len(anchored_lines);
	for ( size_t i = 0; i < nanchs; i++ ) {
		int lidx = python::extract<int>(anchored_lines[i][0]);
		bool lbeg = python::extract<bool>(anchored_lines[i][1]);
		
		line_anchor la = std::make_pair( lidx, lbeg );
		
		
		if ( lidx < 0 or size_t(lidx) >= this->lines.size() ) {
			if ( geomodelr_verbose ) {
				std::wcerr << L"Wrong input to anchored_lines\n";
			}
			continue;
		}
		
		line& ln = this->lines[lidx];
		try {
			extend_line( lbeg, this->bbox, ln );
			this->anchored_lines.insert(la);
		} catch ( const GeomodelrException& e ) {
			if ( geomodelr_verbose ) {
				std::cerr << e.what() << std::endl;
			}
		}
		
		
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



std::tuple<map<wstring, vector<triangle_pt>>, 
	   map<wstring, vector<size_t>>> Section::last_lines( bool is_back, double end ) {
	map<wstring, vector<triangle_pt>> ret;
	map<wstring, vector<size_t>> extended;
	for ( size_t i = 0; i < this->lines.size(); i++ )
	{
		const wstring& name = this->lnames[i];
		if ( name == L"" )
		{
			continue;
		}
		// Avoid extending more than one fault plane per name.
		if ( ret.find( name ) != ret.end() ) 
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
		if ( this->anchored_lines.find( std::make_pair(i, true) ) != this->anchored_lines.end() ) {
			extended[name].push_back(0);
			extended[name].push_back(1);
		}
		if ( this->anchored_lines.find( std::make_pair(i, false) ) != this->anchored_lines.end() ) {
			size_t n = ret[name].size();
			extended[name].push_back(n-2);
			extended[name].push_back(n-1);
		}
	}
	return std::make_tuple(ret, extended);
}

void extend_line( bool beg, const bbox2& bbox, line& l ) {

	if ( l.size() <= 1 ) {
		throw GeomodelrException("An input line has a single point.");
	}
	point2 vct;
	point2 pt;
	if ( beg ) {
		vct = l[0];
		pt = l[1];
		// Obtain the vector.
		geometry::subtract_point(vct, pt);
		pt = l[0];
	} else {
		vct = l[l.size()-1];
		pt = l[l.size()-2];
		// Obtain the vector.
		geometry::subtract_point(vct, pt);
		pt = l[l.size()-1];
	}
	double x;
	double minx = std::numeric_limits<double>::infinity();
	if ( std::fabs( gx( vct ) ) > tolerance ) {
		x = ( g0( g0(bbox) ) - gx(pt) )/gx(vct);
		if ( x > 0 ) {
			minx = std::min(minx, x);
		}
		x = ( g0( g1(bbox) ) - gx(pt) )/gx(vct);
		if ( x > 0 ) {
			minx = std::min(minx, x);
		}
	}
	
	if ( std::fabs( gy( vct ) ) > tolerance ) {
		x = ( g1( g0(bbox) ) - gy(pt) )/gy(vct);
		if ( x > 0 ) {
			minx = std::min(minx, x);
		}
		x = ( g1( g1(bbox) ) - gy(pt) )/gy(vct);
		if ( x > 0 ) {
			minx = std::min(minx, x);
		}
	}
	
	if ( not std::isfinite( minx ) ) {
		if ( not ( std::fabs( gy( vct ) ) > tolerance or std::fabs( gx( vct ) ) > tolerance ) ) {
			throw GeomodelrException("Could not determine direction of line.");
		}
		throw GeomodelrException("Could not extend line.");
	}
	
	if ( std::fabs( minx ) < tolerance ) {
		throw GeomodelrException("The line actually goes to the bounds and does not need modification.");
	}
	
	geometry::multiply_value(vct, minx);
	geometry::add_point(pt, vct);
	
	// Check that the point falls inside (or very close to) the bounding box.
	if ( not ( gx( pt )-g0(g0(bbox)) >= -tolerance and gx(pt)-g0(g1(bbox)) <= tolerance and
	     gy(pt)-g1(g0(bbox)) >= -tolerance and gy(pt)-g1(g1(bbox)) <= tolerance ) ) {
		throw GeomodelrException("The point is not inside the bounding box.");
	}
	if ( beg ) {
		l.insert( l.begin(), pt );
	} else {
		l.push_back( pt );
	}
}

pylist test_extend_line( bool beg, const pyobject& bbox, const pylist& pl ) {
	auto b = std::make_tuple( std::tuple<double, double>( python::extract<double>(bbox[0]), python::extract<double>(bbox[1]) ),
				  std::tuple<double, double>( python::extract<double>(bbox[2]), python::extract<double>(bbox[3]) ) );
	
	line l;
	for ( int i = 0; i < python::len( pl ); i++ ) {
		l.push_back( point2( python::extract<double>(pl[i][0]), python::extract<double>(pl[i][1]) ) );
	}
	extend_line( beg, b, l );
	pylist ret;
	for ( size_t i = 0; i < l.size(); i++ ) {
		ret.append( python::make_tuple( gx( l[i] ), gy( l[i] ) ) );
	}
	return ret;
}
