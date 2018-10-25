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
#include <cmath>
#include <iostream>

Section::~Section()
{
	if ( this->polidx != nullptr ) {
		delete this->polidx;
	}

	if ( this->fault_lines != nullptr ) {
		delete this->fault_lines;
	}

	for ( Polygon * p: this->poly_trees ) {
		delete p;
	}
}

Section::Section( const wstring& name, double cut, const bbox2& bbox ): name(name), cut(cut), bbox( bbox ), polidx(nullptr), fault_lines(nullptr), params(nullptr)
{

}

std::pair<int, double> Section::closest( const point2& pt ) const {
	return this->closest(pt, always_true);
}

void Section::set_params( const map<wstring, wstring> * params ) {
	this->params = params;
	auto kv = this->params->find( L"faults" );
	// Make all polygons use the corresponding function.
	if ( kv != this->params->end() ) {
		wstring fault_method = kv->second;
		for ( Polygon * p: this->poly_trees ) {
			p->set_distance_function( fault_method );
		}
	}
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
	vector<value_l> f_segs;
	
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
			point2 aux = point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1]));
			if ( outer.size() ) {
				if ( geometry::distance(outer.back(), aux) < boost_tol ) {
					continue;
				}
			}
			outer.push_back(aux);
		}
		if ( outer.size() < 3 ) {
			if ( geomodelr_verbose ) {
				std::wcerr << L"non valid polygon in section " << name << L" from with unit " << unit << " has too few points\n";
				continue;
			}
		}
		// Then fill the rest of the rings.
		if ( nrings > 1 ) { 
			vector<ring>& inners = pol.inners();
			for ( size_t j = 1; j < nrings; j++ ) {
				ring inner;
				size_t nnodes = python::len(polygons[i][j]);
				for ( size_t k = 0; k < nnodes; k++ ) {
					pylist pypt = pylist(points[polygons[i][j][k]]);
					point2 aux = point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1]));
					if ( outer.size() ) {
						if ( geometry::distance(outer.back(), aux) < boost_tol ) {
							continue;
						}
					}
					inner.push_back(aux);
				}
				if ( inner.size() > 2 ) {
					inners.push_back(inner);
				}
			}
		}
		
		geometry::validity_failure_type failure;
		if ( not geometry::is_valid(pol, failure) ) {
			geometry::correct(pol);
			string reason;
			if ( not geometry::is_valid(pol, reason) ) {
				if ( geomodelr_verbose ) {
					std::wstring wide(reason.size(), L' ');
					std::copy(reason.begin(), reason.end(), wide.begin());
					std::wcerr << L"non valid polygon in section " << name << L" from with unit " << unit << "\nreason is: " << wide << "\n";
				}
				// continue; not avoiding non valid polygons, as they have been validated by shapely. Somehow these polygons get wronged.
			}
		}
		if ( not geometry::is_simple(pol) ) {
			if ( geomodelr_verbose ) {
				std::wcerr << L"ignored non simple polygon in section " << name << L" from with unit " << unit << "\n";
			}
			continue;
		}
		
		// Calculate the envelope and add it to build the rtree layer.
		box env;
		geometry::envelope(pol, env);
		envelopes.push_back(std::make_tuple(env, unit, this->poly_trees.size()));
		
		// Now add the actual polygon and its unit.
		this->poly_trees.push_back(new Polygon(pol, env, this));
		this->units.push_back(unit);
	}

	// Build the rtree.
	if ( envelopes.size() > 0 ) {
		this->polidx = new rtree_f( envelopes );
	}
	// Add the lines.
	size_t nlines = python::len(lines);
	map<size_t, bool> ancl;
	// Add which are the anchors, (for signaling later), but also extend the previously added lines.
	size_t nanchs = python::len(anchored_lines);
	for ( size_t i = 0; i < nanchs; i++ ) {
		int lidx = python::extract<int>(anchored_lines[i][0]);
		bool lbeg = python::extract<bool>(anchored_lines[i][1]);
		line_anchor la = std::make_pair( size_t(lidx), lbeg );
		if ( lidx < 0 or size_t(lidx) >= nlines ) {
			if ( geomodelr_verbose ) {
				std::wcerr << L"Wrong input to anchored_lines\n";
			}
			continue;
		}
		ancl.insert( la );
	}
	
	int count_lines = 0;
	for ( size_t i = 0; i < nlines; i++ ) {
		line lin;
		size_t nnodes = python::len(lines[i]);
		for ( size_t j = 0; j < nnodes; j++ ) {
			
			pylist pypt = pylist(points[lines[i][j]]);
			point2 aux = point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1]));
			if ( lin.size() ) {
				if ( geometry::distance(lin.back(), aux) < boost_tol ) {
					continue;
				}
			}
			lin.push_back(aux);
		}
		if ( not geometry::is_valid(lin) or not geometry::is_simple(lin) ) {
			if ( geomodelr_verbose ) {
				wstring name = python::extract<wstring>( lnames[i] );
				std::wcerr << L"ignoring line " << name << L" is valid: " 
					   << geometry::is_valid(lin) << L" is simple "<< geometry::is_simple(lin) << L"\n" ;
			}
			continue; // This should guarantee that the line has at least 2 points.
		}
		
		for ( auto it = ancl.find( i ); it != ancl.end(); it++ ) {
			if ( it->first != i ) {
				break;
			}
			try {
				// First extend line.
				extend_line( it->second, this->bbox, lin );
				// Then insert to the given anchors.
				this->anchored_lines.insert( *it );
			} catch ( const GeomodelrException& e ) {
				if ( geomodelr_verbose ) {
					std::cerr << e.what() << std::endl;
				}
			}
		}
		
		this->lines.push_back(lin);
		// Create the ends by projecting these a little bit.
		point2& end = lin[lin.size()-1], beg = lin[0];
		point2& pend = lin[lin.size()-2], nbeg = lin[1];
		point2 ve = end, vb = beg;
		geometry::subtract_point( ve, pend );
		geometry::subtract_point( vb, nbeg );
		geometry::divide_value( ve, std::sqrt( gx(ve)*gx(ve) + gy(ve)*gy(ve) ) );
		geometry::multiply_value( ve, 2.0*boost_tol );
		geometry::divide_value( vb, std::sqrt( gx(vb)*gx(vb) + gy(vb)*gy(vb) ) );
		geometry::multiply_value( vb, 2.0*boost_tol );
		
		geometry::add_point( ve, end );
		geometry::add_point( vb, beg );
		
		this->line_ends.push_back(std::make_pair(vb, ve));
		this->lnames.push_back(python::extract<wstring>( lnames[i] ));
		
		vector<value_l> fault_segments;
		
		for ( size_t i = 0; i < lin.size()-1; i++ ) {	
			fault_segments.push_back(std::make_tuple(line_segment(lin[i],lin[i+1]),count_lines));
		}
		
		f_segs.reserve(f_segs.size() + distance(fault_segments.begin(),fault_segments.end()));
		f_segs.insert(f_segs.end(), fault_segments.begin(),fault_segments.end());
		
		count_lines++;
	}
	this->fault_lines =  new rtree_l( f_segs );
}

pydict SectionPython::info() 
const {
	pydict res;
	res["polygons"] = this->poly_trees.size();
	res["lines"] = this->lines.size();
	return res;
}

pytuple SectionPython::closest( const pyobject& pypt ) 
const {
	double x = python::extract<double>(pypt[0]);
	double y = python::extract<double>(pypt[1]);
	point2 p(x, y);
	
	std::pair<int, double> cls = Section::closest(p);
	
	if ( cls.first == -1 ) {
		return python::make_tuple(wstring(L"NONE"), cls.second);
	}
	return python::make_tuple(this->units[cls.first], cls.second);
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
			throw GeomodelrException("fault not extended: could not determine direction of line.");
		}
		throw GeomodelrException("fault not extended: could not extend line.");
	}
	
	if ( std::fabs( minx ) < tolerance ) {
		throw GeomodelrException("fault not extended: the line actually goes to the bounds and does not need modification.");
	}
	
	geometry::multiply_value(vct, minx);
	geometry::add_point(pt, vct);
	
	// Check that the point falls inside (or very close to) the bounding box.
	if ( not ( gx( pt )-g0(g0(bbox)) >= -tolerance and gx(pt)-g0(g1(bbox)) <= tolerance and
	     gy(pt)-g1(g0(bbox)) >= -tolerance and gy(pt)-g1(g1(bbox)) <= tolerance ) ) {
		throw GeomodelrException("fault not extended: the point is not inside the bounding box.");
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

double SectionPython::distance_poly(const pylist& pypt, int idx) const{

	double x = python::extract<double>(pypt[0]);
	double y = python::extract<double>(pypt[1]);
	return this->poly_trees[idx]->distance_point(point2(x,y));
}

void SectionPython::set_params(const pydict& params) {
	pylist keys = params.keys();
	for ( int i = 0; i < python::len( keys ); i++ ) {
		this->local_params[python::extract<wstring>(keys[i])] = python::extract<wstring>(params[keys[i]]);
	}
	Section::set_params(&(this->local_params));
}

pydict SectionPython::get_params() const {
	pydict out;
	for ( auto& kv: this->local_params ) {
		out[kv.first] = kv.second;
	}
	return out;
}

GeologicalMapPython::GeologicalMapPython( const pyobject& bbox, const pylist& points, 
					  const pylist& polygons, const pylist& units, const pylist& lines,
					  const pylist& lnames, const pylist& anchored_lines ):SectionPython( wstring(L"Geological Map"), 0.0, 
					  								      bbox, points, polygons, 
													      units, lines, lnames, 
													      anchored_lines )
{
	
}
