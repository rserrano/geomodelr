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

#include "match.hpp"
#include <cmath>


Match::Match( const Section * a, const Section * b )
:a(a), b(b), faultidx(nullptr)
{
}

Match::~Match( )
{
	if ( this->faultidx != nullptr ) 
	{
		delete this->faultidx;
	}
}

void Match::set( const vector<std::pair<int, int>>& match ){
	this->a_to_b.clear();
	this->b_to_a.clear();
	
	for ( size_t i = 0; i < match.size(); i++ ) {
		a_to_b[match[i].first].push_back(match[i].second);
		b_to_a[match[i].second].push_back(match[i].first);
	}
}



void Match::match_polygons() {
	map<wstring, vector<int>> units_a;
	map<wstring, vector<int>> units_b;
	
	for ( size_t i = 0; i < this->a->polygons.size(); i++ ) {
		units_a[this->a->units[i]].push_back(i);
	}
	
	for ( size_t i = 0; i < this->b->polygons.size(); i++ ) {
		units_b[this->b->units[i]].push_back(i);
	}
	
	vector<std::pair<int, int>> m;
	for ( auto it = units_a.begin(); it != units_a.end(); it++ ) {
		if ( units_b.find(it->first) != units_b.end() ) {
			vector<int>& pols_a = it->second;
			vector<int>& pols_b = units_b[it->first];
			for ( size_t i = 0; i < pols_a.size(); i++ )
			{
				for ( size_t j = 0; j < pols_b.size(); j++ ) {
					if ( geometry::intersects(this->a->polygons[pols_a[i]], this->b->polygons[pols_b[j]]) ) {
						m.push_back(std::make_pair(pols_a[i], pols_b[j]));
					}
				}
			}
		}
	}

	this->set( m );
}



std::tuple<int, int, int> Match::crosses_triangles(const point2& point, double cut) const {
	if ( this->faultidx == nullptr ) {
		return std::make_tuple(0, -1, -1);
	}
	for ( auto it = this->faultidx->qbegin( geometry::index::contains(point) ); it != this->faultidx->qend(); it++ ) {
		const auto fl = this->faults.find(g1(*it));
		const AlignedTriangle& tr = fl->second[g2(*it)];
		int side = tr.crosses_triangle(point, cut);
		if ( side != 0 ) 
		{
			const auto lni = this->rel_faults.find(g1(*it));
			return std::make_tuple(side, g0(lni->second), g1(lni->second));
		}
	}
	return std::make_tuple(0, -1, -1);
}

AlignedTriangle::AlignedTriangle(const std::tuple<point3, point3, point3>& triangle) {
	const point3& p0 = g0(triangle);
	const point3& p1 = g1(triangle);
	const point3& p2 = g2(triangle);
	point3 v1 = p1;
	point3 v2 = p2;
	// Find vectors and its cross product.
	geometry::subtract_point(v1, p0);
	geometry::subtract_point(v2, p0);
	this->normal = point3(gy(v1)*gz(v2)-gz(v1)*gy(v2), gz(v1)*gx(v2)-gx(v1)*gz(v2), gx(v1)*gy(v2)-gy(v1)*gx(v2));
	// Find norm of the vector.
	double norm = std::sqrt(gx(this->normal)*gx(this->normal) + gy(this->normal)*gy(this->normal) + gz(this->normal)*gz(this->normal));
	// Reverse if below 0.0
	if ( gz(normal) < 0 ) {
		norm *= -1.0;
	}
	// divide by norm.
	geometry::divide_value(this->normal, norm);
	
	this->point  = p0;
	ring& outer = this->triangle.outer();
	outer.push_back(point2(gx(p0), gy(p0)));
	outer.push_back(point2(gx(p1), gy(p1)));
	outer.push_back(point2(gx(p2), gy(p2)));
}

int AlignedTriangle::crosses_triangle(const point2& point, double cut ) const {
	if ( geometry::within(point, this->triangle) ) {
		point3 pt(gx(point), gy(point), cut);
		geometry::subtract_point(pt, this->point);
		double d = geometry::dot_product(pt, this->normal);
		if ( d < 0 ) {
			return -1;
		} else {
			return 1;
		}
	}
	return 0;
}

point3 angles( const point3& t ) 
{
	/*
	Given a triangle it finds the angles of the triangle.
	*/
	point3 st(gx(t)*gx(t), gy(t)*gy(t), gz(t)*gz(t));
	point3 cs(std::acos((gy(st) + gz(st) - gx(st))/(2.0*gy(t)*gz(t))),
		  std::acos((gx(st) + gy(st) - gz(st))/(2.0*gx(t)*gy(t))),
		  0.0);
	geometry::set<2>(cs, M_PI - (gx(cs) + gy(cs)));
	return cs;
}

std::tuple<edge, edge> next_and_check(const triangle& tri, const edge& edg, size_t na)
{
	/*
	Gets the next edge and also the edge that now has been finished.
	*/
	if ( g1(tri) < int(na) ) {
		// Case g0(tri), g1(tri) are in the same edge g2(tri) the contrary.
		// (g0(tri), g1(tri)) need to be removed.
		if ( g0(edg) == g0(tri) ) {
			// Case g0(edg) is g0(tri), then next should be (g1(tri), g2(tri))
			return std::make_tuple( std::make_tuple(g1(tri), g2(tri)), std::make_tuple(g0(tri), g1(tri)) );
		} else {
			// Case g0(edg) is g1(tri), then next should be (g0(tri), g2(tri))
			return std::make_tuple( std::make_tuple(g0(tri), g2(tri)), std::make_tuple(g0(tri), g1(tri)) );
		}
	} else {
		// Case g1(tri), g2(tri) are in the same edge, g0(tri) the contrary.
		// (g1(tri), g2(tri)) need to be removed.
		if ( g1(edg) == g1(tri) ) {
		    // Case g0(edg) is g0(tri), then next should be (g1(tri), g2(tri))
		    return std::make_tuple( std::make_tuple(g0(tri), g2(tri)), std::make_tuple(g1(tri), g2(tri)) );
		} else {
		    // Case g0(edg) is g1(tri), then next should be (g0(tri), g2(tri))
		    return std::make_tuple( std::make_tuple(g0(tri), g1(tri)), std::make_tuple(g1(tri), g2(tri)) );
		}
	}
}

vector<triangle> test_start( const vector<point3>& pa, const vector<point3>& pb, bool binv )
{
	auto next_a = [&] ( int i ) -> int {
		return size_t(i+1) < pa.size() ? i+1:-1;
	};
	
	auto next_b = [&] ( int i ) -> int {
		if ( binv ) {
			return i-1;
		} else {
			return size_t(i+1) < pb.size() ? i+1:-1;
		}
	};
	
	auto b_idx = [&] ( int i ) -> int {
		return pa.size() + i;
	};
	
	double length_a = 0.0;
	double length_b = 0.0;
	
	for ( size_t i = 0; i < pa.size()-1; i++ ) {
		length_a += geometry::distance(pa[i], pa[i+1]);
	}
	
	for ( size_t i = 0; i < pb.size()-1; i++ ) {
		length_b += geometry::distance(pb[i], pb[i+1]);
	}
	
	vector<triangle> result;
	
	int n_a, n_b;
	int c_a = 0;
	int c_b = binv ? pb.size()-1 : 0;
	n_a = next_a(c_a);
	n_b = next_b(c_b);
	
	double prop_a = 0.0;
	double prop_b = 0.0;
	
	while ( n_a != -1 || n_b != -1 ) { // Get next points if advance a, and if advance b, and check if they are not -1 both.
		double next_prop_a = std::numeric_limits<double>::infinity();
		double next_prop_b = std::numeric_limits<double>::infinity();
		
		if ( n_a != -1 ) {
			next_prop_a = prop_a + (geometry::distance( pa[n_a], pa[c_a] )/length_a);
		}
		
		if ( n_b != -1 ) {
			next_prop_b = prop_b + (geometry::distance( pb[n_b], pb[c_b] )/length_b);
		}
		
		// Check angle if advance b.
		if ( next_prop_b < next_prop_a ) {
			result.push_back( triangle( c_a, b_idx(c_b), b_idx(n_b) ) );
			c_b = n_b;
			prop_b = next_prop_b;
		} else {
			result.push_back( triangle( c_a, n_a, b_idx(c_b) ) );
			c_a = n_a;
			prop_a = next_prop_a;
		}

		n_a = next_a(c_a);
		n_b = next_b(c_b);
	}
	return result;
}

vector<triangle> faultplane_for_lines(const vector<point3>& l_a, const vector<point3>& l_b)
{
	// Get the faults plane between lines la, lb.
	
	point3 va = l_a.back();
	geometry::subtract_point(va, l_a[0]);
	geometry::divide_value( va, std::sqrt( geometry::dot_product( va, va ) ) );
	
	point3 vb = l_b.back();
	geometry::subtract_point(vb, l_b[0]);
	geometry::divide_value( vb, std::sqrt( geometry::dot_product( vb, vb ) ) );
	
	double angle = std::acos(geometry::dot_product( va, vb ));
	//std::cerr << "angle " << angle << "\n";
	if ( angle <= M_PI/2.0 ) {
		return test_start( l_a, l_b, false );
	}
	return test_start( l_a, l_b, true  );
}

std::tuple<map<wstring, vector<triangle_pt>>, map<wstring, vector<size_t>>> Match::match_lines( const map<wstring, wstring>& feature_types )
{
	/*
	Creates the fault planes given the cross sections with faults with the same name.
	*/
	// Create map of related faults.
	map<wstring, std::tuple<int, int>>& rel_faults = this->rel_faults;
	for ( size_t i = 0; i < this->a->lines.size(); i++ )
	{
		const wstring& name = this->a->lnames[i];
		if ( name != L"" ) {
			rel_faults[name] = std::make_tuple(i, -1);
		}
	}
	
	for ( size_t i = 0; i < this->b->lines.size(); i++ ) {
		const wstring& name = this->b->lnames[i];
		if ( name != L"" and rel_faults.find(name) != rel_faults.end() ) {
			g1(rel_faults[name]) = i;
		}
	}
	
	for( auto it = rel_faults.begin(); it != rel_faults.end(); ) {
		if( g1(it->second) == -1 ) {
			it = rel_faults.erase(it);
		} else {
			++it;
		}
	}
	map<wstring, vector<triangle_pt>> retfaults;
	map<wstring, vector<size_t>> extended;
	for ( auto it = rel_faults.begin(); it != rel_faults.end(); it++ ) {
		int fa = g0(it->second);
		int fb = g1(it->second);
		const line& la = this->a->lines[fa];
		const line& lb = this->b->lines[fb];
		const auto& ancha = this->a->anchored_lines;
		const auto& anchb = this->b->anchored_lines;
		
		const wstring& name = it->first;
		vector<point3> pa;
		vector<point3> pb;
		std::transform(la.begin(), la.end(), std::back_inserter(pa), [&]( const point2& p ) { return point3(gx(p), gy(p), this->a->cut); });
		std::transform(lb.begin(), lb.end(), std::back_inserter(pb), [&]( const point2& p ) { return point3(gx(p), gy(p), this->b->cut); });
		
		try {
			vector<triangle> fplane = faultplane_for_lines(pa, pb);
			
			size_t na = la.size();
			size_t nb = lb.size();
			
			// Find which of the triangles contains an extended line.
			vector<int> exts;
			if (ancha.find(std::make_pair(fa, true)) != ancha.end()) {
				exts.push_back(0);
			}
			if (ancha.find(std::make_pair(fa, false)) != ancha.end()) {
				exts.push_back(na-1);
			}
			if (anchb.find(std::make_pair(fb, true)) != anchb.end()) {
				exts.push_back(na);
			}
			if (anchb.find(std::make_pair(fb, false)) != anchb.end()) {
				exts.push_back(na + nb - 1);
			}

			for ( size_t j = 0; j < fplane.size(); j++ ) {
				for ( size_t i = 0; i < exts.size(); i++ ) {
					if ( g0(fplane[j]) == exts[i] or
					     g1(fplane[j]) == exts[i] or
					     g2(fplane[j]) == exts[i] ) {
						extended[name].push_back(j);
						break;
					}
				}
			}
			
			// Transform the triangles from idx to points or to aligned triangles that can evaluate line intersection fast.
			auto pt = [&] ( size_t n ) {
				if ( n < na ) {
					return pa[n]; 
				} else {
					return pb[n-na]; 
				}
			};
			
			std::transform(fplane.begin(), fplane.end(), std::back_inserter(retfaults[name]),
				[&] ( const triangle& t ) {
					return std::make_tuple(pt(g0(t)), pt(g1(t)), pt(g2(t))); 
				} );
			
			std::transform(retfaults[name].begin(), retfaults[name].end(), std::back_inserter(this->faults[name]),
				[&] ( const triangle_pt& t ) -> AlignedTriangle {
					return AlignedTriangle(std::make_tuple(g0(t), g1(t), g2(t))); 
				} );
		
		} catch ( const GeomodelrException& e ) {
			if ( geomodelr_verbose ) {
				string aname(this->a->name.begin(), this->a->name.end());
				string bname(this->b->name.begin(), this->b->name.end());
				string sname(name.begin(), name.end());
				std::cerr << "could not interpolate fault " << sname << " between " 
					  << aname << " and " << bname << " " << e.what() << "\n"; 
			}
		}
	}
	vector<value_f> envelopes;
	for ( auto it = this->faults.begin(); it != this->faults.end(); it++ ) {
		wstring ft = feature_types.find(it->first)->second;
		if ( ft != L"FAULT" ) {
			// Don't check other types that are not faults.
			continue;
		}
		const auto& triangles = it->second;
		for ( size_t i = 0; i < triangles.size(); i++ ) {
			const auto& tr = triangles[i];
			box trbox;
			geometry::envelope(tr.triangle, trbox);
			envelopes.push_back(std::make_tuple(trbox, it->first, i));
		}
	}
	this->faultidx = new rtree_f(envelopes.begin(), envelopes.end());
	return std::make_tuple(retfaults, extended);
}

void MatchPython::set( const pylist& match ) {
	vector<std::pair<int, int>> vmatch;
	size_t nmatch = python::len(match);
	for ( size_t i = 0; i < nmatch; i++ ) {
		int a = python::extract<int>(match[i][0]);
		int b = python::extract<int>(match[i][1]);
		vmatch.push_back(std::make_pair(a, b));
	}
	((Match *)this)->set(vmatch);
}

pylist MatchPython::get() const {
	std::set<std::pair<int, int>> ret;
	for ( auto it = this->a_to_b.begin(); it != this->a_to_b.end(); it++ ) {
		const vector<int>& m = it->second;
		for ( size_t j = 0; j < m.size(); j++ ) {
			ret.insert(std::make_pair(it->first, m[j]));
		}
	}
	pylist tret;
	for ( auto it = ret.begin(); it != ret.end(); it++ ) {
		tret.append(python::make_tuple(it->first, it->second));
	}
	return tret;
}

pylist test_faultplane_for_lines(const pylist& pyla, const pylist& pylb) {
	auto pypoint = []( const pyobject& pt ) {
		return point3(python::extract<double>(pt[0]), 
			      python::extract<double>(pt[1]), 
			      python::extract<double>(pt[2]));
	};
	auto pytovct = [pypoint]( const pylist& pyl ) {
		vector<point3> l;
		for ( int i = 0; i < python::len(pyl); i++ ) {
			l.push_back(pypoint(python::extract<pyobject>(pyl[i])));
		}
		return l;
	};
	auto vcttopy = []( const vector<triangle>& l ) {
		pylist pyl;
		for ( size_t i = 0; i < l.size(); i++ ) {
			pyl.append(python::make_tuple(g0(l[i]),g1(l[i]),g2(l[i])));
		}
		return pyl;
	};
	vector<point3> la = pytovct(pyla);
	vector<point3> lb = pytovct(pylb);
	
	vector<triangle> res = faultplane_for_lines(la, lb);
	return vcttopy(res);
}
