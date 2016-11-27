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

#include "geomodel.hpp"
#include <algorithm>
#include <functional>
#include <cmath>

/*
template<typename P>
void print_possible( const P& possible ) {
	for ( const auto& a: possible ) {
		std::cerr << "p " << a.a_match << " " << a.b_match << " " << a.a_dist << " " << a.b_dist << std::endl;
	}
	std::cerr << "\n";
}*/

Match::Match( const vector<std::pair<int, int>>& match, size_t sa, size_t sb ): a_free(sa, true), b_free(sb, true) {
	for ( size_t i = 0; i < match.size(); i++ ) {
		this->a_free[match[i].first]  = false;
		this->b_free[match[i].second] = false;
		a_to_b[match[i].first].push_back(match[i].second);
		b_to_a[match[i].second].push_back(match[i].first);
	}
}

pylist Match::get() const {
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

Match Model::make_match( const Section& a, const Section& b ) {
	map<wstring, vector<int>> units_a;
	map<wstring, vector<int>> units_b;
	for ( size_t i = 0; i < a.polygons.size(); i++ ) {
		units_a[a.units[i]].push_back(i);
	}
	for ( size_t i = 0; i < b.polygons.size(); i++ ) {
		units_b[b.units[i]].push_back(i);
	}
	vector<std::pair<int, int>> m;
	for ( auto it = units_a.begin(); it != units_a.end(); it++ ) {
		if ( units_b.find(it->first) != units_b.end() ) {
			vector<int>& pols_a = it->second;
			vector<int>& pols_b = units_b[it->first];
			for ( size_t i = 0; i < pols_a.size(); i++ )
			{
				for ( size_t j = 0; j < pols_b.size(); j++ ) {
					if ( geometry::intersects(a.polygons[pols_a[i]], b.polygons[pols_b[j]]) ) {
						m.push_back(std::make_pair(pols_a[i], pols_b[j]));
					}
				}
			}
		}
	}
	return Match( m, a.polygons.size(), b.polygons.size() );
}

Match Model::load_match( const pylist& match, size_t sa, size_t sb ) {
	vector<std::pair<int, int>> vmatch;
	size_t nmatch = python::len(match);
	for ( size_t i = 0; i < nmatch; i++ ) {
		int a = python::extract<int>(match[i][0]);
		int b = python::extract<int>(match[i][1]);
		vmatch.push_back(std::make_pair(a, b));
	}
	return Match(vmatch, sa, sb);
}

Model::Model( const pyobject& base_point, const pyobject& direction, const pylist& sections):
base_point(python::extract<double>(base_point[0]), python::extract<double>(base_point[1]) ),
direction(python::extract<double>(direction[0]), python::extract<double>(direction[1]) )
{
	size_t nsects = python::len(sections);
	for ( size_t i = 0; i < nsects; i++ ) {
		const wstring& name    = python::extract<wstring>(sections[i][0]);
		double cut             = python::extract<double>(sections[i][1]);
		const pylist& points   = python::extract<pylist>(sections[i][2]);
		const pylist& polygons = python::extract<pylist>(sections[i][3]);
		const pylist& units    = python::extract<pylist>(sections[i][4]); 
		const pylist& lines    = python::extract<pylist>(sections[i][5]);
		const pylist& lnames   = python::extract<pylist>(sections[i][6]);
		this->sections.push_back(new Section(name, cut, points, polygons, units, lines, lnames));
	}
	std::sort(this->sections.begin(), this->sections.end(), [](const Section* a, const Section* b){ return a->cut < b->cut; });
	for ( size_t i = 0; i < this->sections.size(); i++ ) {
		this->cuts.push_back(this->sections[i]->cut);
	}
}

Model::~Model(){
	for ( auto it = this->sections.begin(); it != this->sections.end(); it++ ) {
		delete *it;
	}
}

void Model::make_matches() {
	this->match.clear();
	for ( size_t i = 1; i < this->sections.size(); i++ ) {
		this->match.push_back(this->make_match(*(this->sections[i-1]), *(this->sections[i])));
	}
}

void Model::set_matches( const pylist& matching ) {
	this->match.clear();
	size_t nmatch = python::len(matching);
	map<std::pair<wstring, wstring>, int> match;
	for ( size_t i = 0; i < nmatch; i++ ) {
		const wstring& name1 = python::extract<wstring>(matching[i][0][0]);
		const wstring& name2 = python::extract<wstring>(matching[i][0][1]);
		auto key = std::make_pair(name1,name2);
		match[key] = i;
	}
	for ( size_t i = 1; i < this->sections.size(); i++ ) {
		const wstring& name1 = this->sections[i-1]->name;
		const wstring& name2 = this->sections[i]->name;
		if ( match.find(std::make_pair(name1, name2)) == match.end() ) {
			string sname1(name1.begin(), name1.end());
			string sname2(name2.begin(), name2.end());
			string error("match does not contain sections: ");
			error += sname1;
			error += ", ";
			error += sname2;
			throw GeomodelrException(error.c_str());
		}
		const pylist& m = python::extract<pylist>(matching[match[std::make_pair(name1, name2)]][1]);
		this->match.push_back(this->load_match(m, this->sections[i-1]->polygons.size(), this->sections[i]->polygons.size()));
	}
}

pylist Model::get_matches( ) const {
	pylist ret;
	for ( size_t i = 0; i < this->match.size(); i++ ) {
		const wstring& name1 = this->sections[i]->name;
		const wstring& name2 = this->sections[i+1]->name;
		ret.append(python::make_tuple(python::make_tuple(name1, name2), this->match[i].get()));
	}
	return ret;
}

std::pair<int, double> Model::closest_match( bool a, int a_idx, int pol_idx, const point2& pt ) const {
	const Match& m = this->match[a_idx];
	const Section& s = (a) ? *(this->sections[a_idx+1]) : *(this->sections[a_idx]);
	const vector<bool> fr = (a) ? m.a_free : m.b_free;
	if ( fr[pol_idx] ) {
		return std::make_pair(-1, std::numeric_limits<double>::infinity());
	}
	const vector<int>& op = (a) ? m.a_to_b.find(pol_idx)->second: m.b_to_a.find(pol_idx)->second;
	
	double mindist = std::numeric_limits<double>::infinity();
	int minidx = -1;
	for ( size_t i = 0; i < op.size(); i++ ) {
		double dist = geometry::distance(s.polygons[op[i]], pt);
		if ( dist < mindist ) {
			mindist = dist;
			minidx = op[i];
		}
	}
	return std::make_pair( minidx, mindist );
}

Model::Possible::Possible(int a_match, int b_match, double a_dist, double b_dist):
a_match(a_match), b_match(b_match), a_dist(a_dist), b_dist(b_dist) 
{
}

bool Model::Possible::operator<(const Model::Possible& other) const
{
	if ( this->a_match < other.a_match ) {
		return true;
	}
	return this->b_match < other.b_match;
}
double Model::Possible::distance(double c) const {
	return this->a_dist*c + this->b_dist*(1.0-c);
}

vector<Model::Possible> Model::get_candidates(size_t a_idx, const point2& pt ) const {
	using geometry::index::satisfies;
	const Section& s_a = *(this->sections[a_idx]);
	const Section& s_b = *(this->sections[a_idx+1]);
	
	auto always_true = [](const value& v){ return true; };
	double s_dist = this->cuts[a_idx+1]-this->cuts[a_idx];
	
	// Simple lambda for obtaining the match as a Possible.
	double inf = std::numeric_limits<double>::infinity();
	auto a_possible = [&](const std::pair<int, double>& cls) { 
		if ( cls.first == -1 ) {
			return Model::Possible(-1, -1, inf, inf);
		}
		std::pair<int, double> match = this->closest_match( true,  a_idx, cls.first, pt );
		if ( match.first != -1 ) {
			return Model::Possible(cls.first, match.first, cls.second, match.second);
		} else {
			return Model::Possible(cls.first, -1, cls.second, cls.second+s_dist);
		}
	};
	
	auto b_possible = [&](const std::pair<int, double>& cls) { 
		if ( cls.first == -1 ) {
			return Model::Possible(-1, -1, inf, inf);
		}
		std::pair<int, double> match = this->closest_match( false,  a_idx, cls.first, pt );
		if ( match.first != -1 ) {
			return Model::Possible(match.first, cls.first, match.second, cls.second);
		} else {
			return Model::Possible(-1, cls.first, cls.second+s_dist, cls.second);
		}
	};
	
	std::set<Model::Possible> possible;
	
	// Find the closest polygon in s_a => clst_a.
	Model::Possible clst_a = a_possible(s_a.closest_to(pt, satisfies(always_true)));
	
	if ( clst_a.a_match != -1 ) {
		possible.insert(clst_a);
		auto not_first = [&](const value& v) { return v.second != clst_a.b_match; };
		// Find the pols in s_b that are closer mindist => all_cls_a
		vector<std::pair<int, double>> clsr_b = s_b.closer_than( pt, clst_a.b_dist, satisfies(not_first) );
		if ( not clsr_b.size() ) 
		{
			return vector<Model::Possible>{ clst_a };
		}
		for ( const auto& c: clsr_b ) {
			possible.insert(b_possible(c));
		}
	}
	
	// Repeat the operation from s_b to s_a.
	Model::Possible clst_b = b_possible(s_b.closest_to(pt, satisfies(always_true)));
	
	if ( clst_b.b_match != -1 ) {
		possible.insert(clst_b);
		auto not_first = [&](const value& v) { return v.second != clst_b.b_match; };
		// Find the pols in s_b that are closer mindist => all_cls_a.
		vector<std::pair<int, double>> clsr_a = s_a.closer_than( pt, clst_b.a_dist, satisfies(not_first) );
		for ( const auto& c: clsr_a ){
			possible.insert(a_possible(c));
		}
	}
	return vector<Model::Possible>(possible.begin(), possible.end());
}

vector<Model::Possible> Model::all_closest( size_t a_idx, const point2& pt ) const {
	vector<Model::Possible> possible = this->get_candidates(a_idx, pt);
	// First drop the non possible.
	std::set<double> intersections;
	vector<bool> drop(possible.size(), false);
	for ( size_t idx = 0; idx < possible.size(); idx++ ) {
		if ( drop[idx] ) {
			continue;
		}
		for ( size_t jdx = idx+1; jdx < possible.size(); jdx++ ) {
			if ( drop[jdx] ) {
				continue;
			}
			
			double da1 = possible[idx].a_dist;
			double db1 = possible[idx].b_dist;
			double da2 = possible[jdx].a_dist;
			double db2 = possible[jdx].b_dist; 
			
			if ( da1 < da2 and db1 < db2 ) {
				drop[jdx] = true;
			}
			if ( da1 >= da2 and db1 >= db2 ) {
				drop[idx] = true;
				break;
			}
		}
	}
	// Then get the intersections for the possible.
	for ( size_t idx = 0; idx < possible.size(); idx++ ) {
		if ( drop[idx] ) {
			continue;
		}
		for ( size_t jdx = idx+1; jdx < possible.size(); jdx++ ) {
			if ( drop[jdx] ) {
				continue;
			}
			
			double da1 = possible[idx].a_dist;
			double db1 = possible[idx].b_dist;
			double da2 = possible[jdx].a_dist;
			double db2 = possible[jdx].b_dist; 
			
			double da = da1-da2;
			double db = db1-db2;
			
			intersections.insert(da/(da-db));
		}
	}
	vector<Model::Possible> reduced;
	for ( size_t idx = 0; idx < possible.size(); idx++ ) {
		if ( not drop[idx] ) {
			reduced.push_back(possible[idx]);
		}
	}
	intersections.insert(0.0);
	intersections.insert(1.0);
	vector<double> mids;
	for ( auto it = intersections.begin(); std::next(it) != intersections.end(); it++ ) {
		mids.push_back((*it + *std::next(it))/2.0);
	}
	vector<std::pair<int, double>> minimum(intersections.size(), std::make_pair(-1, std::numeric_limits<double>::infinity()));
	for ( size_t i = 0; i < reduced.size(); i++ ) {
		int jdx = 0;
		for ( auto jt = intersections.begin(); jt != intersections.end(); jt++ ){
			double d = reduced[i].distance(*jt);
			if ( d < minimum[jdx].second ) {
				minimum[jdx].first = i;
				minimum[jdx].second = d;
			}
			jdx++;
		}
	}
	vector<bool> final_possible(reduced.size(), false);
	for ( size_t i = 0; i < minimum.size(); i++ ) {
		if ( minimum[i].first != -1 ) {
			final_possible[minimum[i].first] = true;
		}
	}
	vector<Model::Possible> ret;
	for ( size_t i = 0; i < reduced.size(); i++ ){
		if ( final_possible[i] ) {
			ret.push_back(reduced[i]);
		}
	}
	return ret;
}

std::tuple<int, int, double> Model::closest_to( size_t a_idx, const point2& pt, double cut ) const {
	double s_dist = this->cuts[a_idx+1]-this->cuts[a_idx];
	double cut_rem = cut - this->cuts[a_idx];
	double mult_b = cut_rem/s_dist;
	double mult_a = 1.0-mult_b;
	
	vector<Model::Possible> possible = this->get_candidates(a_idx, pt);
	double mindist = std::numeric_limits<double>::infinity();
	double minidx = -1;
	
	for ( size_t i = 0; i < possible.size(); i++ ) {
		double dist = possible[i].distance(mult_a);
		if ( dist < mindist ) {
			mindist = dist;
			minidx = i;
		}
	}
	if ( minidx == -1 ) {
		return std::make_tuple(-1, -1, std::numeric_limits<double>::infinity());
	}
	return std::make_tuple(possible[minidx].a_match, possible[minidx].b_match, mindist);
}


std::pair<point2, double> Model::to_model_point( const point3& pt ) const {

	point2 ppoint(gx(pt), gy(pt));
	const double& z = gz(pt);
	geometry::subtract_point(ppoint, this->base_point);
	
	double norm = std::sqrt(gx(ppoint)*gx(ppoint) + gy(ppoint)*gy(ppoint));

	if ( norm < 1e-10 ) {
		return std::make_pair(point2(0.0, z), 0.0);
	}
	
	geometry::divide_value(ppoint, norm);

	const double& d0 = gx(this->direction);
	const double& d1 = gy(this->direction);

	double c = d0*gy(ppoint) - d1*gx(ppoint);
	double d = d0*gx(ppoint) + d1*gy(ppoint);
	return std::make_pair(point2(d*norm, z), c*norm);
}

point3 Model::to_inverse_point( const point2& pt, double cut ) const {
	const double& d0 = gx(this->direction);
	const double& d1 = gy(this->direction);
	
	double dsq = d0*d0 + d1*d1;
	
	const double& p0 = gx(pt);
	const double& p1 = gy(pt);
	
        double v1 = (cut*d0 + p0*d1)/dsq;
        double v0 = ( std::fabs(d1) < 1e-9 ) ? p0/d0 : (d0*v1-cut)/d1;
        
	return point3(gx(this->base_point) + v0, gy(this->base_point) + v1, p1);
}

pytuple Model::model_point( const pyobject& pypt ) const {
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& p2 = python::extract<double>(pypt[2]);
	point3 pt(p0, p1, p2);
	std::pair<point2, double> mp = this->to_model_point( pt );
	const double& m0 = gx(mp.first);
	const double& m1 = gy(mp.first);
	return python::make_tuple(m0, m1, mp.second);
}

pytuple Model::inverse_point( const pyobject& pypt ) const {
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& cut = python::extract<double>(pypt[2]);
	point2 pt(p0, p1);
	point3 mp = this->to_inverse_point( pt, cut );
	const double& m0 = gx(mp);
	const double& m1 = gy(mp);
	const double& m2 = gz(mp);
	return python::make_tuple(m0, m1, m2);
}

pylist Model::possible_closest(const pyobject& pypt) const {
	if ( not this->sections.size() ) {
		return pylist();
	}
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& p2 = python::extract<double>(pypt[2]);
	point3 pt(p0, p1, p2);
	std::pair<point2, double> mp = this->to_model_point(pt);
	auto it = std::upper_bound(this->cuts.begin(), this->cuts.end(), mp.second);
	size_t a_idx = it - this->cuts.begin();
	// If it's behind the last or above the first, return the closest in the section.
	if ( a_idx <= 0 or a_idx >= this->sections.size() ) {
		const Section& s = ( a_idx <=0 ) ? *(this->sections.front()) : *(this->sections.back());
		std::pair<int, double> cls = s.closest_to(mp.first, geometry::index::satisfies(s.valid));
		wstring unit = s.units[cls.first];
		pylist ret;
		ret.append(python::make_tuple(unit, cls.second, cls.second));
		return ret;
	}
	a_idx--;
	vector<Model::Possible> ap = this->all_closest(a_idx, mp.first);
	pylist ret;
	for ( size_t i = 0; i < ap.size(); i++ ) {
		int a_match = ap[i].a_match;
		int b_match = ap[i].b_match;
		wstring unit = ( a_match != -1 ) ? this->sections[a_idx]->units[a_match] : this->sections[a_idx+1]->units[b_match];
		ret.append(python::make_tuple(unit, ap[i].a_dist, ap[i].b_dist));
	}
	return ret;
}

pytuple Model::closest(const pyobject& pypt) const {
	if ( not this->sections.size() ) {
		return python::make_tuple(wstring(L"NONE"), std::numeric_limits<double>::infinity());
	}
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& p2 = python::extract<double>(pypt[2]);
	point3 pt(p0, p1, p2);
	std::pair<point2, double> mp = this->to_model_point(pt);
	auto it = std::upper_bound(this->cuts.begin(), this->cuts.end(), mp.second);
	size_t a_idx = it - this->cuts.begin();
	// If it's behind the last or above the first, return the closest in the section.
	if ( a_idx <= 0 or a_idx >= this->sections.size() ) {
		const Section& s = ( a_idx <=0 ) ? *(this->sections.front()) : *(this->sections.back());
		std::pair<int, double> cls = s.closest_to(mp.first, geometry::index::satisfies(s.valid));
		wstring unit = s.units[cls.first];
		return python::make_tuple(unit, cls.second);
	}
	a_idx--;
	std::tuple<int, int, double> clst = this->closest_to(a_idx, mp.first, mp.second);
	int a_match = std::get<0>(clst);
	int b_match = std::get<1>(clst);
	wstring unit;
	if ( a_match == -1 and b_match == -1 ) {
		unit = L"NONE";
	} else if ( a_match == -1 ) {
		unit = this->sections[a_idx+1]->units[b_match];
	} else {
		unit = this->sections[a_idx]->units[a_match];
	}
	return python::make_tuple(unit, std::get<2>(clst));
}
