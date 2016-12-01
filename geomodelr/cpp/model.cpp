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

Model::Model( const pyobject& base_point, const pyobject& direction, const pyobject& map, const pyobject& topography, const pylist& sections):
base_point(python::extract<double>(base_point[0]), python::extract<double>(base_point[1]) ),
direction(python::extract<double>(direction[0]), python::extract<double>(direction[1]) ),
topography(nullptr)
{
	if ( python::len(topography) > 0 ) {
		const pyobject& point = python::extract<pyobject>(topography["point"]);
		const pyobject& sample = python::extract<pyobject>(topography["sample"]);
		const pyobject& dims = python::extract<pyobject>(topography["dims"]);
		const pylist& heights = python::extract<pylist>(topography["heights"]);
		this->topography = new Topography(point, sample, dims,heights);
	}
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
	Match::load_triangulate();
}

Model::~Model(){
	if ( this->topography != nullptr ) {
		delete this->topography;
	}
	this->clear_matches();
	for ( auto it = this->sections.begin(); it != this->sections.end(); it++ ) {
		delete *it;
	}
}
void Model::clear_matches() {
	for ( Match * c: this->match ) {
		delete c;
	}
	this->match.clear();
}

pydict Model::make_matches() {
	this->clear_matches();
	map<wstring, vector<triangle_pt>> faults;
	for ( size_t i = 1; i < this->sections.size(); i++ ) {
		this->match.push_back(new Match(this->sections[i-1], this->sections[i]));
		this->match.back()->match_polygons();
		map<wstring, vector<triangle_pt>> m = this->match.back()->match_lines();
		for ( auto it = m.begin(); it != m.end(); it++ ) {
			auto& f = faults[it->first];
			f.reserve(f.size() + it->second.size());
			f.insert(f.end(), it->second.begin(), it->second.end());
		}
	}
	auto point_coords = [&]( const point3& pt ) {
		point2 p(gx(pt), gy(pt));
		point3 glp = this->to_inverse_point(p, gz(pt));
		return python::make_tuple(gx(glp), gy(glp), gz(glp));
	};
	auto triangle_coords = [point_coords]( const triangle_pt& tr ) {
		return python::make_tuple(point_coords(g0(tr)), point_coords(g1(tr)), point_coords(g2(tr)));
	};
	// Now convert to python and return.
	pydict ret;
	for ( auto it = faults.begin(); it != faults.end(); it++ ) {
		pylist tris;
		for ( const triangle_pt& t: it->second ) {
			tris.append( triangle_coords( t ) );
		}
		ret[it->first] = tris;
	}
	return ret;
}

void Model::set_matches( const pylist& matching ) {
	this->clear_matches();
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
		this->match.push_back(new Match(this->sections[i-1], this->sections[i]));
		this->match.back()->load_polygons_match(m);
	}
}

pylist Model::get_matches( ) const {
	pylist ret;
	for ( size_t i = 0; i < this->match.size(); i++ ) {
		const wstring& name1 = this->sections[i]->name;
		const wstring& name2 = this->sections[i+1]->name;
		ret.append(python::make_tuple(python::make_tuple(name1, name2), this->match[i]->get()));
	}
	return ret;
}

std::pair<int, double> Model::closest_match( bool a, int a_idx, int pol_idx, const point2& pt ) const {
	const Match& m = *(this->match[a_idx]);
	const Section& s = (a) ? *(this->sections[a_idx+1]) : *(this->sections[a_idx]);
	const vector<bool>& fr = (a) ? m.a_free : m.b_free;
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
	if ( this->match.size()+1 != this->sections.size() ) {
		throw GeomodelrException("You need to call make_matches before using this function.");
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
		std::pair<int, double> cls = s.closest_to(mp.first, geometry::index::satisfies(always_true));
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

pydict Model::info() const {
	pydict ret;
	for ( const Section * s: this->sections ) {
		ret[s->name] = s->info();
	}
	return ret;
}

pytuple Model::closest(const pyobject& pypt) const {
	
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& p2 = python::extract<double>(pypt[2]);
	point3 pt(p0, p1, p2);
	if ( this->topography != nullptr ) {
		if ( this->topography->height(pt) < gz(pt) ) {
			return python::make_tuple(wstring(L"AIR"), std::numeric_limits<double>::infinity());
		}
	}
	if ( not this->sections.size() ) {
		return python::make_tuple(wstring(L"NONE"), std::numeric_limits<double>::infinity());
	}
	
	if ( this->match.size()+1 != this->sections.size() ) {
		throw GeomodelrException("You need to call make_matches before using this function.");
	}
	
	std::pair<point2, double> mp = this->to_model_point(pt);
	auto it = std::upper_bound(this->cuts.begin(), this->cuts.end(), mp.second);
	size_t a_idx = it - this->cuts.begin();
	// If it's behind the last or above the first, return the closest in the section.

	auto closest_single =[&](const Section& s) {
		std::pair<int, double> cls = s.closest_to(mp.first, geometry::index::satisfies(always_true));
		if ( cls.first == -1 ) {
			return python::make_tuple(wstring(L"NONE", cls.second));
		}
		wstring unit = s.units[cls.first];
		return python::make_tuple(unit, cls.second);
	};
	// For a cut below the lowest or above the highest.
	if ( a_idx <= 0 or a_idx >= this->sections.size() ) {
		const Section& s = ( a_idx <=0 ) ? *(this->sections.front()) : *(this->sections.back());
		return closest_single(s);
	}
	// Check if it crosses a fault and then only evaluate the cross in the actual site.
	a_idx--;
	int crosses = this->match[a_idx]->crosses_triangles(mp.first, mp.second);
	if ( crosses < 0 ) {
		return closest_single(*(this->sections[a_idx]));
	} else if ( crosses > 0 ) {
		return closest_single(*(this->sections[a_idx+1]));
	}
	// Finally evaluate the full transition.
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

Topography::Topography( const pyobject& point, const pyobject& sample, const pyobject& dims, const pylist& heights ):
point(python::extract<double>(point[0]), python::extract<double>(point[1])),
sample(python::extract<double>(sample[0]), python::extract<double>(sample[1])),
heights(python::extract<int>(dims[0]) * python::extract<int>(dims[1]))
{
	this->dims[0] = python::extract<size_t>(dims[0]);
	this->dims[1] = python::extract<size_t>(dims[1]);
	size_t rows = python::len(heights);
	if ( rows != this->dims[0] ) {
		throw GeomodelrException("topography rows does not correspond with dims.");
	}
	for ( size_t i = 0; i < rows; i++ ) {
		size_t cols = python::len(heights[i]);
		if ( cols != this->dims[1] ) {
			throw GeomodelrException("topography columns does not correspond with dims.");
		}
		for ( size_t j = 0; j < cols; j++ ) {
			this->heights[i*this->dims[1]+j] = python::extract<double>(heights[i][j]);
		}
	}
}

double Topography::height(const point3& pt) const {
	point2 pos(gx(pt), gy(pt));
	geometry::subtract_point(pos, this->point);
	geometry::divide_point(pos, this->sample);
	size_t x = int(gx(pos));
	size_t y = int(gy(pos));

        if ( x < 0 )
            x = 0;
        
        if ( y < 0 )
            y = 0;
        
        if ( x >= this->dims[0] )
            x = this->dims[0]-1;
        
        if ( y >= this->dims[1] )
            y = this->dims[1]-1;
        
	return this->heights[x*dims[1]+y];
}

