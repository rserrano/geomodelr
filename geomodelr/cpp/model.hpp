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
#ifndef GEOMODELR_MODEL_HPP
#define GEOMODELR_MODEL_HPP

#include "basic.hpp"
#include "section.hpp"
#include "match.hpp"

class Topography {
protected:
	point2 point;
	point2 sample;
	int dims[2];
	vector<double> heights;
public:
	double height(const point2&) const;
	Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims );
	Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims, const vector<double>& heights );
};

class TopographyPython : public Topography {
public:
	double height(const pyobject&) const;
	TopographyPython( const pyobject& point, const pyobject& sample, const pyobject& dims, const pylist& heights );
};

class Model {
protected:
	std::pair<double, double> cuts_range;
	std::tuple<std::tuple<double, double, double>, 
		   std::tuple<double, double, double>> bbox;
	
	point2 base_point;
	point2 direction;
	vector<Section *> sections;
	vector<Match *> match;
	vector<double> cuts;
	Topography * topography;
	
	std::pair<int, double> closest_match(bool a, int a_idx, int pol_idx, const point2& pt) const;
	
	struct Possible {
		Possible(int a_match, int b_match, double a_dist, double b_dist);
		// The match in the section a.
		int a_match;
		// The match in the section b.
		int b_match;
		// The polygon distance in section a.
		double a_dist;
		// The polygon distance in section b.
		double b_dist;
		// Order to be able to put in a map.
		bool operator<( const Possible& other ) const;
		// The given distance, at the given c.
		double distance( double c ) const;
	};
	
	void clear_matches();
	// Returns all the possible matches of this 2d point, given the distance is unknown.
	template<typename Predicates>
	vector<Possible> get_candidates( size_t a_idx, const point2& pt_a, const point2& pt_b, const Predicates& predicates ) const {
		const Section& s_a = *(this->sections[a_idx]);
		const Section& s_b = *(this->sections[a_idx+1]);
		
		double s_dist = this->cuts[a_idx+1]-this->cuts[a_idx];
		
		// Given a match in section a, find a match in section b
		double inf = std::numeric_limits<double>::infinity();
		auto a_possible = [&](const std::pair<int, double>& cls) {
			// If it didn't find a match, just return no match.
			if ( cls.first == -1 ) {
				return Possible(-1, -1, inf, inf);
			}
			// Find the closest match in the oposite cross section.
			std::pair<int, double> match = this->closest_match( true,  a_idx, cls.first, pt_b );
			if ( match.first != -1 ) {
				return Possible(cls.first, match.first, cls.second, match.second);
			} else {
				return Possible(cls.first, -1, cls.second, cls.second+s_dist);
			}
		};
		
		// Given a match in section b, find a match in section a.
		auto b_possible = [&](const std::pair<int, double>& cls) { 
			if ( cls.first == -1 ) {
				return Possible(-1, -1, inf, inf);
			}
			// Find the match in section a.
			std::pair<int, double> match = this->closest_match( false,  a_idx, cls.first, pt_a );
			if ( match.first != -1 ) {
				return Possible(match.first, cls.first, match.second, cls.second);
			} else {
				return Possible(-1, cls.first, cls.second+s_dist, cls.second);
			}
		};
		
		std::set<Possible> possible;
		// Find the closest polygon to a.
		auto closest_in_a = s_a.closest(pt_a, predicates);
		
		// Find the closest polygon in s_a => clst_a.
		Possible clst_a = a_possible(closest_in_a);
		
		// Find all polygons closer than the one we just found.
		if ( clst_a.a_match != -1 ) {
			possible.insert(clst_a);
			auto not_first = [&](const value_f& v) { return g2(v) != clst_a.b_match; }; 
			// Find the pols in s_b that are closer mindist => all_cls_a.
			vector<std::pair<int, double>> clsr_b = s_b.closer_than( pt_b, clst_a.b_dist, add_funcs<decltype(not_first), decltype(predicates)>(not_first, predicates) );
			if ( not clsr_b.size() ) 
			{
				return vector<Possible>{ clst_a };
			}
			for ( const auto& c: clsr_b ) {
				possible.insert(b_possible(c));
			}
		}
		
		// Repeat the operation from s_b to s_a.
		auto closest_in_b = s_b.closest(pt_b, predicates);
		Possible clst_b = b_possible(closest_in_b);
		
		if ( clst_b.b_match != -1 ) {
			possible.insert(clst_b);
			auto not_first = [&](const value_f& v) { return g2(v) != clst_b.a_match; }; // CHECK it's b_match.
			// Find the pols in s_b that are closer mindist => all_cls_a.
			vector<std::pair<int, double>> clsr_a = s_a.closer_than( pt_a, clst_b.a_dist, add_funcs<decltype(not_first), decltype(predicates)>(not_first, predicates) ); // CHECK it's pt_b
			for ( const auto& c: clsr_a ){
				possible.insert(a_possible(c));
			}
		}
	
		return vector<Possible>(possible.begin(), possible.end());
	}
	
	
	vector<Possible> all_closest( size_t a_idx, const point2& pt_a, const point2& pt_b ) const;
	
	template<typename Predicates>
	std::tuple<int, int, double> closest_between( size_t a_idx, const point2& pt_a, const point2& pt_b, double cut, const Predicates& predicates ) const {
		/* Returns the closest formation to a model. */
		double s_dist = this->cuts[a_idx+1]-this->cuts[a_idx];
		double cut_rem = cut - this->cuts[a_idx];
		double mult_b = cut_rem/s_dist;
		// assert 0.0 <= mult_b <= 1.0
		double mult_a = 1.0-mult_b;
		
		vector<Possible> possible = this->get_candidates(a_idx, pt_a, pt_b, predicates);
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
	std::tuple<int, int, double> closest_between( size_t a_idx, const point2& pt_a, const point2& pt_b, double cut ) const;
public:
	
	Model(const std::tuple<std::tuple<double, double, double>, std::tuple<double, double, double>>& bbox, const point2& basepoint, const point2& direction);
	virtual ~Model();
	
	double signed_distance( const wstring& unit, const point3& pt ) const;
	// In this case the bounding box bounds all the solids.
	double signed_distance_bounded( const wstring& unit, const point3& pt ) const;
	// In this case the solids are not bounded by the bounding box, only by the topography.
	double signed_distance_unbounded( const wstring& unit, const point3& pt ) const;

	// Methods to create matches or load them from files.
	map<wstring, vector<triangle_pt>> make_matches(); // Returns the faults in global coordinates, (at least until moving plane-fault intersection to C++).
	void set_matches(const vector< std::tuple< std::tuple<wstring, wstring>, vector<std::pair<int, int>> > >& matching);
	
	// Methods to query matches.
	vector<std::tuple<wstring, double, double>> possible_closest(const point3& pt) const;
	std::pair<point2, double> model_point(const point3& pt) const;
	point3 inverse_point(const point2& pt, double cut) const;
	template<typename Predicates>
	std::tuple<wstring, double> closest( const point3& pt, const Predicates& predicates ) const {
		
		if ( not this->sections.size() ) {
			return std::make_tuple(wstring(L"NONE"), std::numeric_limits<double>::infinity());
		}
		
		if ( this->match.size()+1 != this->sections.size() ) {
			throw GeomodelrException("You need to call make_matches before using this function.");
		}
		
		std::pair<point2, double> mp = this->model_point(pt);
		
		auto it = std::upper_bound(this->cuts.begin(), this->cuts.end(), mp.second);
		
		size_t a_idx = it - this->cuts.begin();
		// If it's behind the last or above the first, return the closest in the section.
		
		auto closest_single = [&](const Section& s) {
			std::pair<int, double> cls = s.closest(mp.first, predicates);
			if ( cls.first == -1 ) {
				return std::make_tuple(wstring(L"NONE"), cls.second);
			}
			wstring unit = s.units[cls.first];
			return std::make_tuple(unit, cls.second);
		};
		// For a cut below the lowest or above the highest.
		if ( a_idx <= 0 ) {
			return closest_single(*this->sections.front());
		}
		if ( a_idx  >= this->sections.size() ) {
			return closest_single(*this->sections.back());
		}
		auto closest_middle = [&]( const point2& pt_a, const point2& pt_b ) {
			// Finally evaluate the full transition.
			std::tuple<int, int, double> clst = this->closest_between(a_idx, pt_a, pt_b, mp.second, predicates);
			int a_match = std::get<0>(clst);
			int b_match = std::get<1>(clst);
			wstring unit;
			
			// Get the closest unit.
			if ( a_match == -1 and b_match == -1 ) {
				unit = L"NONE";
			} else if ( a_match == -1 ) { // If there's no unit located in one cross section. return the other.
				unit = this->sections[a_idx+1]->units[b_match];
			} else { // There's either a unit in the cross section or it's the same in both.
				unit = this->sections[a_idx]->units[a_match];
			}
			return std::make_tuple(unit, std::get<2>(clst));
		};
	
	
		// Check if it crosses a fault and then, evaluate the cross section in the side normally, but move the point to the fault in the other.
		a_idx--;
		std::tuple<int, int, int> crosses = this->match[a_idx]->crosses_triangles(mp.first, mp.second);
		if ( g0(crosses) < 0 ) {
			std::tuple<point2, double> closest_in_line = point_line_projection( mp.first, this->sections[a_idx+1]->lines[g2(crosses)] );
			return closest_middle( mp.first, g0(closest_in_line) );
		} else if ( g0(crosses) > 0 ) {
			std::tuple<point2, double> closest_in_line = point_line_projection( mp.first, this->sections[a_idx]->lines[g1(crosses)] );
			return closest_middle( g0(closest_in_line), mp.first );
		}
		return closest_middle( mp.first, mp.first );
	}
	std::tuple<wstring, double> closest(const point3& pt) const;
	
	std::tuple<wstring, double> closest_topo(const point3& pt) const;
	double height(const point2& pt) const;
};

class ModelPython : public Model {
public:
	ModelPython(const pyobject& bbox,
	            const pyobject& basepoint, 
	            const pyobject& direction, 
	            const pyobject& map, 
	            const pyobject& topography,
	            const pylist& sections);
	
	// Methods to create matches or load them from files.
	pydict make_matches(); // Returns the faults in global coordinates, (at least until moving plane-fault intersection to C++).
	void set_matches(const pylist& matching);
	pylist get_matches() const;
	// Methods to query matches.
	pylist possible_closest(const pyobject& pt) const;
	pylist pybbox() const;
	pytuple model_point(const pyobject& pt) const;
	pytuple inverse_point(const pyobject& pt) const;
	pytuple closest(const pyobject& pt) const;
	pytuple closest_topo(const pyobject& pt) const;
	double signed_distance( const wstring& unit, const pyobject& pt ) const;
	double signed_distance_bounded( const wstring& unit, const pyobject& pt ) const;
	double signed_distance_unbounded( const wstring& unit, const pyobject& pt ) const;
	pydict info() const;
	double height(const pyobject& pt) const;
};

#endif
