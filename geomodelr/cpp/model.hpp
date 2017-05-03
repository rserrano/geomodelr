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
	point2 base_point;
	point2 direction;
	vector<Section *> sections;
	vector<Match *> match;
	vector<double> cuts;
	Topography * topography;
	
	std::pair<int, double> closest_match(bool a, int a_idx, int pol_idx, const point2& pt) const;
	
	struct Possible {
		Possible(int a_match, int b_match, double a_dist, double b_dist);
		int a_match;
		int b_match;
		double a_dist;
		double b_dist;
		bool operator<( const Possible& other ) const;
		double distance( double c ) const;
	};
	
	void clear_matches();
	// Returns all the possible matches of this 2d point, given the distance is unknown.
	vector<Possible> get_candidates(size_t a_idx, const point2& pt_a, const point2& pt_b ) const;
	vector<Possible> all_closest( size_t a_idx, const point2& pt_a, const point2& pt_b ) const;
	
public:
	Model(const std::pair<double, double>& cuts_range, const point2& basepoint, const point2& direction);
	virtual ~Model();
	std::tuple<int, int, double> closest( size_t a_idx, const point2& pt_a, const point2& pt_b, double cut ) const;
	// Methods to create matches or load them from files.
	map<wstring, vector<triangle_pt>> make_matches(); // Returns the faults in global coordinates, (at least until moving plane-fault intersection to C++).
	void set_matches(const vector< std::tuple< std::tuple<wstring, wstring>, vector<std::pair<int, int>> > >& matching);
	
	// Methods to query matches.
	vector<std::tuple<wstring, double, double>> possible_closest(const point3& pt) const;
	std::pair<point2, double> model_point(const point3& pt) const;
	point3 inverse_point(const point2& pt, double cut) const;
	std::tuple<wstring, double> closest(const point3& pt) const;
	std::tuple<wstring, double> closest_topo(const point3& pt) const;
	double height(const point2& pt) const;
};

class ModelPython : public Model {
public:
	ModelPython(const pyobject& cuts_range,
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
	pytuple model_point(const pyobject& pt) const;
	pytuple inverse_point(const pyobject& pt) const;
	pytuple closest(const pyobject& pt) const;
	pytuple closest_topo(const pyobject& pt) const;
	pydict info() const;
	double height(const pyobject& pt) const;
};

#endif
