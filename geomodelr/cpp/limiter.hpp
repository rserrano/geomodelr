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

#ifndef GEOMODELR_LIMITER_HPP
#define GEOMODELR_LIMITER_HPP

#include "basic.hpp"
#include <assert.h>
class Model;
class ModelPython;

// Nodes of the binary tree of the topography
struct node {
	int i_min;
	int i_max;
	int j_min;
	int j_max;
	double z_min;
	double z_max;

	node *left;
	node *right;
};

class Topography {
	friend Model;
protected:
	point2 point;
	point2 sample;
	int dims[2];
	vector<double> heights;
	node * topography_tree;
public:
	double max;
	double min;
	double height(const point2&) const;
	double height_new(const point2& pt) const;
	std::map<double, point2> get_points_section(const point2& p0, const point2& pf);
	std::pair<point3, double> intersection(const point3& pt, const point3& projection);

	void destroy_tree(node * tree);
	void set_tree(int i_min, int i_max, int j_min, int j_max, node * leaf);
	void print_level(int k, const node * leaf, int level);
	std::tuple<int, int, int, int> square_limits(const point3& pt, const point3& projection) const;
	std::pair<point3, double> intersects(const point3& pt, const point3& projection);
	void tree_intersection(std::pair<point3, double>& output, const point3& pt, const point3& vec,
		const point3& abs_vec, const node * leaf);

	Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims );
	Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims, const vector<double>& heights );
	~Topography();
};

class TopographyPython : public Topography {
public:
	double height(const pyobject&) const;
	pylist get_points_section(const pyobject& p0_py, const pyobject& pf_py) const;
	void info();
	void print_level(int level);
	TopographyPython( const pyobject& point, const pyobject& sample, const pyobject& dims, const pylist& heights );

	pytuple inter_1(const pyobject& pt, const pyobject& pro);
	pytuple inter_2(const pyobject& pt, const pyobject& pro);
};

// Abstract limiters to send to signed distance function.

class Limiter {
public:
  virtual double limit_signed_distance(const point3& pt, double sds) const = 0;
  virtual ~Limiter() = 0;
};


// Polygon limiters, the main reason I do this.
class PolygonLimiter: public Limiter {

polygon limit;
line lpoly;
double bottom;
const Model * model;

public:
  PolygonLimiter(const polygon& poly, double bottom, const Model * model);
  virtual double limit_signed_distance(const point3& pt, double sds) const;
  virtual ~PolygonLimiter();
};
// BBox limiter, it limits to a given bbox. (it can be different from the bbox of the model).
class BBoxLimiter: public Limiter {
bbox3 bbox;
const Model * model;

public:
  BBoxLimiter(const bbox3& bbox, const Model * model);
  virtual ~BBoxLimiter();
  virtual double limit_signed_distance(const point3& pt, double sds) const;
};

class RestrictedFunction {
std::shared_ptr<Limiter> limit;
const ModelPython * model;
public:
  RestrictedFunction( const pyobject& model, const wstring& restype, const pyobject& data );
  double signed_distance( const wstring& unit, const pyobject& point ) const;
};

#endif
