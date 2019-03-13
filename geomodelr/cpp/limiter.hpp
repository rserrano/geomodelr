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
class Model;
class ModelPython;

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

// Abstract limiters to send to signed distance function.

class Limiter {
public:
  virtual double limit_signed_distance(const point3& pt, double sds) const = 0;
  virtual ~Limiter() = 0;
};

// Specific class to pass to sd aligned.
class AlignedLimiter {
public:
  virtual double limit_signed_distance(const point3& pt, double sds) const = 0;
  virtual ~AlignedLimiter() = 0;
};

// Polygon limiters, the main reason I do this.
class PolygonLimiter: public Limiter {
polygon limit;
line lpoly;
const Model * model;

public:
  PolygonLimiter(const polygon& poly, const Model * model);
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

// Class specific to aligned sdf.
class BBoxAlignedLimiter: public AlignedLimiter {
bbox3 abbox;
const Model * model;

public:
  BBoxAlignedLimiter(const bbox3& abbox, const Model * model);
  virtual ~BBoxAlignedLimiter();
  virtual double limit_signed_distance(const point3& pt, double sds) const;
};

class RestrictedFunction {
std::shared_ptr<Limiter> limit;
const ModelPython * model;
public:
  RestrictedFunction( const pyobject& model, const wstring& restype, const pyobject& data );
  double signed_distance( const wstring& unit, const pyobject& point ) const;
};

class AlignedRestrictedFunction {
std::shared_ptr<AlignedLimiter> limit;
const ModelPython * model;
public:
  AlignedRestrictedFunction( const pyobject& model, const wstring& restype, const pyobject& data );
  double signed_distance( const wstring& unit, const pyobject& point ) const;
};

#endif
