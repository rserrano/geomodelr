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

#include "model.hpp"
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <array>

Model::Model( const std::tuple<std::tuple<double,double,double>, std::tuple<double,double,double>>& bbox, const point2& base_point, const point2& direction )
: bbox(bbox), base_point(base_point), direction(direction), topography(nullptr), horizontal(false)
{
	point2 b0(g0(g0(bbox)), g1(g0(bbox)));
	point2 b1(g0(g1(bbox)), g1(g1(bbox)));
	point2 b2(g0(g0(bbox)), g1(g1(bbox)));
	point2 b3(g0(g1(bbox)), g1(g0(bbox)));
	auto cut = [base_point, direction]( point2& v ) {
		geometry::subtract_point(v, base_point);
		double norm = std::sqrt(gx(v)*gx(v) + gy(v)*gy(v));
		if ( norm < tolerance ) {
			return 0.0;
		}
		geometry::divide_value(v, norm);
		double s = gx(direction)*gy(v) - gy(direction)*gx(v);
		return norm*s;
	};
	vector<double> cuts = { cut(b0), cut(b1), cut(b2), cut(b3) };
	std::sort(cuts.begin(), cuts.end());
	this->cuts_range = std::make_pair(cuts[0], cuts[3]);
}

Model::Model( const std::tuple<std::tuple<double,double,double>, std::tuple<double,double,double>>& bbox )
: bbox(bbox), base_point(), direction(), topography(nullptr), horizontal(true)
{
	this->cuts_range = std::make_pair(g2(g0(bbox)), g2(g1(bbox)));
}

Model::~Model(){
	/* clears matches and sections from memory */
	if ( this->topography != nullptr ) {
		delete this->topography;
	}
	
	this->clear_matches();
	
	for ( auto it = this->sections.begin(); it != this->sections.end(); it++ ) {
		delete *it;
	}
}

void Model::clear_matches() {
	/* clears all the matches from memory */
	for ( Match * c: this->match ) {
		delete c;
	}
	this->match.clear();
}

void Model::make_matches() {
	this->clear_matches();
	
	auto& faults = this->global_faults;
	// Converts to global points.
	auto global_point = [&]( const point3& pt ) {
		point2 p(gx(pt), gy(pt));
		return this->inverse_point(p, gz(pt));
	};
	
	// Converts to global triangles.
	auto global_triangle = [global_point]( const triangle_pt& tr ) {
		return triangle_pt(global_point(g0(tr)), global_point(g1(tr)), global_point(g2(tr)));
	};
	
	// Converts to global coordinates and pushes to faults.
	auto add_to_faults = [&]( std::tuple<map<wstring, vector<triangle_pt>>, map<wstring, vector<size_t>>>& tup ) {
		auto& m = g0( tup );
		auto& e = g1( tup );
		for ( auto it = m.begin(); it != m.end(); it++ ) {
			// Transform to global coordinates.
			vector<triangle_pt>& trs = it->second;
			std::transform( trs.begin(), trs.end(), trs.begin(), global_triangle );
			// Push into the faults.
			auto& f = faults[it->first];
			auto& o = this->extended_faults[it->first];
			for ( size_t& c: e[it->first] ) {
				o.push_back( c + f.size() );
			}
			f.reserve(f.size() + it->second.size());
			f.insert(f.end(), it->second.begin(), it->second.end());
		}
	};
	
	for ( size_t i = 1; i < this->sections.size(); i++ ) {
		this->match.push_back(new Match(this->sections[i-1], this->sections[i]));
		// Match the polygons.
		this->match.back()->match_polygons();
		// Get the matching faults.
		auto m = this->match.back()->match_lines( this->feature_types );
		add_to_faults(m);
	}
	
	// Get the extended faults from the begining.
	if ( this->sections.size() and (this->sections[0]->cut - this->cuts_range.first) > tolerance ) {
		auto m = this->sections[0]->last_lines(true, this->cuts_range.first);
		add_to_faults(m);
	}
	
	// Get the extended faults from the end.
	if ( this->sections.size() and (this->cuts_range.second - this->sections.back()->cut) > tolerance ) {
		// Get the extended faults from the end.
		auto m = this->sections.back()->last_lines(false, this->cuts_range.second);
		add_to_faults(m);
	}
}

void Model::set_matches( const vector< std::tuple< std::tuple<wstring, wstring>, vector<std::pair<int, int>> > >& matching ) 
{
	this->clear_matches();
	map<std::pair<wstring, wstring>, int> match;
	
	for ( size_t i = 0; i < matching.size(); i++ ) {
		const wstring& name1 = g0(g0(matching[i]));
		const wstring& name2 = g1(g0(matching[i]));
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
		const vector<std::pair<int, int>>& m = g1(matching[match[std::make_pair(name1, name2)]]);
		this->match.push_back(new Match(this->sections[i-1], this->sections[i]));
		this->match.back()->set(m);
	}
}

std::pair<int, double> Model::closest_match( bool a, int a_idx, int pol_idx, const point2& pt ) const {
	// Reference to match.
	const Match& m = *(this->match[a_idx]);
	// Reference to section s.
	const Section& s = (a) ? *(this->sections[a_idx+1]) : *(this->sections[a_idx]);
	// Get match depending on the direction.
	const map<int, vector<int>>& mp = (a) ? m.a_to_b : m.b_to_a;
	auto it = mp.find(pol_idx);
	// This polygon does not have matches.
	if ( it == mp.end() ) {
		return std::make_pair(-1, std::numeric_limits<double>::infinity());
	}
	// In all the matches of this polygon, check the minimum distance.
	const vector<int>& op = it->second;
	double mindist = std::numeric_limits<double>::infinity();
	int minidx = -1;
	for ( size_t i = 0; i < op.size(); i++ ) {
		size_t pl = op[i];
		double dist = geometry::distance(s.polygons[pl], pt);
		if ( dist < mindist ) {
			mindist = dist;
			minidx = pl;
		}
	}
	// Return minimum distance and index of match.
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

double Model::signed_distance( const wstring& unit, const point3& pt ) const{
	
	auto all_except = [unit](const value_f& v) {
		const wstring& s = g1(v);
		return s != unit;
	};
	
	auto just = [unit](const value_f& v) {
		const wstring& s = g1(v);
		return s == unit;
	};

	std::tuple<wstring, double> inside = this->closest( pt, just );
	std::tuple<wstring, double> outside = this->closest( pt, all_except );
	if ( g0(inside) == L"NONE" ) {
		return g1(inside);
	}
	return g1(inside) - g1(outside);
}

double Model::signed_distance_bounded( const wstring& unit, const point3& pt ) const {
	double sdist = this->signed_distance( unit, pt );
	bool outside = false;
	double odist = 0.0;
	double idist = -std::numeric_limits<double>::infinity(); // In the future, fix distance below too.
	
	double x = gx(pt);
	double y = gy(pt);
	double z = gz(pt);
	
	double minx = g0(g0(this->bbox));
	double miny = g1(g0(this->bbox));
	double minz = g2(g0(this->bbox));
	
	double maxx = g0(g1(this->bbox));
	double maxy = g1(g1(this->bbox));
	double maxz = this->height( point2( x, y ) );
	
	double dists[6] = { minx - x, miny - y, minz - z, x - maxx, y - maxy, z - maxz };
	
	for ( size_t i = 0; i < 6; i++ ) {
		if ( dists[i] >= 0 ) {
			outside = true;
			odist += dists[i]*dists[i];
		} else {
			idist = std::max(idist, dists[i]);
		}
	}
	
	if ( outside ) {
		return std::max(sdist, std::sqrt(odist));
	}
	return std::max(sdist, idist);
}

double Model::signed_distance_unbounded( const wstring& unit, const point3& pt ) const {
	double sdist = this->signed_distance( unit, pt );
	
	double x = gx(pt);
	double y = gy(pt);
	double z = gz(pt);
	
	double maxz = this->height( point2( x, y ) );
	double distz = z - maxz;
	
	if ( distz >= 0 ) {
		return std::max(sdist, distz);
	}
	return std::max(sdist, distz);
}


vector<Model::Possible> Model::all_closest( size_t a_idx, const point2& pt_a, const point2& pt_b ) const {
	vector<Model::Possible> possible = this->get_candidates(a_idx, pt_a, pt_b, always_true);
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


std::pair<point2, double> Model::model_point( const point3& pt ) const {
	// If it's a horizontal cross section, just translate it.
	if ( this->horizontal ) {
		return std::make_pair(point2(gx(pt), gy(pt)), gz(pt));
	}
	point2 ppoint(gx(pt), gy(pt));
	const double& z = gz(pt);
	geometry::subtract_point(ppoint, this->base_point);
	
	double norm = std::sqrt(gx(ppoint)*gx(ppoint) + gy(ppoint)*gy(ppoint));

	if ( norm < tolerance ) {
		return std::make_pair(point2(0.0, z), 0.0);
	}
	
	geometry::divide_value(ppoint, norm);

	const double& d0 = gx(this->direction);
	const double& d1 = gy(this->direction);
	
	double c = d0*gy(ppoint) - d1*gx(ppoint);
	double d = d0*gx(ppoint) + d1*gy(ppoint);
	return std::make_pair(point2(d*norm, z), c*norm);
}

point3 Model::inverse_point( const point2& pt, double cut ) const {
	// If it's a horizontal cross section, just translate it.
	if ( this->horizontal ) {
		return point3( gx(pt), gy(pt), cut );
	}
	const double& d0 = gx(this->direction);
	const double& d1 = gy(this->direction);
	
	double dsq = d0*d0 + d1*d1;
	
	const double& p0 = gx(pt);
	const double& p1 = gy(pt);
	
        double v1 = (cut*d0 + p0*d1)/dsq;
        double v0 = ( std::fabs(d1) < tolerance ) ? p0/d0 : (d0*v1-cut)/d1;
        
	return point3(gx(this->base_point) + v0, gy(this->base_point) + v1, p1);
}

map<wstring,vector<line>> Model::intersect_planes(const vector<line_3d>& planes) const{
	return find_faults_multiple_planes_intersection(this->global_faults, planes);
}

map<wstring,vector<line>> Model::intersect_topography(const vector<vector<double>>& topography_array,
	double z_max, double z_min,double x_inf, double y_inf, double dx, double dy, int rows, int cols) const{

	return find_faults_topography_intersection(this->global_faults,topography_array, z_max, z_min,
		x_inf, y_inf, dx, dy, rows, cols);
}

map<wstring,vector<line>> Model::intersect_plane(const line_3d& plane) const{
	return find_faults_multiple_planes_intersection(this->global_faults, vector<line_3d>(1, plane));
}

std::pair<double, bool> Model::find_unit_limits(double xp, double yp,double z_max, double z_min, double eps) const{
	return  find_unit_limits_cpp(this, xp, yp, z_max, z_min, eps);
}

vector<std::tuple<wstring, double, double>> Model::possible_closest( const point3& pt ) const {
	if ( not this->sections.size() ) {
		return vector<std::tuple<wstring, double, double>>();
	}
	if ( this->match.size()+1 != this->sections.size() ) {
		throw GeomodelrException("You need to call make_matches before using this function.");
	}
	std::pair<point2, double> mp = this->model_point(pt);
	auto it = std::upper_bound(this->cuts.begin(), this->cuts.end(), mp.second);
	size_t a_idx = it - this->cuts.begin();
	// If it's behind the last or above the first, return the closest in the section.
	
	if ( a_idx <= 0 or a_idx >= this->sections.size() ) {
		const Section& s = ( a_idx <=0 ) ? *(this->sections.front()) : *(this->sections.back());
		std::pair<int, double> cls = s.closest(mp.first, always_true);
		wstring unit = s.units[cls.first];
		vector<std::tuple<wstring, double, double>> ret = { std::make_tuple(unit, cls.second, cls.second) };
		return ret;
	}
	
	a_idx--;
	
	vector<Model::Possible> ap = this->all_closest(a_idx, mp.first, mp.first);
	
	vector<std::tuple<wstring, double, double>> ret;
	
	for ( size_t i = 0; i < ap.size(); i++ ) {
		int a_match = ap[i].a_match;
		int b_match = ap[i].b_match;
		wstring unit = ( a_match != -1 ) ? this->sections[a_idx]->units[a_match] : this->sections[a_idx+1]->units[b_match];
		ret.push_back(std::make_tuple(unit, ap[i].a_dist, ap[i].b_dist));
	}
	
	return ret;

}

std::tuple<wstring, double> Model::closest_topo( const point3& pt ) const {
	if ( this->topography != nullptr ) {
		if ( this->height(point2(gx(pt), gy(pt))) < gz(pt) ) {
			return std::make_tuple(wstring(L"AIR"), std::numeric_limits<double>::infinity());
		}
	}
	return this->closest(pt);
}

std::tuple<wstring, double> Model::closest( const point3& pt ) const {
	return this->closest(pt, always_true);
}

double Model::height( const point2& pt ) const {
	if ( this->topography != nullptr ) {
		return this->topography->height(pt);
	}
	return std::numeric_limits<double>::infinity();
}

pylist ModelPython::possible_closest(const pyobject& pypt) const {
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& p2 = python::extract<double>(pypt[2]);
	point3 pt(p0, p1, p2);
	vector<std::tuple<wstring, double, double>> possible = ((Model *)this)->possible_closest(pt);
	pylist ret;
	for ( size_t i = 0; i < possible.size(); i++ ) {
		ret.append(python::make_tuple(g0(possible[i]), g1(possible[i]), g2(possible[i])));
	}
	return ret;
}


pytuple ModelPython::closest_topo( const pyobject& pypt ) const {
	double p0 = python::extract<double>(pypt[0]);
	double p1 = python::extract<double>(pypt[1]);
	double p2 = python::extract<double>(pypt[2]);
	std::tuple<wstring, double> ret = ((Model *)this)->closest_topo( point3( p0, p1, p2 ) );
	return python::make_tuple(g0(ret), g1(ret));
}

pytuple ModelPython::model_point( const pyobject& pypt ) const {
	double p0 = python::extract<double>(pypt[0]);
	double p1 = python::extract<double>(pypt[1]);
	double p2 = python::extract<double>(pypt[2]);
	
	point3 pt(p0, p1, p2);
	std::pair<point2, double> mp = ((Model *)this)->model_point( pt );
	
	double m0 = gx(mp.first);
	double m1 = gy(mp.first);
	return python::make_tuple(m0, m1, mp.second);
}

pytuple ModelPython::inverse_point( const pyobject& pypt ) const {
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& cut = python::extract<double>(pypt[2]);
	point2 pt(p0, p1);
	point3 mp = ((Model *)this)->inverse_point( pt, cut );
	const double& m0 = gx(mp);
	const double& m1 = gy(mp);
	const double& m2 = gz(mp);
	return python::make_tuple(m0, m1, m2);
}

pydict ModelPython::info() const {
	pydict ret;
	for ( const Section * s: this->sections ) {
		ret[((SectionPython *)s)->name] = ((SectionPython *)s)->info();
	}
	return ret;
}

pytuple ModelPython::closest(const pyobject& pypt) const {
	// python exposed functions to search for the closest formation.
	const double& p0 = python::extract<double>(pypt[0]);
	const double& p1 = python::extract<double>(pypt[1]);
	const double& p2 = python::extract<double>(pypt[2]);
	point3 pt(p0, p1, p2);
	std::tuple<wstring, double> res = ((Model *)this)->closest( pt );
	return python::make_tuple(g0(res), g1(res));
}

pydict ModelPython::intersect_planes(const pylist& planes) const{
	
	vector<line_3d> planes_cpp;
    for (int i=0; i<python::len(planes); i++){
        int N = python::len(python::extract<pylist>(planes[i]));
        line_3d plane;
        for (int j=0; j<N; j++){
            plane.push_back(point3(python::extract<double>(planes[i][j][0]),python::extract<double>(planes[i][j][1]),
                python::extract<double>(planes[i][j][2])) );
        }
        planes_cpp.push_back(plane);
	}
		
	// Now convert intersections to python and return.
	return (map_to_pydict(((Model *)this)->intersect_planes(planes_cpp)));
	
}

pydict ModelPython::intersect_topography(const pydict& topography_info) const{
	
	double x_inf = python::extract<double>(topography_info["point"][0]);
    double y_inf = python::extract<double>(topography_info["point"][1]);

    double dx = python::extract<double>(topography_info["sample"][0]);
    double dy = python::extract<double>(topography_info["sample"][1]);

    int rows = python::extract<int>(topography_info["dims"][1]);
    int cols = python::extract<int>(topography_info["dims"][0]);

    double z_max, z_min;
    vector<vector<double>> topography_array = topography_to_vector(python::extract<pylist>(topography_info["heights"]),
        rows, cols, z_max, z_min);
		
	// Now convert intersections to python and return.
	return map_to_pydict(((Model *)this)->intersect_topography(topography_array, z_max, z_min, x_inf, y_inf, dx, dy,
        rows, cols));
	
}

pydict ModelPython::intersect_plane(const pylist& plane) const{
	
    line_3d plane_cpp;
    for (int j=0; j<python::len(plane); j++){
        plane_cpp.push_back(point3(python::extract<double>(plane[j][0]),python::extract<double>(plane[j][1]),
            python::extract<double>(plane[j][2])));
    }
	// Now convert intersections to python and return.
	return (map_to_pydict(((Model *)this)->intersect_plane(plane_cpp)));
}

pytuple ModelPython::find_unit_limits(double xp, double yp, double z_max, double z_min, double eps) const{

	std::pair<double, bool> output = ((Model *)this)->find_unit_limits(xp, yp, z_max, z_min, eps);
	return python::make_tuple(output.first,output.second);
}


double ModelPython::height( const pyobject& pt ) const {
	double d0 = python::extract<double>(pt[0]);
	double d1 = python::extract<double>(pt[1]);
	return ((Model *)this)->height( point2( d0, d1 ) );
}

Topography::Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims ):
point(point), sample(sample), heights(dims[0]*dims[1])
{
	this->dims[0] = dims[0];
	this->dims[1] = dims[1];
}


Topography::Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims, const vector<double>& heights ):
Topography( point, sample, dims ) {
	int rows = dims[0];
	int cols = dims[1];
	for ( int i = 0; i < rows; i++ ) {
		for ( int j = 0; j < cols; j++ ) {
			this->heights[i*this->dims[1]+j] = heights[i*this->dims[1]+j];
		}
	}
}

double Topography::height(const point2& pt) const {
	point2 pos(gx(pt), gy(pt));
	geometry::subtract_point(pos, this->point);
	geometry::divide_point(pos, this->sample);
	int x = int(gx(pos));
	int y = int(gy(pos));

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



TopographyPython::TopographyPython( const pyobject& point, const pyobject& sample, const pyobject& dims, const pylist& heights ):
Topography(point2(python::extract<double>(point[0]), python::extract<double>(point[1])),
	   point2(python::extract<double>(sample[0]), python::extract<double>(sample[1])),
           {python::extract<int>(dims[0]), python::extract<int>(dims[1])})
{
	int rows = python::len(heights);
	if ( rows != this->dims[0] ) {
		throw GeomodelrException("topography rows does not correspond with dims.");
	}
	for ( int i = 0; i < rows; i++ ) {
		int cols = python::len(heights[i]);
		if ( cols != this->dims[1] ) {
			throw GeomodelrException("topography columns does not correspond with dims.");
		}
		for ( int j = 0; j < cols; j++ ) {
			this->heights[i*this->dims[1]+j] = python::extract<double>(heights[i][j]);
		}
	}
}

double TopographyPython::height( const pyobject& pypt ) const {
	return ((Topography *)this)->height( point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])) );
}

void ModelPython::fill_model( const pyobject& topography, const pylist& sections, const pydict& feature_types ) {
	if ( python::len(topography) > 0 ) {// Get all the information sent from python to build the topography.
		
		const pyobject& point = python::extract<pyobject>(topography["point"]);
		const pyobject& sample = python::extract<pyobject>(topography["sample"]);
		const pyobject& dims = python::extract<pyobject>(topography["dims"]);
		const pylist& heights = python::extract<pylist>(topography["heights"]);
		this->topography = new TopographyPython(point, sample, dims, heights);
	}
	size_t nsects = python::len(sections);
	for ( size_t i = 0; i < nsects; i++ ) {
		// Get all the information sent from python to create each cross section.
		if ( python::len( sections[i] ) != 9 ) {
			throw GeomodelrException("Every section needs 8 parameters.");
		}
		wstring name    = python::extract<wstring>(sections[i][0]);
		double cut             = python::extract<double>(sections[i][1]);
		const pyobject& bbox   = python::extract<pyobject>(sections[i][2]);
		const pylist& points   = python::extract<pylist>(sections[i][3]);
		const pylist& polygons = python::extract<pylist>(sections[i][4]);
		const pylist& units    = python::extract<pylist>(sections[i][5]); 
		const pylist& lines    = python::extract<pylist>(sections[i][6]);
		const pylist& lnames   = python::extract<pylist>(sections[i][7]);
		const pylist& anchored_lines = python::extract<pylist>(sections[i][8]);
		
		// Pass all the information to the cross section, including the bounding box of the given section.
		this->sections.push_back( new SectionPython( name, cut, bbox, points, polygons, units, lines, lnames, anchored_lines ) );
	}
	
	std::sort(this->sections.begin(), this->sections.end(), [](const Section* a, const Section* b){ return a->cut < b->cut; });
	for ( size_t i = 0; i < this->sections.size(); i++ ) {
		this->cuts.push_back(this->sections[i]->cut);
	}
	pylist keys = feature_types.keys();
	for ( int i = 0; i < python::len( keys ); i++ ) {
		this->feature_types[python::extract<wstring>(keys[i])] = python::extract<wstring>(feature_types[keys[i]]);
	}
}

ModelPython::ModelPython( const pyobject& bbox, const pyobject& base_point, 
		    	 const pyobject& direction, const pyobject& map, 
		    	 const pyobject& topography, const pylist& sections, 
			 const pydict& feature_types ): 
	Model(make_tuple(std::tuple<double, double, double>(python::extract<double>(bbox[0]),python::extract<double>(bbox[1]),python::extract<double>(bbox[2])), 
			 std::tuple<double, double, double>(python::extract<double>(bbox[3]),python::extract<double>(bbox[4]),python::extract<double>(bbox[5]))),
		point2( python::extract<double>(base_point[0]), python::extract<double>(base_point[1]) ),
		point2( python::extract<double>(direction[0]),  python::extract<double>(direction[1])  ))
{
	for ( int i = 0; i < python::len( sections ); i++ ) {
		double cut = python::extract<double>( sections[i][1] );
		pytuple sbbox = calculate_section_bbox( bbox, base_point, direction, cut );
		pylist pl = python::extract<pylist>( sections[i] );
		pl.insert( 2, sbbox );
	}
	this->fill_model( topography, sections, feature_types );
}


ModelPython::ModelPython( const pyobject& bbox, const pyobject& map, 
			  const pyobject& topography, const pylist& sections,
			  const pydict& feature_types ): 
	Model(make_tuple(std::tuple<double, double, double>(python::extract<double>(bbox[0]),python::extract<double>(bbox[1]),python::extract<double>(bbox[2])), 
			 std::tuple<double, double, double>(python::extract<double>(bbox[3]),python::extract<double>(bbox[4]),python::extract<double>(bbox[5]))))
{
	pytuple sbbox = python::make_tuple( bbox[0], bbox[1], bbox[3], bbox[4] );
	
	for ( int i = 0; i < python::len( sections ); i++ ) {
		pylist pl = python::extract<pylist>( sections[i] );
		pl.insert( 2, sbbox );
	}
	
	this->fill_model( topography, sections, feature_types );
}


void ModelPython::make_matches() {
	// map<wstring, vector<triangle_pt>> faults = ((Model *)this)->make_matches();
	((Model *)this)->make_matches();
	
}

pydict ModelPython::filter_lines( bool ext, const wstring& ft ) const {
	auto python_point = [&]( const point3& pt ) {
		return python::make_tuple(gx(pt), gy(pt), gz(pt));
	};
	
	auto python_triangle = [python_point]( const triangle_pt& tr ) {
		return python::make_tuple(python_point(g0(tr)), python_point(g1(tr)), python_point(g2(tr)));
	};
	
	// Now convert faults to python and return.
	pydict ret;
	auto add_line_triangles = [&]( map<wstring, vector<triangle_pt>>::const_iterator it ) {
		pylist tris;
		for ( const triangle_pt& t: it->second ) {
			pytuple tr = python_triangle( t );
			tris.append( tr );
		}
		ret[it->first] = tris;
	};
	
	auto add_line_triangles_not_ext = [&]( map<wstring, vector<triangle_pt>>::const_iterator it ) {
		pylist tris;
		// Find the extended lines of the given plane.
		const auto ef = this->extended_faults.find(it->first);
		for ( size_t i = 0; i < it->second.size(); i++ ) {
			if ( ef != this->extended_faults.end() ) {
				const vector<size_t>& vct = ef->second;
				if ( std::find(vct.begin(), vct.end(), i) != vct.end() ) {
					continue;
				}
			}
			const triangle_pt& t = it->second[i];
			pytuple tr = python_triangle( t );
			tris.append( tr );
		}
		ret[it->first] = tris;
	};
	
	for ( auto it = this->global_faults.begin(); it != this->global_faults.end(); it++ ) {
		if ( ft != L"ALL" ) 
		{
			wstring lft = (this->feature_types.find(it->first))->second;
			if ( lft != ft ) {
				continue;
			}
		}
		if ( not ext ) {
			add_line_triangles_not_ext( it );
		} else {
			add_line_triangles( it );
		}
	}
	
	return ret;
}

// Direct methods to get different kinds of extruded lines.
pydict ModelPython::get_lines() const { return this->filter_lines( true, L"ALL" ); }
pydict ModelPython::get_not_extended_lines() const { return this->filter_lines( false, L"ALL" ); }

pydict ModelPython::get_fracts( ) const { return this->filter_lines( true, L"FRACT" ); }
pydict ModelPython::get_not_extended_fracts() const { return this->filter_lines( false, L"FRACT" ); }

pydict ModelPython::get_faults( ) const { return this->filter_lines( true, L"FAULT" ); }
pydict ModelPython::get_not_extended_faults() const { return this->filter_lines( false, L"FAULT" ); }

pydict ModelPython::get_veins( ) const { return this->filter_lines( false, L"VEIN" ); }
pydict ModelPython::get_not_extended_veins() const { return this->filter_lines( false, L"VEIN" ); }

void ModelPython::set_matches( const pylist& matching ) {
	this->clear_matches();
	size_t nmatch = python::len(matching);
	vector< std::tuple< std::tuple<wstring, wstring>, vector< std::pair<int, int>>>> cppmatching;
	for ( size_t i = 0; i < nmatch; i++ ) {
		const wstring& name1 = python::extract<wstring>(matching[i][0][0]);
		const wstring& name2 = python::extract<wstring>(matching[i][0][1]);
		vector<std::pair<int, int>> vmatch;
		size_t mmatch = python::len(matching[i][1]);
		for ( size_t j = 0; j < mmatch; j++ ) {
			int a = python::extract<int>(matching[i][1][j][0]);
			int b = python::extract<int>(matching[i][1][j][1]);
			
			vmatch.push_back(std::make_pair(a, b));
		}
		cppmatching.push_back(std::make_tuple( std::make_tuple(name1, name2), vmatch ));
	}
	((Model *)this)->set_matches( cppmatching );
}

pylist ModelPython::get_matches( ) const {
	pylist ret;
	for ( size_t i = 0; i < this->match.size(); i++ ) {
		const wstring& name1 = this->sections[i]->name;
		const wstring& name2 = this->sections[i+1]->name;
		ret.append(python::make_tuple(python::make_tuple(name1, name2), ((MatchPython *)(this->match[i]))->get()));
	}
	return ret;
}

pylist ModelPython::pybbox( ) const {
	pylist p;
	p.append( g0( g0( this->bbox ) ) );
	p.append( g1( g0( this->bbox ) ) );
	p.append( g2( g0( this->bbox ) ) );
	p.append( g0( g1( this->bbox ) ) );
	p.append( g1( g1( this->bbox ) ) );
	p.append( g2( g1( this->bbox ) ) );
	return p;
}

double ModelPython::signed_distance( const wstring& unit, const pyobject& pt ) const {
	return ((Model *)this)->signed_distance(unit, point3(python::extract<double>(pt[0]), python::extract<double>(pt[1]), python::extract<double>(pt[2])));
}
	
double ModelPython::signed_distance_bounded( const wstring& unit, const pyobject& pt ) const {
	return ((Model *)this)->signed_distance_bounded(unit, point3(python::extract<double>(pt[0]), python::extract<double>(pt[1]), python::extract<double>(pt[2])));
}

double ModelPython::signed_distance_unbounded( const wstring& unit, const pyobject& pt ) const {
	return ((Model *)this)->signed_distance_unbounded(unit, point3(python::extract<double>(pt[0]), python::extract<double>(pt[1]), python::extract<double>(pt[2])));
}

pytuple calculate_section_bbox( const pyobject& bbox, const pyobject& point, const pyobject& direction, double cut ) {
	
	// Get the point of the cross section.
	double x = python::extract<double>( direction[0] ), y = python::extract<double>( direction[1] );
	point2 d( x, y );
	x = python::extract<double>( point[0] );
	y = python::extract<double>( point[1] );
	point2 bp( x, y );
	
	double minx = python::extract<double>( bbox[0] );
	double miny = python::extract<double>( bbox[1] );
	double minz = python::extract<double>( bbox[2] );
	
	double maxx = python::extract<double>( bbox[3] );
	double maxy = python::extract<double>( bbox[4] );
	double maxz = python::extract<double>( bbox[5] );
	
	point2 p( gy( d ), -gx(d) );
	geometry::multiply_value( p, cut );
	geometry::add_point( p, bp );
	
	// Get all the cuts ( values of x in the equation p + x*v = bbox ).
	vector<double> xs;
	// Add the cuts of X.
	if ( std::fabs( gx( d ) ) > tolerance ) {
		xs.push_back( ( minx - gx( p ) ) / gx( d ) );
		xs.push_back( ( maxx - gx( p ) ) / gx( d ) );
	}
	
	// Add the cuts of Y.
	if ( std::fabs( gy( d ) ) > tolerance ) {
		xs.push_back( ( miny - gy( p ) ) / gy( d ) );
		xs.push_back( ( maxy - gy( p ) ) / gy( d ) );
	}
	
	// Sort them.
	std::sort( xs.begin(), xs.end() );
	
	if ( xs.size() == 4 ) { // In case the xs cut four times, return the two in the middle.
		return python::make_tuple( xs[1], minz, xs[2], maxz );
	} else if ( xs.size() == 2 ) { // In case the xs cut two times, return them both.
		return python::make_tuple( xs[0], minz, xs[1], maxz );
	}
	throw GeomodelrException("The vector direction is invalid.");
}


