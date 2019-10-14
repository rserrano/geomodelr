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

Model::Model( const std::tuple<std::tuple<double,double,double>, std::tuple<double,double,double>>& bbox, 
			  const std::tuple<std::tuple<double,double,double>, std::tuple<double,double,double>>& abbox,
		  const point2& base_point, const point2& direction)
: bbox(bbox), abbox(abbox), base_point(base_point), direction(direction), geomap(nullptr), topography(nullptr), horizontal(false), faults_disabled(false)
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

	point3 min = point3( g0(g0(bbox)) , g1(g0(bbox)) , g2(g0(bbox)));
	point3 max = point3( g0(g1(bbox)) , g1(g1(bbox)) , g2(g1(bbox)));
	geometry::subtract_point(max,min);
	this->bbox_diag = std::sqrt(geometry::dot_product(max,max));

	// this->set_signed_distance(L"classic");
	// this->set_soil_depth(-std::numeric_limits<double>::infinity());
}

Model::Model( const std::tuple<std::tuple<double,double,double>, std::tuple<double,double,double>>& bbox,
			  const std::tuple<std::tuple<double,double,double>, std::tuple<double,double,double>>& abbox)
: bbox(bbox), abbox(abbox), base_point(), direction(), geomap(nullptr), topography(nullptr), horizontal(true), faults_disabled(false)
{
	this->cuts_range = std::make_pair(g2(g0(bbox)), g2(g1(bbox)));

	point3 min = point3( g0(g0(bbox)) , g1(g0(bbox)) , g2(g0(bbox)));
	point3 max = point3( g0(g1(bbox)) , g1(g1(bbox)) , g2(g1(bbox)));
	geometry::subtract_point(max,min);
	this->bbox_diag = std::sqrt(geometry::dot_product(max,max));
	this->set_signed_distance(L"classic");
}

Model::~Model() {
	/* clears matches and sections from memory */
	if ( this->topography != nullptr ) {
		delete this->topography;
	}
	/* clears geological map from memory */
	if ( this->geomap != nullptr ) {
		delete this->geomap;
	}
	/* deletes every match individually */
	this->clear_matches();
	/* deletes every cross section individually */
	for ( Section * s: this->sections ) {
		delete s;
	}
}

void Model::clear_matches() {
	/* clears all the matches from memory */
	for ( Match * c: this->match ) {
		delete c;
	}
	this->match.clear();
};

void Model::set_params( const map<wstring, wstring>& params ) {
	this->params = params;

	for ( Section * s: this->sections ) {
		s->set_params(&(this->params));
	}
	
	for ( Match * m: this->match ) {
		m->set_params(&(this->params));
	}
	if ( this->geomap != nullptr ) {
		this->geomap->set_params(&(this->params));
	}
	auto it = this->params.find( L"faults" );
	if ( it != this->params.end() ) {
		if ( it->second == L"disabled" ) {
			this->faults_disabled = true;
		} else {
			this->faults_disabled = false;
		}
	} else {
		this->faults_disabled = false;
	}
	it = this->params.find( L"map" );
	if ( it != this->params.end() ) {
		if ( it->second == L"soils" ) {
			this->check_soils = true;
		} else {
			this->check_soils = false;
		}
	} else {
		this->check_soils = false;
	}
}

vector< std::pair<wstring, double> > Model::get_closests(const Section * s1, const Section * s2, 
	const point2& pt1, const point2& pt2, double d1, double d2, wstring unit, bool check_depth, bool two_sections,
	const point3& real_pt) const{

	auto aux_pair = std::make_pair(L"NONE", std::numeric_limits<double>::infinity() );
	vector< std::pair<wstring, double> > output { aux_pair, aux_pair, aux_pair };
	double theta = 1.0 - d1/(d1+d2);

	for (wstring uni: this->units){

		if (check_depth && (this->soils.find(uni) != this->soils.end())){
			continue;
		}

		double dis_s1 = s1->unit_distance(uni, pt1);
		double dis_s2;
		if (two_sections){
			dis_s2 = s2->unit_distance(uni, pt2);
		} else{
			dis_s2 = 0.0;
			theta = 1.0;
		}

		double distance;
		double r = 0.999;
		if (dis_s1 == std::numeric_limits<double>::infinity() and dis_s2 == std::numeric_limits<double>::infinity()){
			distance = this->bbox_diag;
		} else if(dis_s1 == std::numeric_limits<double>::infinity()){
			if (theta<=r and two_sections){
				distance = dis_s2 + (d1+d2)*(theta);
			} else{
				distance = this->bbox_diag;
			}
		} else if(dis_s2 == std::numeric_limits<double>::infinity()){
			if (theta>=(1.0-r)){
				distance = dis_s1 + (d1+d2)*(1.0 - theta);
			} else{
				distance = this->bbox_diag;
			}			
		} else{
			distance = theta*dis_s1 + (1-theta)*dis_s2;
		}

		if (distance < output[1].second){
			output[2] = output[1];
			output[1] = std::make_pair(uni, distance);
		} else if (distance < output[2].second){
			output[2] = std::make_pair(uni, distance);
		}
		if (uni == unit){
			output[0] = std::make_pair(unit, distance);
		}
	}

	if (output[0].first != unit){
		output[0] = std::make_pair(unit, this->bbox_diag);
	}


	//Check extrusion soils
	// std::wcerr << std::boolalpha << this->check_soils << std::endl;
	if (this->check_soils){
		point2 real_pt2 = point2(gx(real_pt), gy(real_pt));
		for (auto it = this->soil_depths.begin(); it != this->soil_depths.end(); ++it){

			// std::wcerr << it->first << std::endl;
			double dis = this->geomap->unit_distance(it->first, real_pt2);
			if (dis > 0.0){
				continue;
			}
			double h = this->height( real_pt2 );
			h -= (it->second + gz(real_pt));
			if (h>0){
				continue;
			}
			if (output[0].first == it->first){
				output[0].second = 0.0;
			}
			if (output[2].first == it->first){
				output[2] = output[1];
				output[1] = std::make_pair(it->first, 0.0);
			} else{
				// output[1].second = 0.0;
				output[1] = std::make_pair(it->first, 0.0);
			}
		}
	}

	return output;
}


vector<std::pair<wstring, double>> Model::closest_projection_vector( const point3& real_pt, wstring unit) const{
	
	std::pair<point2, double> mp = ((Model *)this)->model_point( real_pt );
	const point2& ft = mp.first;
	double sd = mp.second;

	if ( not this->sections.size() ) {
		auto aux_pair = std::make_pair(wstring(L"NONE"), std::numeric_limits<double>::infinity());
		return vector<std::pair<wstring, double>> {aux_pair, aux_pair, aux_pair};
	}

	auto it = std::upper_bound(this->cuts.begin(), this->cuts.end(), sd);	
	size_t a_idx = it - this->cuts.begin();
	double height = this->height( point2(gx(real_pt), gy(real_pt)));
	bool check_depth = (height - gz(real_pt)) > this->soil_depth;

	// Distance for single section
	auto closest_single = [&](const Section * s, double cut) {
		point2 pt = this->project_point(s, ft, sd, height);
		// std::pair<point3, double> topo_point = this->topography->intersection(real_pt, this->projection);
		std::pair<point3, double> topo_point = this->topography->intersects(real_pt, this->projection);

		double dx = gx(pt) - gx(ft);
		double dy = gy(pt) - gy(ft);
		double dz = cut - sd;
		double dist_pt = std::sqrt(dx*dx + dy*dy + dz*dz);

		if (topo_point.second == std::numeric_limits<double>::infinity()){
			return this->get_closests(s, s, pt, pt, dist_pt, 1.0, unit, check_depth, false, real_pt);
		}else{
			double dx = gx(pt) - gx(ft);
			double dy = gy(pt) - gy(ft);
			double dz = cut - sd;
			double dist_pt = std::sqrt(dx*dx + dy*dy + dz*dz);
			return this->get_closests(s, this->geomap, pt, point2( gx(topo_point.first), gy(topo_point.first)),
				dist_pt, topo_point.second, unit, check_depth, true, real_pt);
		}
	};

	if ( a_idx <= 0 ) {
		return closest_single(this->sections.front(), this->cuts.front());
	}
	
	if ( a_idx  >= this->sections.size() ) {
		return closest_single(this->sections.back(), this->cuts.back());
	}

	a_idx--;
	// Projected points over sections A and B.
	point2 pa = this->project_point(this->sections[a_idx], ft, sd, height);
	point2 pb = this->project_point(this->sections[a_idx+1], ft, sd, height);

	// Distance to sections A and B.
	double dx = gx(pa) - gx(ft);
	double dy = gy(pa) - gy(ft);
	double dz = this->cuts[a_idx] - sd;
	double dist_a = std::sqrt(dx*dx + dy*dy + dz*dz);

	dx = gx(pb) - gx(ft);
	dy = gy(pb) - gy(ft);
	dz = this->cuts[a_idx + 1] - sd;
	double dist_b = std::sqrt(dx*dx + dy*dy + dz*dz);
	
	// Topography point
	// std::pair<point3, double> topo_point = this->topography->intersection(real_pt, this->projection);
	std::pair<point3, double> topo_point = this->topography->intersects(real_pt, this->projection);

	// Topography and section B
	if (( dist_a > dist_b) and  (dist_a > topo_point.second)){
		return this->get_closests(this->geomap, this->sections[a_idx+1], point2( gx(topo_point.first), gy(topo_point.first)),
			pb, topo_point.second, dist_b, unit, check_depth, true, real_pt);
	// Topography and section A
	}else if (( dist_b > dist_a) and  (dist_b > topo_point.second)){
		return this->get_closests(this->geomap, this->sections[a_idx], point2( gx(topo_point.first), gy(topo_point.first)),
			pa, topo_point.second, dist_a, unit, check_depth, true, real_pt);
	// Section A and B.
	}else{
		return this->get_closests(this->sections[a_idx], this->sections[a_idx+1], pa, pb, dist_a, dist_b, unit, check_depth, true,
			real_pt);
	}
}

std::pair<wstring, double> Model::closest_projection( const point3& real_pt, wstring unit) const{
	return this->closest_projection_vector(real_pt, unit)[1];
}

void Model::set_projection(const point3& pt){

	// Get projection in model coordinates
	double x = gx(pt);
	double y = gy(pt);
	double z = gz(pt);

	double norm = std::sqrt(x*x + y*y + z*z);
	double val;
	if (horizontal){
		val = z;
	} else{
		val = gx(this->direction)*y - gy(this->direction)*x;
	}
	if (std::abs(val/norm) < epsilon){
		throw GeomodelrException("Th projection vector is parallel to the cross sections.");
	}

	this->projection = point3(x/norm, y/norm, z/norm);
}

void Model::set_projection_aligned(const point3& pt){
	double x = gx(pt);
	double y = gy(pt);
	double z = gz(pt);
	double norm = std::sqrt(x*x + y*y + z*z);
	if (std::abs(z/norm) < epsilon){
		throw GeomodelrException("Th projection vector is parallel to the cross sections.");
	}

	point3 aux = this->inverse_point( point2(x, y), z );
	x = gx(aux);
	y = gy(aux);
	z = gz(aux);

	if (not(horizontal)){
		x -= gx(this->base_point);
		y -= gy(this->base_point);
	}
	norm = std::sqrt(x*x + y*y + z*z);
	this->projection = point3(x/norm, y/norm, z/norm);
}

point3 Model::get_projection() const{
	return this->projection;
}

point2 Model::project_point(const Section * sec, const point2& pt, double z, double height) const{
	
	double cut = sec->cut;
	// const point3& projection = this->projection_aligned;

	double xp, yp, zp, delta;
	if (this->horizontal){
		xp = gx(this->projection);
		yp = gy(this->projection);
		zp = gz(this->projection);
		// delta = height - z;
		delta = std::numeric_limits<double>::infinity();
	} else{
		xp = gx(this->direction)*gx(this->projection) + gy(this->direction)*gy(this->projection);
		yp = gz(this->projection);
		zp = gx(this->direction)*gy(this->projection) - gy(this->direction)*gx(this->projection);
		delta = height - gy(pt);
	}

	double t = (cut - z)/zp;
	double x = gx(pt) + t*xp;

	if (this->soils.empty() or delta > this->soil_depth){
		return point2(x , gy(pt) + t*yp );
	} else{
		point3 aux_pt = this->inverse_point(point2(x, 0), cut);
		return point2(x , this->height( point2(gx(aux_pt), gy(aux_pt)) ) - delta );
	}
	
}

// Model--------------------------- ModelPthon

pydict ModelPython::get_params() const {
	pydict ret;
	for ( auto p: this->params ) {
		ret[p.first] = p.second;
	}
	return ret;
}

void ModelPython::set_params( const pydict& params ) {
	// pylist keys = params.keys();
	// map<wstring, wstring> out;
	// for ( int i = 0; i < python::len(keys); i++ ) {
	// 	out[python::extract<wstring>(keys[i])] = python::extract<wstring>(params[keys[i]]);
	// }

	map<wstring, wstring> out;
	if (params.has_key(wstring(L"faults"))){
		out[L"faults"] = python::extract<wstring>(params[wstring(L"faults")]);
	}
	if (params.has_key(wstring(L"map"))){
		out[L"map"] = python::extract<wstring>(params[wstring(L"map")]);
	}

	Model::set_params( out );

	if (params.has_key(wstring(L"projection"))){
		double x = python::extract<double>( params[wstring(L"projection")][0] );
		double y = python::extract<double>( params[wstring(L"projection")][1] );
		double z = python::extract<double>( params[wstring(L"projection")][2] );
		((Model* )this)->set_projection( point3(x, y, z) );
	} else{
		((Model* )this)->set_projection_aligned( point3(0, 0, 1) );
	}

	if (params.has_key(wstring(L"soil_depth"))){
		this->soil_depth = python::extract<double>( params[wstring(L"soil_depth")] );
	} else{
		this->soil_depth = -std::numeric_limits<double>::infinity();
	}

	if (params.has_key(wstring(L"mode"))){
		this->set_signed_distance( python::extract<wstring>( params[wstring(L"mode")] ) );
	} else{
		this->set_signed_distance( L"classic");
	}
}

pydict ModelPython::get_soil_depths() const {
	pydict ret;
	for ( auto p: this->soil_depths ) {
		ret[p.first] = p.second;
	}

	for (auto p: this->soils){
		ret[p] = pyobject();
	}
	
	return ret;
}

pylist ModelPython::get_soils() const {
	pylist ret;
	for ( auto p: this->soils ) {
		ret.append(p);
	}
	return ret;
}

void ModelPython::set_soil_depths( const pydict& depths ) {
	pylist keys = depths.keys();
	map<wstring, double> out;
	std::set<wstring> out_soils;

	for ( int i = 0; i < python::len( keys ); i++ ) {
		python::extract<double> value(depths[keys[i]]);
		if (value.check()){
			out[python::extract<wstring>(keys[i])] = python::extract<double>(depths[keys[i]]);	
		} else{
			out_soils.insert(python::extract<wstring>(keys[i]));
		}		
	}

	this->soil_depths = out;
	this->soils = out_soils;

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
		// Get the matching faults.
		auto m = this->match.back()->match_lines( this->feature_types );
		// Match the polygons.
		this->match.back()->match_polygons();
		add_to_faults(m);
		this->match.back()->set_params(&(this->params));
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
		double dist = s.poly_trees[pl]->distance_point(pt);
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

double Model::geomodelr_distance( const wstring& unit, const point3& pt ) const{
	
	auto just = [unit](const value_f& v) {
		const wstring& s = g1(v);
		return s == unit;
	};

	std::tuple<wstring, double> inside = this->closest_predicates( pt, just );
	return g1(inside);
}
vector<point2>  Model::get_polygon(const wstring sec, int poly_idx) const{

	int pos = 0;
	int n = this->sections.size();
	for (int k=0; k<n;k++){
		if (this->sections[k]->name==sec){
			pos = k;
			break;
		}
	}

	ring outer = (this->sections[pos]->poly_trees[poly_idx])->boost_poly.outer();
	size_t nnodes = outer.size();

	vector<point2> pts;
	for ( size_t k = 0; k < nnodes; k++ ) {
		pts.push_back(outer[k]);
	}
	return pts;

}

vector<point2> Model::get_fault(const wstring sec, int fault_idx) const{

	int pos = 0;
	int n = this->sections.size();
	for (int k=0; k<n;k++){
		if (this->sections[k]->name==sec){
			pos = k;
			break;
		}
	}

	line outer = (this->sections[pos]->lines[fault_idx]);
	size_t nnodes = outer.size();

	vector<point2> pts;
	for ( size_t k = 0; k < nnodes; k++ ) {
		pts.push_back(outer[k]);
	}
	return pts;

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
	
	std::tuple<wstring, double> inside = this->closest_predicates( pt, just );
	std::tuple<wstring, double> outside = this->closest_predicates( pt, all_except );

	// std::wcerr << g0(inside) << L"\t" << g1(inside) << std::endl;
	// std::wcerr << g0(outside) << L"\t" << g1(outside) << std::endl;

	double d;
	if ( g0(inside) == L"NONE" ) {
		d = g1(inside);
	} else if ( g0(outside) == L"NONE" ) {
		d = -g1(outside);
	} else {
		d = g1(inside) - g1(outside);
	}

	// Check if d is +-infinity.
	if (std::isinf(d)){
		if (d > 0){
			d = bbox_diag;
		} else{
			d = -bbox_diag;
		}
	}

	if ( this->check_soils ) {
		std::tuple<wstring, double> sl = this->soil( pt );
		// std::wcerr << "val: " << g1(sl) << L"\n";
		if ( g0(sl) == unit ) {
			return std::min( d, g1(sl) );
		} else {
			return std::max( d, -g1(sl) );
		}
	}
	return d;
}

double Model::signed_distance_projection( const wstring& unit, const point3& pt ) const {

	auto it = this->soil_depths.find(unit);
	if (this->check_soils and ( it != this->soil_depths.end())){
		point2 pt2 = point2(gx(pt), gy(pt));
		double dis = this->geomap->unit_distance(unit, pt2);
		double h = this->height( pt2 );
		h -= (it->second + gz(pt));
		if (dis <= 0.0){
			return h;
		}
	}

	auto res = this->closest_projection_vector( pt, unit);
	// for (auto a: res){
	// 	std::wcerr << a.first << L"\t" << a.second << std::endl;
	// }

	if (res[0].first != res[1].first){
		return (res[0].second - res[1].second);
	}else{
		return (res[0].second - res[2].second);
	}
}

double Model::signed_distance_bounded( const wstring& unit, const point3& pt ) const {
	// double sdist = this->signed_distance( unit, pt );
	double sdist = this->signed_distance_model(unit, pt);
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

double Model::signed_distance_bounded_restricted( const wstring& unit, const std::shared_ptr<Limiter> limit, const point3& pt ) const {
	// double sdist = this->signed_distance( unit, pt );
	// double sdist = this->signed_distance_projection( unit, pt );
	double sdist = this->signed_distance_model(unit, pt);
	
  return limit->limit_signed_distance(pt, sdist);
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


void Model::set_soil_depth(double depth){
	this->soil_depth = depth;
}

std::tuple<wstring, double> Model::soil( const point3& pt ) const {
	if ( this->geomap == nullptr ) {
		throw GeomodelrException("To check a soil you need a geological map.");
	}
	point2 xy( gx(pt), gy(pt) );
	double h = this->height( xy );
	std::pair<int, double> cls = this->geomap->closest(xy);
	if ( cls.second > 0.0 ) {
		return std::make_tuple( L"NONE", std::numeric_limits<double>::infinity() );
	}
	std::wstring unit = this->geomap->units[cls.first];
	auto it = this->soil_depths.find(unit);
	if ( it == this->soil_depths.end() ) {
		return std::make_tuple( L"NONE", std::numeric_limits<double>::infinity() );
	}
	h -= it->second;
	return std::make_tuple( unit, h-gz(pt) );
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

bbox3 Model::get_bbox() const{
	return this->bbox;
}

bbox3 Model::get_abbox() const{
	return this->abbox;
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
	if ( this->check_soils ) {
		std::tuple<wstring, double> sl = this->soil( pt );
		if ( g1( sl ) <= 0.0 ) {
			return std::make_tuple( g0(sl), 0.0 );
		}
	}
	return this->closest_predicates(pt, always_true);
}

void Model::set_signed_distance(const wstring& sdf_mode){

	if (sdf_mode == L"classic"){
		this->signed_distance_model = std::bind(&Model::signed_distance, this, std::placeholders::_1,
			std::placeholders::_2);

		this->closest_model = std::bind(&Model::closest, this, std::placeholders::_1);

	} else if (sdf_mode == L"vectorial"){
		this->signed_distance_model = std::bind(&Model::signed_distance_projection, this, std::placeholders::_1,
			std::placeholders::_2);

		this->closest_model = std::bind(&Model::closest_projection, this, std::placeholders::_1, L"NONE");
		
	} else{
		throw GeomodelrException("Signed distance mode have to be 'classic' or 'vectorial'.");
	}

	this->mode = sdf_mode;
}


double Model::height( const point2& pt ) const {
	if ( this->topography != nullptr ) {
		return this->topography->height_new(pt);
	}
	return std::numeric_limits<double>::infinity();
}

void ModelPython::set_signed_distance(const wstring& sdf_mode){
	Model::set_signed_distance(sdf_mode);
}

void ModelPython::set_projection(const pyobject& pypt){
	double p0 = python::extract<double>(pypt[0]);
	double p1 = python::extract<double>(pypt[1]);
	double p2 = python::extract<double>(pypt[2]);
	((Model *)this)->set_projection(point3(p0, p1, p2));
}

void ModelPython::set_projection_aligned(const pyobject& pypt){
	double p0 = python::extract<double>(pypt[0]);
	double p1 = python::extract<double>(pypt[1]);
	double p2 = python::extract<double>(pypt[2]);
	((Model *)this)->set_projection_aligned(point3(p0, p1, p2));
}

pytuple ModelPython::closest_projection(const pyobject& pypt) const {
	// python exposed functions to search for the closest formation.
	double p0 = python::extract<double>(pypt[0]);
	double p1 = python::extract<double>(pypt[1]);
	double p2 = python::extract<double>(pypt[2]);
	auto res = ((Model *)this)->closest_projection( point3(p0, p1, p2), L"NONE");
	return python::make_tuple(res.first, res.second);
}

void ModelPython::set_soil_depth(double depth){
	this->soil_depth = depth;
}

double ModelPython::get_soil_depth() const{
	return this->soil_depth;
}

pytuple ModelPython::get_projection() const{
	const point3& pt = ((Model* )this)->get_projection();
	return python::make_tuple(gx(pt), gy(pt), gz(pt));
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
	double p0 = python::extract<double>(pypt[0]);
	double p1 = python::extract<double>(pypt[1]);
	double cut = python::extract<double>(pypt[2]);
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
	double p0 = python::extract<double>(pypt[0]);
	double p1 = python::extract<double>(pypt[1]);
	double p2 = python::extract<double>(pypt[2]);
	point3 pt(p0, p1, p2);
	std::tuple<wstring, double> res = ((Model *)this)->closest_model( pt );
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

double ModelPython::height( const pyobject& pt ) const {
	double d0 = python::extract<double>(pt[0]);
	double d1 = python::extract<double>(pt[1]);
	return ((Model *)this)->height( point2( d0, d1 ) );
}

void ModelPython::fill_model( const pylist& geomap, const pyobject& topography, const pylist& sections,
	const pydict& feature_types, const pydict& params, const pydict& pyunits) {

	this->set_params( params );	
	if ( python::len(topography) > 0 ) {// Get all the information sent from python to build the topography.
		pyobject point = python::extract<pyobject>(topography["point"]);
		pyobject sample = python::extract<pyobject>(topography["sample"]);
		pyobject dims = python::extract<pyobject>(topography["dims"]);
		pylist heights = python::extract<pylist>(topography["heights"]);
		this->topography = new TopographyPython(point, sample, dims, heights);
	}
	// Soils info
	pylist units = pyunits.keys();
	size_t num = python::len(units);
	this->units.resize(num);

	map<wstring, double> soils_val;
	std::set<wstring> soils_none;
	for (size_t k=0; k<num; k++){
		wstring uni = python::extract<wstring>(units[k]);
		this->units[k] = uni;
		bool is_soil = python::extract<bool>(pyunits[uni][wstring(L"soil")]);
		if (is_soil){
			python::extract<double> value( pyunits[uni][wstring(L"depth")] );
			if (value.check()){
				soils_val[uni] = python::extract<double>( pyunits[uni][wstring(L"depth")] );
			} else{
				soils_none.insert(uni);
			}
		}
	}
	this->soil_depths = soils_val;
	this->soils = soils_none;
	
	size_t nsects = python::len(sections);	
	int lgeomap = python::len( geomap );
	if ( lgeomap > 1 ) {
		// Get all the information sent from python to create each cross section.
		if ( lgeomap != 7 ) {
			throw GeomodelrException("Every map needs 6 parameters.");
		}
		pyobject sbbox   = python::extract<pyobject>(geomap[0]);
		pylist points   = python::extract<pylist>(geomap[1]);
		pylist polygons = python::extract<pylist>(geomap[2]);
		pylist units	= python::extract<pylist>(geomap[3]); 
		pylist lines	= python::extract<pylist>(geomap[4]);
		pylist lnames   = python::extract<pylist>(geomap[5]);
		pylist anchored_lines = python::extract<pylist>(geomap[6]);
		
		// Pass all the information to the cross section, including the bounding box of the given section.
		this->geomap = (Section *) ( new GeologicalMapPython( sbbox, points, polygons, units, lines, lnames, anchored_lines ) );
		this->geomap->set_params( &(this->params) );
	}
	
	for ( size_t i = 0; i < nsects; i++ ) {
		// Get all the information sent from python to create each cross section.
		if ( python::len( sections[i] ) != 9 ) {
			throw GeomodelrException("Every section needs 8 parameters.");
		}
		
		wstring name	= python::extract<wstring>(sections[i][0]);
		double cut			 = python::extract<double>(sections[i][1]);
		pyobject bbox   = python::extract<pyobject>(sections[i][2]);
		pylist points   = python::extract<pylist>(sections[i][3]);
		pylist polygons = python::extract<pylist>(sections[i][4]);
		pylist units	= python::extract<pylist>(sections[i][5]); 
		pylist lines	= python::extract<pylist>(sections[i][6]);
		pylist lnames   = python::extract<pylist>(sections[i][7]);
		pylist anchored_lines = python::extract<pylist>(sections[i][8]);
		
		// Pass all the information to the cross section, including the bounding box of the given section.
		this->sections.push_back( new SectionPython( name, cut, bbox, points, polygons, units, lines, lnames, anchored_lines ) );
		this->sections.back()->set_params( &(this->params) );
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

// Model python for VERTICAL cross sections
ModelPython::ModelPython( const pyobject& bbox, const pyobject& abbox, const pyobject& base_point, 
				 const pyobject& direction, const pylist& geomap,
				 const pyobject& topography, const pylist& sections, 
			 const pydict& feature_types, const pydict& params, const pydict& pyunits): 
	Model(make_tuple(std::tuple<double, double, double>(python::extract<double>(bbox[0]),python::extract<double>(bbox[1]),python::extract<double>(bbox[2])), 
			 std::tuple<double, double, double>(python::extract<double>(bbox[3]),python::extract<double>(bbox[4]),python::extract<double>(bbox[5]))),
		  make_tuple(std::tuple<double, double, double>(python::extract<double>(abbox[0]),python::extract<double>(abbox[1]),python::extract<double>(abbox[2])), 
			 std::tuple<double, double, double>(python::extract<double>(abbox[3]),python::extract<double>(abbox[4]),python::extract<double>(abbox[5]))),
		point2( python::extract<double>(base_point[0]), python::extract<double>(base_point[1]) ),
		point2( python::extract<double>(direction[0]),  python::extract<double>(direction[1]) ) )
{

	for ( int i = 0; i < python::len( sections ); i++ ) {
		double cut = python::extract<double>( sections[i][1] );
		pytuple sbbox = calculate_section_bbox( bbox, base_point, direction, cut );
		pylist pl = python::extract<pylist>( sections[i] );
		pl.insert( 2, sbbox );
	}
	{
		pylist pl = geomap;
		pytuple sbbox = python::make_tuple(bbox[0], bbox[1], bbox[3], bbox[4]);
		pl.insert( 0, sbbox );
	}

	this->fill_model( geomap, topography, sections, feature_types, params, pyunits);
}

// Model python for HORIZONTAL cross sections.
ModelPython::ModelPython( const pyobject& bbox, const pyobject& abbox, const pylist& geomap, 
			  const pyobject& topography, const pylist& sections,
			  const pydict& feature_types, const pydict& params, const pydict& pyunits): 
	Model(make_tuple(std::tuple<double, double, double>(python::extract<double>(bbox[0]),python::extract<double>(bbox[1]),python::extract<double>(bbox[2])), 
			 std::tuple<double, double, double>(python::extract<double>(bbox[3]),python::extract<double>(bbox[4]),python::extract<double>(bbox[5]))),
		  make_tuple(std::tuple<double, double, double>(python::extract<double>(abbox[0]),python::extract<double>(abbox[1]),python::extract<double>(abbox[2])), 
			 std::tuple<double, double, double>(python::extract<double>(abbox[3]),python::extract<double>(abbox[4]),python::extract<double>(abbox[5]))))
{
	pytuple sbbox = python::make_tuple( bbox[0], bbox[1], bbox[3], bbox[4] );
	for ( int i = 0; i < python::len( sections ); i++ ) {
		pylist pl = python::extract<pylist>( sections[i] );
		pl.insert( 2, sbbox );
	}
	{
		pylist pl = geomap;
		pl.insert( 0, sbbox );
	}
	this->fill_model( geomap, topography, sections, feature_types, params, pyunits );
}


void ModelPython::make_matches() {
	// map<wstring, vector<triangle_pt>> faults = ((Model *)this)->make_matches();
	((Model *)this)->make_matches();
	
}

wstring ModelPython::mode(){
	return ((Model *)this)->mode;
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
			auto itf = this->feature_types.find(it->first);
			// Avoid seg fault, but should not happen.
			if ( itf == this->feature_types.end() ) {
				continue;
			}
			wstring lft = itf->second;
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
		wstring name1 = python::extract<wstring>(matching[i][0][0]);
		wstring name2 = python::extract<wstring>(matching[i][0][1]);
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

pylist ModelPython::pyabbox( ) const {
	pylist p;
	p.append( g0( g0( this->abbox ) ) );
	p.append( g1( g0( this->abbox ) ) );
	p.append( g2( g0( this->abbox ) ) );
	p.append( g0( g1( this->abbox ) ) );
	p.append( g1( g1( this->abbox ) ) );
	p.append( g2( g1( this->abbox ) ) );
	return p;
}

double ModelPython::signed_distance( const wstring& unit, const pyobject& pt ) const {
	// return ((Model *)this)->signed_distance(unit, point3(python::extract<double>(pt[0]), python::extract<double>(pt[1]), python::extract<double>(pt[2])));
	return ((Model *)this)->signed_distance_model(unit, point3(python::extract<double>(pt[0]), python::extract<double>(pt[1]), python::extract<double>(pt[2])));
}

double ModelPython::signed_distance_bounded( const wstring& unit, const pyobject& pt ) const {
	return ((Model *)this)->signed_distance_bounded(unit, point3(python::extract<double>(pt[0]), python::extract<double>(pt[1]), python::extract<double>(pt[2])));
}


double ModelPython::signed_distance_unbounded( const wstring& unit, const pyobject& pt ) const {
	return ((Model *)this)->signed_distance_unbounded(unit, point3(python::extract<double>(pt[0]), python::extract<double>(pt[1]), python::extract<double>(pt[2])));
}


double ModelPython::geomodelr_distance( const wstring& unit, const pylist& point ) const{

	double x = python::extract<double>( point[0] );
	double y = python::extract<double>( point[1] );
	double z = python::extract<double>( point[2] );
	return ((Model *)this)->geomodelr_distance(unit,point3(x,y,z));
	
}

pylist ModelPython::get_polygon(const wstring sec, int pol_idx){
	auto input = ((Model *)this)->get_polygon(sec,pol_idx);

	pylist output;
	for (auto& it_point: input){
		pylist point; point.append(gx(it_point)); point.append(gy(it_point));
		output.append(point);
	}
	return output;
}

pylist ModelPython::get_fault(const wstring sec, int pol_idx){
	auto input = ((Model *)this)->get_fault(sec,pol_idx);

	pylist output;
	for (auto& it_point: input){
		pylist point; point.append(gx(it_point)); point.append(gy(it_point));
		output.append(point);
	}
	return output;
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

