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

Match::Match( const vector<std::pair<int, int>>& match, size_t sa, size_t sb ): a_free(sa, true), b_free(sb, true) {
	for ( size_t i = 0; i < match.size(); i++ ) {
		this->a_free[match[i].first]  = false;
		this->b_free[match[i].second] = false;
		a_to_b[match[i].first].push_back(match[i].second);
		b_to_a[match[i].second].push_back(match[i].first);
	}
}

Match Model::make_match( const Section& a, const Section& b ) {
	map<string, vector<int>> units_a;
	map<string, vector<int>> units_b;
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
		int a = python::extract<int>(match[0]);
		int b = python::extract<int>(match[1]);
		vmatch.push_back(std::make_pair(a, b));
	}
	return Match(vmatch, sa, sb);
}

Model::Model( const pylist& basepoint, const pylist& direction, const pylist& sections) 
{
}

void Model::make_matches() {
	for ( size_t i = 1; i < this->sections.size(); i++ ) {
		this->match.push_back(this->make_match(this->sections[i-1], this->sections[i]));
	}
}

void Model::load_matches(const pylist& matching) {
	size_t nmatch = python::len(matching);
	for ( size_t i = 0; i < nmatch; i++ ) {
		const pylist& m = python::extract<pylist>(matching[i]);
		this->match.push_back(this->load_match(m, this->sections[i].polygons.size(), this->sections[i+1].polygons.size()));
	}
}

