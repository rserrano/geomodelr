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

#include "section.hpp"

polymatch make_match( const Section& a, const Section& b ) {
	map<string, vector<int>> units_a;
	map<string, vector<int>> units_b;
	vector<bool> sel_a(a.polygons.size(), false);
	vector<bool> sel_b(b.polygons.size(), false);
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
						sel_a[pols_a[i]] = true;
						sel_b[pols_b[j]] = true;
						m.push_back(std::make_pair(pols_a[i], pols_b[j]));
					}
				}
			}
		}
	}
	vector<int> sa;
	for ( size_t i = 0; i < sel_a.size(); i++ ) {
		if ( not sel_a[i] ) {
			sa.push_back(i);
		}
	}
	vector<int> sb;
	for ( size_t i = 0; i < sel_b.size(); i++ ) {
		if ( not sel_b[i] ) {
			sb.push_back(i);
		}
	}
	return std::make_tuple( m, sa, sb );
}

class Model {
public:
	Model(const pyobject& geojson) {
	}
};

BOOST_PYTHON_MODULE(cpp)
{
	python::class_<Section>("Section", python::init<double, const pylist&, 
							const pylist&, const pylist&, 
							const pylist&, const pylist&>())
							.def("info", &Section::info)
							.def("closest", &Section::closest);
}

