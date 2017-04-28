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

class Section;
class Model;
class Topography;

struct GeomodelrException : std::runtime_error
{
	GeomodelrException(const string&);
};

bool always_true( const value& v ); 


pylist test_faultplane_for_lines(const pylist& pyla, const pylist& pylb);
std::tuple<point2, double> point_line_projection( const point2& p, const line& l );

