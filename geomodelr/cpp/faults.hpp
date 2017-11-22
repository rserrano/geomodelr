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
#ifndef GEOMODELR_FAULTS_HPP
#define GEOMODELR_FAULTS_HPP

#include "basic.hpp"

// pylist test_faultplane_for_lines(const pylist& pyla, const pylist& pylb);

pydict find_faults_plane_intersection(const pydict& fplanes, const pylist& plane);
pydict find_faults_multiple_planes_intersection(const pydict& fplanes, const pylist& planes);
 
#endif // GEOMODELR_FAULTS_HPP
