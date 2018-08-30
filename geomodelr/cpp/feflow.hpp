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
#ifndef GEOMODELR_FEFLOW
#define GEOMODELR_FEFLOW

#include "basic.hpp"
#include <boost/dynamic_bitset.hpp>
#include <iostream>

class Model;

std::pair<triangMesh2D, vectorLayers > prismaticMesh(const Model* geo_model, const polygon& domain,
    map<wstring, std::pair<point2, double> >& points, const map<wstring, multi_line>& constraints,
    const vector<point2>& riverCorners, double triSize, double edgeSize, int num_layers, double rate,
    bool optimization, bool dist_alg);

#endif //GEOMODELR_FEFLOW
