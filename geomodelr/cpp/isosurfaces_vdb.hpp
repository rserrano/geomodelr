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
#ifndef GEOMODELR_ISOSURFACES_VDB
#define GEOMODELR_ISOSURFACES_VDB

#include "basic.hpp"
#include <iostream>
#include <functional>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/GridTransformer.h>
// #include <openvdb/tools/ChangeBackground.h>

typedef openvdb::FloatGrid GridType;
typedef openvdb::v3_1::Int32 Int32;

class Model;

unitMesh getIsosurface(const Model* geo_model, const wstring& unit, bool bounded, bool aligned, int grid_divisions,
    bool activeResampler);

#endif //GEOMODELR_ISOSURFACES_VDB
