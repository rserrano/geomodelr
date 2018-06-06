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
#include "isosurfaces_vdb.hpp"
#include "model.hpp"

int sgn(float v) {
  return (v > 0) - (v < 0);
}

void set_sdfVoxels(GridType::Accessor& accessor, const vector<float>& vals,
    const openvdb::CoordBBox& indexBB, float bg){

    int C=0;
    for (Int32 i = indexBB.min().x(); i <= indexBB.max().x(); i+=2) {
        for (Int32 j = indexBB.min().y(); j <= indexBB.max().y(); j+=2) {
            for (Int32 k = indexBB.min().z(); k <= indexBB.max().z(); k+=2) {

                if (std::abs(vals[C])< bg){
                    accessor.setValue(openvdb::Coord(i,j,k),vals[C]);
                } else{
                    accessor.setValueOff(openvdb::Coord(i,j,k),vals[C]);
                }
                C++;
            }
        }
    }
}

vector<float> get_sdfVoxels(std::function<double (const point3&)>& signed_distance,
    double dx, double dy, double dz, double xi, double yi, double zi, const openvdb::CoordBBox& indexBB){

    vector<float> vals; vals.reserve(64);
    for (Int32 i = indexBB.min().x(); i <= indexBB.max().x(); i+=2) {
        for (Int32 j = indexBB.min().y(); j <= indexBB.max().y(); j+=2) {
            for (Int32 k = indexBB.min().z(); k <= indexBB.max().z(); k+=2) {
                vals.push_back(float( signed_distance( point3(dx*i + xi, dy*j + yi, dz*k + zi) )));
            }
        }
    }

    return vals;
}

GridType::Ptr get_openvdbGrid(std::function<double (const point3&)>& signed_distance, int grid_divisions, bbox3 xyzBB, int ndelta){

    // Initial parameters
    double dx = (g0(g1(xyzBB)) - g0(g0(xyzBB)))/double(grid_divisions-1);
    double dy = (g1(g1(xyzBB)) - g1(g0(xyzBB)))/double(grid_divisions-1);
    double dz = (g2(g1(xyzBB)) - g2(g0(xyzBB)))/double(grid_divisions-1);

    double xi = g0(g0(xyzBB))-ndelta*dx;
    double yi = g1(g0(xyzBB))-ndelta*dy;
    double zi = g2(g0(xyzBB))-ndelta*dz;

    grid_divisions += 2*ndelta;
    float bg = float(1.5*std::sqrt(dx*dx + dy*dy + dz*dz));

    // Define openvdb grid
    openvdb::CoordBBox indexBB(openvdb::Coord(0, 0, 0),
        openvdb::Coord(grid_divisions-1, grid_divisions-1, grid_divisions-1));

    GridType::Ptr grid = openvdb::FloatGrid::create(bg);
    typename GridType::Accessor accessor = grid->getAccessor();
    accessor.getTree()->fill(indexBB, bg);
    bg = grid->background();

    // Iterates over all active voxels.
    bool change;
    do {
        change = false;
        for (GridType::ValueOnIter iter = grid->beginValueOn(); iter; ++iter){
            switch (iter.getDepth()) {
                
                // level 1: InternalNode
                case 2:{

                    iter.getBoundingBox(indexBB);
                    vector<float> vals = get_sdfVoxels(signed_distance, dx, dy, dz, xi, yi, zi, indexBB);

                    if (std::any_of(vals.cbegin(), vals.cend(), [bg](float v){ return std::abs(v)<bg; })){
                        set_sdfVoxels(accessor, vals, indexBB, bg);
                        change = true;
                    } else{
                        iter.setValue(sgn(vals[32])*bg);
                        iter.setValueOff();
                    }
                    break;
                }

                // level 0: LeafNode
                case 3:{
                    if (iter.getValue() == bg){
                        openvdb::Coord coord  = iter.getCoord();
                        float dist = signed_distance( point3(dx*coord.x() + xi, dy*coord.y() + yi, dz*coord.z() + zi) );
                        iter.setValue(dist);

                        if (std::abs(dist)>=bg){
                            iter.setValueOff();
                        }
                    }
                    break;
                }
            }
        }

    } while (change);

    grid->insertMeta("deltas", openvdb::Vec3DMetadata(openvdb::Vec3d(dx, dy, dz)));
    grid->insertMeta("corner", openvdb::Vec3DMetadata(openvdb::Vec3d(xi, yi, zi)));

    return grid;
}

void setLinear_Transforms(GridType::Ptr& grid){

    openvdb::math::Transform::Ptr linearTransform =
        openvdb::math::Transform::createLinearTransform( openvdb::math::Mat4d::identity() );
    linearTransform->postScale(grid->metaValue<openvdb::Vec3d>("deltas"));
    linearTransform->postTranslate(grid->metaValue<openvdb::Vec3d>("corner"));
    grid->setTransform(linearTransform);
}

void gridToMesh(GridType::Ptr& grid, vector<openvdb::Vec3s>& points, vector<openvdb::Vec3I>& triangles){

    int av = grid->activeVoxelCount();
    double adap =  std::max(0.0,0.05*(double(av)/512000.0 - 1.0));

    vector<openvdb::Vec4I> quads;
    openvdb::tools::volumeToMesh(*grid, points, triangles, quads, tolerance, adap);

    triangles.reserve(triangles.size() + 2*distance(quads.begin(),quads.end()));

    for (auto& quad: quads){
        triangles.push_back(openvdb::Vec3I(quad.x(),quad.y(),quad.z()));
        triangles.push_back(openvdb::Vec3I(quad.z(),quad.w(),quad.x()));
    }
}

unitMesh getIsosurface(const Model* geo_model, const wstring& unit, bool bounded, bool aligned, int grid_divisions){

    openvdb::initialize();

    std::function<double (const point3&)> signed_distance;
    bbox3 xyzBB;
    int ndelta = 4;

    if (aligned){
        if (bounded){
            signed_distance = std::bind(&Model::signed_distance_bounded_aligned, geo_model, unit, std::placeholders::_1);
        } else{
            signed_distance = std::bind(&Model::signed_distance_unbounded_aligned, geo_model, unit, std::placeholders::_1);
        }
        xyzBB = geo_model->get_abbox();
    } else{
        if (bounded){
            signed_distance = std::bind(&Model::signed_distance_bounded, geo_model, unit, std::placeholders::_1);
        } else{
            signed_distance = std::bind(&Model::signed_distance_unbounded, geo_model, unit, std::placeholders::_1);
        }
        xyzBB = geo_model->get_bbox();
    }

    GridType::Ptr grid = get_openvdbGrid(signed_distance, grid_divisions, xyzBB, ndelta);
    setLinear_Transforms(grid);

    vector<openvdb::Vec3s> points;
    vector<openvdb::Vec3I> triangles;   
    gridToMesh(grid, points, triangles);

    return std::make_pair(points,triangles);

}

