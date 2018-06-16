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
                    accessor.setValueOnly(openvdb::Coord(i,j,k),vals[C]);
                } else{
                    // accessor.setValueOnly(openvdb::Coord(i,j,k),sgn(vals[C])*bg);
                    accessor.setValueOnly(openvdb::Coord(i,j,k),vals[C]);
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

bool checkCorrectTile(vector<float>& vals, float bg){

    bool active = std::any_of(vals.cbegin(), vals.cend(), [bg](float v){ return std::abs(v)<bg; });
    bool any_positive = std::any_of(vals.cbegin(), vals.cend(), [bg](float v){ return v>0; });
    bool any_negative = std::any_of(vals.cbegin(), vals.cend(), [bg](float v){ return v<0; });
    return active || (any_positive && any_negative);
}


template<typename func>
GridType::Ptr get_openvdbGrid(std::function<double (const point3&)>& signed_distance, int grid_divisions, bbox3& xyzBB, int ndelta,
    const func& setSignedDistance){

    // Initial parameters
    double xi = g0(g0(xyzBB));
    double yi = g1(g0(xyzBB));
    double zi = g2(g0(xyzBB));

    double xf = g0(g1(xyzBB));
    double yf = g1(g1(xyzBB));
    double zf = g2(g1(xyzBB));

    double dx = (xf - xi)/double(grid_divisions-1);
    double dy = (yf - yi)/double(grid_divisions-1);
    double dz = (zf - zi)/double(grid_divisions-1);

    xi -= ndelta*dx; xf += ndelta*dx;
    yi -= ndelta*dy; yf += ndelta*dy;
    zi -= ndelta*dz; zf += ndelta*dz;
    xyzBB = std::make_tuple(std::make_tuple(xi + dx,yi + dy,zi + dz),std::make_tuple(xf - dx, yf - dy, zf - dz));
    setSignedDistance(xyzBB, signed_distance);

    grid_divisions += 2*ndelta;
    float bg = float(1.5*std::sqrt(dx*dx + dy*dy + dz*dz));

    // Define openvdb grid
    openvdb::CoordBBox indexBB(openvdb::Coord(0, 0, 0),
        openvdb::Coord(grid_divisions-1, grid_divisions-1, grid_divisions-1));

    GridType::Ptr grid = openvdb::FloatGrid::create(bg);
    grid->setGridClass(openvdb::GRID_LEVEL_SET);
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

                    if (checkCorrectTile(vals,bg)){
                        set_sdfVoxels(accessor, vals, indexBB, bg);
                        change = true;
                    } else{
                        iter.setValue(vals[32]);
                        // iter.setValue(sgn(vals[32])*bg);
                        // iter.setValueOff();
                    }
                    break;
                }

                // level 0: LeafNode
                case 3:{
                    if (iter.getValue() == bg){
                        openvdb::Coord coord  = iter.getCoord();
                        float dist = signed_distance( point3(dx*coord.x() + xi, dy*coord.y() + yi, dz*coord.z() + zi) );

                        if (std::abs(dist)<bg){ 
                            iter.setValue(dist);
                        } else{
                            // iter.setValue(sgn(dist)*bg);
                            iter.setValue(dist);
                            // iter.setValueOff();
                        }
                    }
                    break;
                }
                // Other cases
                default:{
                    openvdb::Coord coord = iter.getBoundingBox().min();
                    std::cerr << "Depth: " << iter.getDepth() << std::endl;
                    float dist = signed_distance( point3(dx*coord.x() + xi, dy*coord.y() + yi, dz*coord.z() + zi) );
                    if (std::abs(dist)<bg){
                        accessor.setValueOnly(coord,dist);
                    } else{
                        // accessor.setValueOnly(coord,sgn(dist)*bg);
                        accessor.setValueOnly(coord,dist);
                    }
                    change = true;
                    break;
                }
            }
        }

    } while (change);

    grid->insertMeta("deltas", openvdb::Vec3DMetadata(openvdb::Vec3d(dx, dy, dz)));
    grid->insertMeta("corner", openvdb::Vec3DMetadata(openvdb::Vec3d(xi, yi, zi)));
    grid->insertMeta("resample", openvdb::Vec3IMetadata(openvdb::Vec3i(1, 1, 1)));

    return grid;
}

void setLinear_Transforms(GridType::Ptr& grid){
    /*
    Assign the scale and translate transformation to the grid.

    Args:
        grid: tree of openvdb that represents the unit surface.
    */

    openvdb::math::Transform::Ptr linearTransform =
        openvdb::math::Transform::createLinearTransform( openvdb::math::Mat4d::identity() );
    linearTransform->postScale(grid->metaValue<openvdb::Vec3d>("deltas"));
    linearTransform->postTranslate(grid->metaValue<openvdb::Vec3d>("corner"));
    grid->setTransform(linearTransform);
}

void gridToMesh(const GridType::Ptr& grid, vector<openvdb::Vec3s>& points, vector<openvdb::Vec3I>& triangles){
    /*
    Creates a triangulation of the unit using a tree of openvdb.
    
    Args:
        grid: The tree of openvdb that represents the unit surface.

        points: points of the mesh (empty).

        triangles: elements of the mesh (empty).
    */

    int av = grid->activeVoxelCount();
    openvdb::Vec3i rs = grid->metaValue<openvdb::Vec3i>("resample");
    double adaptivity =  std::max(0.0,0.2*(double(av)/681472.0 - 1.0))/double(rs.x()*rs.y()*rs.z());
    std::cerr << "adaptivity: " << adaptivity << std::endl;

    vector<openvdb::Vec4I> quads;
    openvdb::tools::volumeToMesh(*grid, points, triangles, quads, tolerance, adaptivity);

    triangles.reserve(triangles.size() + 2*distance(quads.begin(),quads.end()));

    for (auto& quad: quads){
        triangles.push_back(openvdb::Vec3I(quad.x(),quad.y(),quad.z()));
        triangles.push_back(openvdb::Vec3I(quad.z(),quad.w(),quad.x()));
    }
}


vector<float> check_BboxSurface(const GridType::Ptr& grid, int ndelta, int grid_divisions){
    /*
    Checks if the bounding box of the object modeled is very small. When a geological unit covers a very
    small part of the model, it needs to be refined. The new bbox of the unit is returned to check that
    case.

    Args:
        grid: tree of openvdb that represents the unit surface.

        ndelta: number of voxels to use to increase the bounding box.

        grid_divisions: divisions of the grid in the X, Y and Z directions.

    */

    float inf = std::numeric_limits<float>::infinity();
    vector<float> output = {inf, inf, inf, -inf, -inf, -inf};
    
    for (GridType::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter){
        if (iter.getValue()<0.0){

            openvdb::Coord coord  = iter.getCoord();
            output[0] = std::min(float(coord.x()) , output[0]);
            output[1] = std::min(float(coord.y()) , output[1]);
            output[2] = std::min(float(coord.z()) , output[2]);

            output[3] = std::max(float(coord.x()) , output[3]);
            output[4] = std::max(float(coord.y()) , output[4]);
            output[5] = std::max(float(coord.z()) , output[5]);            
        }
    }

    if (output[0] > output[3]){
        throw GeomodelrException("This model does not contain the unit or the sample is too coarse");
    }
    output[0] = std::max(output[0]-2,float(ndelta));
    output[1] = std::max(output[1]-2,float(ndelta));
    output[2] = std::max(output[2]-2,float(ndelta));

    output[3] = std::min(output[3]+2,float(ndelta + grid_divisions-1));
    output[4] = std::min(output[4]+2,float(ndelta + grid_divisions-1));
    output[5] = std::min(output[5]+2,float(ndelta + grid_divisions-1));

    return output;
}
// void closeUnbounded()

void deleteBorders(const vector<openvdb::Vec3s>& points, vector<openvdb::Vec3I>& triangles, const bbox3& xyzBB,
    openvdb::Vec3d& dL, bool activeResampler){

    double min_x = g0(g0(xyzBB)) + epsilon + 0.6*dL.x();
    double min_y = g1(g0(xyzBB)) + epsilon + 0.6*dL.y();
    double min_z = g2(g0(xyzBB)) + epsilon + 0.6*dL.z();

    double max_x = g0(g1(xyzBB)) - epsilon - 0.6*dL.x();
    double max_y = g1(g1(xyzBB)) - epsilon - 0.6*dL.y();
    double max_z = g2(g1(xyzBB)) - epsilon - 0.6*dL.z();

    auto check_outsidePT = [min_x,min_y,min_z,max_x,max_y,max_z](const openvdb::Vec3s& pt){
        bool x_out = (pt.x()<=min_x) || (pt.x()>=max_x);
        bool y_out = (pt.y()<=min_y) || (pt.y()>=max_y);
        bool z_out = (pt.z()<=min_z) || (pt.z()>=max_z);

        return (x_out || y_out || z_out);
    };

    std::function<bool (const openvdb::Vec3I&)> check_outsideTR;

    if (activeResampler){
        check_outsideTR = [&points, &check_outsidePT](const openvdb::Vec3I& tri){
            vector<openvdb::Vec3s> pts = {points[tri.x()],points[tri.y()],points[tri.z()]};
            return std::any_of(pts.cbegin(), pts.cend(), check_outsidePT);
        };
    } else{
        check_outsideTR = [&points, &check_outsidePT](const openvdb::Vec3I& tri){
            vector<openvdb::Vec3s> pts = {points[tri.x()],points[tri.y()],points[tri.z()]};
            return std::all_of(pts.cbegin(), pts.cend(), check_outsidePT);
        };
    }

    triangles.erase(std::remove_if(triangles.begin(), triangles.end(), check_outsideTR),triangles.end());  

}

/* --------------------------------------------------
   ---------------- Begin: RESAMPLE -----------------
   -------------------------------------------------- */

bool x_checker(vector<float>& vals,double dl){
    bool b1 = (vals[0]*vals[1]> -epsilon) && (std::abs(vals[0]+vals[1])<dl);
    bool b2 = (vals[2]*vals[3]> -epsilon) && (std::abs(vals[2]+vals[3])<dl);
    bool b3 = (vals[4]*vals[5]> -epsilon) && (std::abs(vals[4]+vals[5])<dl);
    bool b4 = (vals[6]*vals[7]> -epsilon) && (std::abs(vals[6]+vals[7])<dl);
    return b1||b2||b3||b4;
}
bool y_checker(vector<float>& vals,double dl){
    bool b1 = (vals[0]*vals[2]> -epsilon) && (std::abs(vals[0]+vals[2])<dl);
    bool b2 = (vals[1]*vals[3]> -epsilon) && (std::abs(vals[1]+vals[3])<dl);
    bool b3 = (vals[4]*vals[6]> -epsilon) && (std::abs(vals[4]+vals[6])<dl);
    bool b4 = (vals[5]*vals[7]> -epsilon) && (std::abs(vals[5]+vals[7])<dl);
    return b1||b2||b3||b4;
}
bool z_checker(vector<float>& vals,double dl){
    bool b1 = (vals[0]*vals[4]> -epsilon) && (std::abs(vals[0]+vals[4])<dl);
    bool b2 = (vals[1]*vals[5]> -epsilon) && (std::abs(vals[1]+vals[5])<dl);
    bool b3 = (vals[2]*vals[6]> -epsilon) && (std::abs(vals[2]+vals[6])<dl);
    bool b4 = (vals[3]*vals[7]> -epsilon) && (std::abs(vals[3]+vals[7])<dl);
    return b1||b2||b3||b4;
}

bool checkResample(GridType::ConstAccessor& accessor_old, GridType::Accessor& accessor_new, openvdb::Coord& coord,
    double dx, double dy, double dz, int rs_x, int rs_y, int rs_z){

    vector<float> vals; vals.reserve(8);
    bool states = true;

    float aux_val;
    for (Int32 k = coord.z(); k <= coord.z()+1; k++) {
        for (Int32 j = coord.y(); j <= coord.y()+1; j++) {
            for (Int32 i = coord.x(); i <= coord.x()+1; i++) {                
                states = accessor_old.probeValue(openvdb::Coord(i,j,k),aux_val) and states;
                vals.push_back(aux_val);
                accessor_new.setValueOn(openvdb::Coord(i*rs_x, j*rs_y,k*rs_z),aux_val);
            }
        }
    }
    if (states){
        return x_checker(vals,dx)||y_checker(vals,dy)||z_checker(vals,dz);
    }
    return false;
}

GridType::Ptr resample_openvdbGrid(const GridType::Ptr& grid_old, std::function<double (const point3&)>& signed_distance,
    int rs_x, int rs_y, int rs_z){

    double dx = (grid_old->metaValue<openvdb::Vec3d>("deltas")).x();
    double dy = (grid_old->metaValue<openvdb::Vec3d>("deltas")).y();
    double dz = (grid_old->metaValue<openvdb::Vec3d>("deltas")).z();

    double xi = (grid_old->metaValue<openvdb::Vec3d>("corner")).x();
    double yi = (grid_old->metaValue<openvdb::Vec3d>("corner")).y();
    double zi = (grid_old->metaValue<openvdb::Vec3d>("corner")).z();

    double dx_new = dx/rs_x; double dy_new = dy/rs_y; double dz_new = dz/rs_z;

    float bg = grid_old->background();
    GridType::Ptr grid_new = openvdb::FloatGrid::create(bg);
    grid_new->setGridClass(openvdb::GRID_LEVEL_SET);

    typename GridType::Accessor accessor_new = grid_new->getAccessor();
    typename GridType::ConstAccessor accessor_old = grid_old->getConstAccessor();

    for (GridType::ValueOnCIter iter = grid_old->cbeginValueOn(); iter; ++iter){

        openvdb::Coord coord  = iter.getCoord();
        bool need_resample = checkResample(accessor_old, accessor_new, coord, dx, dy, dz, rs_x, rs_y, rs_z);

        if (need_resample){
            for (Int32 i = rs_x*coord.x(); i <= rs_x*coord.x()+rs_x; i++) {
                for (Int32 j = rs_y*coord.y(); j <= rs_y*coord.y()+rs_y; j++) {
                    for (Int32 k = rs_z*coord.z(); k <= rs_z*coord.z()+rs_z; k++) {

                        auto ijk =  openvdb::Coord(i,j,k);
                        if (accessor_new.getValue(ijk) == bg){                            
                            float dist = signed_distance( point3(dx_new*i + xi, dy_new*j + yi, dz_new*k + zi) );
                            accessor_new.setValueOn(ijk,dist);
                        }
                    }
                }
            }
        }
    }

    grid_new->insertMeta("deltas", openvdb::Vec3DMetadata(openvdb::Vec3d(dx_new, dy_new, dz_new)));
    grid_new->insertMeta("corner", openvdb::Vec3DMetadata(openvdb::Vec3d(xi, yi, zi)));
    setLinear_Transforms(grid_new);
    return grid_new;

}

openvdb::Vec3i getNumberResamples(const GridType::Ptr& grid){

    typename GridType::ConstAccessor accessor = grid->getConstAccessor();
    size_t nx = 0, ny = 0, nz = 0;

    auto deltas = grid->metaValue<openvdb::Vec3d>("deltas");
    double dx = deltas.x(), dy = deltas.y(), dz = deltas.z();

    for (GridType::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter){
        
        openvdb::Coord coord  = iter.getCoord();
        vector<bool> states; states.reserve(8);

        for (Int32 k = coord.z(); k <= coord.z()+1; k++) {
            for (Int32 j = coord.y(); j <= coord.y()+1; j++) {
                for (Int32 i = coord.x(); i <= coord.x()+1; i++) {
                    states.push_back(accessor.isValueOn(openvdb::Coord(i,j,k)));
               }
            }
        }

        if (std::all_of(states.cbegin(),states.cend(),[](bool b){return b; })) {
            
            float val0 = accessor.getValue(coord);

            float valx = accessor.getValue(openvdb::Coord(coord.x()+1,coord.y(),coord.z()));
            if ((val0*valx > 0) && (std::abs(val0 + valx)<dx)){
                nx++;
            }
            
            float valy = accessor.getValue(openvdb::Coord(coord.x(),coord.y()+1,coord.z()));
            if ((val0*valy > 0) && (std::abs(val0 + valy)<dy)){
                ny++;
            }

            float valz = accessor.getValue(openvdb::Coord(coord.x(),coord.y(),coord.z()+1));
            if ((val0*valz > 0) && (std::abs(val0 + valz)<dz)){
                nz++;
            }

        }
    }

    size_t total = nx + ny + nz;
    int rs_x = int(std::ceil(double(10*nx)/total));
    int rs_y = int(std::ceil(double(10*ny)/total));
    int rs_z = int(std::ceil(double(10*nz)/total));

    std::cerr << "Resamples: " << rs_x << " " << rs_y << " " << rs_z << "\n";
    return openvdb::Vec3i(rs_x, rs_y, rs_z);

}

GridType::Ptr resample_unitSurface(const GridType::Ptr& sourceGrid, std::function<double (const point3&)>& signed_distance){

    auto rs = getNumberResamples(sourceGrid);
    GridType::Ptr new_grid = resample_openvdbGrid(sourceGrid, signed_distance, rs.x(), rs.y(), rs.z());
    float bg = new_grid->background();
    new_grid->tree().prune();

    GridType::Ptr targetGrid = openvdb::FloatGrid::create(bg);
    targetGrid->setGridClass(openvdb::GRID_LEVEL_SET);
    targetGrid->insertMeta("deltas", openvdb::Vec3DMetadata((new_grid->metaValue<openvdb::Vec3d>("deltas"))));
    targetGrid->insertMeta("corner", openvdb::Vec3DMetadata((new_grid->metaValue<openvdb::Vec3d>("corner"))));
    setLinear_Transforms(targetGrid);

    // Get the source and target grids' index space to world space transforms.
    const openvdb::math::Transform
        &sourceXform = sourceGrid->transform(),
        &targetXform = targetGrid->transform();

    // Compute a source grid to target grid transform.
    openvdb::Mat4R xform =
        sourceXform.baseMap()->getAffineMap()->getMat4() *
        targetXform.baseMap()->getAffineMap()->getMat4().inverse();
    
    // Create the transformer.
    openvdb::tools::GridTransformer transformer(xform);
    
    // Resample using trilinear interpolation.
    transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(*sourceGrid, *targetGrid);
    
    // Prune the target tree for optimal sparsity.
    targetGrid->tree().prune();

    auto def_ValueState = [bg](openvdb::CombineArgs<float>& args){
        if (std::abs(args.a()-bg) > tolerance) {
            // Transfer the A value and its active state.
            args.setResult(args.a());
            args.setResultIsActive(args.aIsActive());
        } else {
            // Preserve the B value and its active state.
            args.setResult(args.b());
            args.setResultIsActive(args.bIsActive());
        }
    };

    typename GridType::ConstAccessor accessor = targetGrid->getConstAccessor();

    targetGrid->setName("target");
    new_grid->setName("new_grid");
    sourceGrid->setName("old_grid");

    openvdb::io::File file("/media/sf_CompartidaVB/OpenVDB/Petroleo/allGrids.vdb");
    // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(targetGrid);
    grids.push_back(new_grid);
    grids.push_back(sourceGrid);
    // Write out the contents of the container.
    file.write(grids);
    file.close();

    // targetGrid->tree().combineExtended(new_grid->tree(), def_ValueState);
    new_grid->tree().combineExtended(targetGrid->tree(), def_ValueState);
    new_grid->insertMeta("resample", openvdb::Vec3IMetadata(rs));

    return new_grid;
}
/* --------------------------------------------------
   ----------------- End: RESAMPLE ------------------
   -------------------------------------------------- */   


unitMesh getIsosurface(const Model* geo_model, const wstring& unit, bool bounded, bool aligned, int grid_divisions,
    bool activeResampler){

    openvdb::initialize();

    std::function<double (const point3&)> signed_distance;
    bbox3 xyzBB;
    int ndelta = 4;

    if (aligned){
        if (bounded){
            signed_distance = std::bind(&Model::signed_distance_bounded_aligned, geo_model, unit, std::placeholders::_1);
        }
        xyzBB = geo_model->get_abbox();
    } else{
        if (bounded){
            signed_distance = std::bind(&Model::signed_distance_bounded, geo_model, unit, std::placeholders::_1);
        }
        xyzBB = geo_model->get_bbox();
    }

    //Define lambda function to get the correct signed function in case of undounded surface.
    auto setSignedDistance = [bounded,aligned,&geo_model,&unit](const bbox3& Bbox, std::function<double (const point3&)>& signed_distance){
        if (!bounded){
            if (aligned){
                signed_distance = std::bind(&Model::signed_distance_unbounded_aligned_restricted,
                    geo_model, unit, Bbox, std::placeholders::_1);
            } else{
                signed_distance = std::bind(&Model::signed_distance_unbounded_restricted,
                    geo_model, unit, Bbox, std::placeholders::_1);
            }
        }
    };

    // Creates a openvdb grid
    GridType::Ptr grid = get_openvdbGrid(signed_distance, grid_divisions, xyzBB, ndelta, setSignedDistance);

    // If the object is at least 8 times smaller than the full bbox, it will benefit lots from a thinner sample.
    vector<float> bb = check_BboxSurface(grid, ndelta, grid_divisions);
    float obj_cells = (bb[3]-bb[0])*(bb[4]-bb[1])*(bb[5]-bb[2]);
    if (float(pow(grid_divisions,3))/obj_cells > 6){

        std::cerr << "Malla burda.\n";
        openvdb::Vec3d ds = grid->metaValue<openvdb::Vec3d>("deltas");
        openvdb::Vec3d Xs = grid->metaValue<openvdb::Vec3d>("corner");
        auto min_Bbox = std::make_tuple(bb[0]*ds.x()+Xs.x(), bb[1]*ds.y()+Xs.y(), bb[2]*ds.z()+Xs.z());
        auto max_Bbox = std::make_tuple(bb[3]*ds.x()+Xs.x(), bb[4]*ds.y()+Xs.y(), bb[5]*ds.z()+Xs.z());
        xyzBB = std::make_tuple(min_Bbox,max_Bbox);
        grid = get_openvdbGrid(signed_distance, grid_divisions, xyzBB , ndelta, setSignedDistance);
    }

    // Assign the corresponding transformations: scale and translate.
    setLinear_Transforms(grid);

    // Resample grid.
    if (activeResampler){
        grid = resample_unitSurface(grid, signed_distance);
    }
    // Creates the triangular mesh.
    vector<openvdb::Vec3s> points;
    vector<openvdb::Vec3I> triangles;
    gridToMesh(grid, points, triangles);

    if (!bounded){
        deleteBorders(points, triangles, xyzBB, grid->metaValue<openvdb::Vec3d>("deltas"), activeResampler);
    }

    // Create a VDB file object.
    grid->setName("algorithm");
    openvdb::io::File file("/media/sf_CompartidaVB/OpenVDB/Petroleo/mygrids.vdb");
    // // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    // // Write out the contents of the container.
    file.write(grids);
    file.close();

    return std::make_pair(points,triangles);
}