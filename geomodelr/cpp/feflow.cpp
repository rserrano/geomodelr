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

#include "feflow.hpp"
#include "model.hpp"

std::vector<double> linspace(double a, double b, size_t N) {
    double h = (b - a) / double(N-1);
    std::vector<double> xs(N);
    typename std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

triangMesh2D cdtToMesh(const CDT& cdt, std::map<std::pair<size_t, size_t>, size_t>& edges_map){

    vector<CDT::Point> points;
    points.reserve(cdt.number_of_vertices());

    for(CDT::Finite_vertices_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit)    { 
        // vertex_index[vit] = idx;
        // idx++;
        points.push_back(vit->point());
    }

    vector<openvdb::Vec3I> triangles;
    size_t C = 0;
    size_t E = 0;
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++){

        if (fit->is_in_domain()){
            triangles.push_back(openvdb::Vec3I(fit->vertex(0)->info(), fit->vertex(1)->info(), fit->vertex(2)->info()));
            E++;
            for (size_t k=0; k<=2; k++){
                size_t idx_min = fit->vertex(k)->info();
                size_t idx_max = fit->vertex( (k+1)%3 )->info();
                if (idx_min > idx_max){
                    std::swap(idx_min, idx_max);
                }
                auto edge_key = std::make_pair(idx_min, idx_max);
                if (edges_map.find( edge_key ) == edges_map.end() ){
                    edges_map[edge_key] = C;
                    bool aux_1 = (idx_min != triangles.back().x() and idx_min != triangles.back().y() and idx_min != triangles.back().z());
                    bool aux_2 = (idx_max != triangles.back().x() and idx_max != triangles.back().y() and idx_max != triangles.back().z());
                    if (aux_1 or aux_2){
                        std::cerr << "Element: " << E << "\t" << triangles.back() << "\n";
                        std::cerr << "edge: " << C+1 << ". " << edge_key.first+1 << " - " << edge_key.second+1 << "\t";
                        std::cerr << idx_min+1 << "-" << idx_max+1 << "\n\n";
                    }
                    C++;
                }
            }
        }
    }
    return std::make_pair(points, triangles);
}

void splitEdge(CDT& cdt,const point2& X0, const point2& X1, const vector<point2>& riverCorners){

    double x0 = gx(X0);
    double delta = gx(X1) - x0;

    auto get_x = [](const point2& pt){
        return gx(pt);
    };
    auto get_y = [](const point2& pt){
        return gy(pt);
    };
    std::function<double (const point2&)> get_xy = std::bind(get_x, std::placeholders::_1);

    if (std::abs(delta)<tolerance){
        x0 = gy(X0);
        delta = gy(X1) - x0;
        get_xy = std::bind(get_y, std::placeholders::_1);
    }

    segment edge = segment(X0,X1);    
    vector<double> value_t;
    vector<point2> pts;
    for (auto& corner: riverCorners){
        double distance = geometry::distance(edge, corner);
        if (distance < epsilon){
            value_t.push_back( (get_xy(corner) - x0)/delta );
            pts.push_back(corner);
        }
    }
    if (pts.size() > 0){

        auto output = sort_indexes(value_t);
        double length = geometry::length(edge);
        double T = value_t[output[0]];
        if (length * T < epsilon){
            throw GeomodelrException("There is a line very close to a doamain point.");
        }

        point2 pt = pts[output[0]];
        cdt.insert_constraint(CDT::Point(gx(X0),gy(X0)), CDT::Point(gx(pt),gy(pt)));

        for (size_t k=1; k<pts.size(); k++){
            if (length * (value_t[output[k]] - T) < epsilon){
                throw GeomodelrException("There is a line very close to a doamain point.");
            }

            T = value_t[output[k]];
            point2 pt2 = pts[output[k]];
            cdt.insert_constraint(CDT::Point(gx(pt),gy(pt)), CDT::Point(gx(pt2),gy(pt2)));
            pt = pt2;
        }

        if (length * (1.0 - T) < epsilon){
            throw GeomodelrException("There is a line very close to a doamain point.");
        }
        cdt.insert_constraint(CDT::Point(gx(pt),gy(pt)), CDT::Point(gx(X1),gy(X1)));
    } else{
        cdt.insert_constraint(CDT::Point(gx(X0),gy(X0)), CDT::Point(gx(X1),gy(X1)));
    }
}

void setDomain(CDT& cdt, const polygon& domain, const vector<point2>& riverCorners){

    int nPoints =  domain.outer().size();
    const ring& outer = domain.outer();
    point2 X0 = outer[0];

    for (int k = 0; k < nPoints; k++){
        point2 X1 = outer[(k+1)%nPoints];
        splitEdge(cdt, X0, X1, riverCorners);
        X0 = X1;
    }
}

void remeshPoint(CDT& cdt, const point2& pt, double radio){

    double sqrt_3 = 0.8660254037844386;
    double sin[] = {0., sqrt_3, sqrt_3, 0, -sqrt_3, -sqrt_3};
    double cos[] = {1., 0.5, -0.5, -1., -0.5, 0.5};

    for (int k=0; k < 6; k++){
        cdt.insert_constraint(CDT::Point(gx(pt),gy(pt)) , CDT::Point( gx(pt) + radio*cos[k], gy(pt) + radio*sin[k] ));
    }
}

void remeshCorrector(map<wstring, std::pair<point2, double> >& points, const polygon& domain,
 const vector<value_s>& constraints){


    if (points.empty()){
        return;
    }
    // Comparison between points.
    for (auto i = points.begin(); i != std::prev(points.end()); i++){
        for (auto j = std::next(i,1); j != points.end(); j++){            
            double di = i->second.second;
            double dj = j->second.second;
            double remeshDis = di + dj + std::min(di,dj);
            double distance = geometry::distance(i->second.first, j->second.first);
            if (remeshDis > distance){
                i->second.second = (distance*di)/remeshDis;
                j->second.second = (distance*dj)/remeshDis;
            }
        }
    }

    for (auto i = points.begin(); i != points.end(); i++){
        double dis = std::numeric_limits<double>::infinity();
        point2 pt = i->second.first;

        // Iterates over multiline constraints.
        for (auto& it : constraints){
            dis =  std::min(dis, geometry::distance(pt, g0(it)));
        }

        // Iterates over segments of the polygon.
        geometry::for_each_segment(domain,
            [&dis, &pt](geometry::model::referring_segment<const point2> edge) {
                dis = std::min(dis, geometry::distance(edge, pt) );
            }
        );

        if (dis< epsilon){
            auto x = std::to_string(double(gx(pt)));
            auto y = std::to_string(double(gy(pt)));

            string exception = "Point [" + x + ", " + y + "] is very close to domain or a line constraint.";
            throw GeomodelrException(exception);
        }
        if (0.5*dis < i->second.second){
            i->second.second = 0.5*dis;
        }
    }
}

void remeshEdge(CDT& cdt, const line& polyLine, double edgeSize){

    int numPts = polyLine.size();
    point2 source = polyLine[0];
    double x0 = gx(source), y0 = gy(source);

    for (int j = 1; j < numPts; j++){

        point2 target = polyLine[j];        
        double xf = gx(target), yf = gy(target);
        double dx = xf-x0, dy =  yf-y0;
        double norm = std::sqrt(dx*dx + dy*dy);
        int div = std::ceil(norm/edgeSize);

        polygonCGAL edge;
        edge.push_back(CDT::Point(x0,y0));
        for (int k=1; k < div; k++){
            edge.push_back( CDT::Point(x0 + (dx*k)/div, y0 + (dy*k)/div) );
        }
        edge.push_back(CDT::Point(xf,yf));
        cdt.insert_constraint(edge.vertices_begin(), edge.vertices_end(), false);
        x0 = xf; y0 = yf;
    }
}
void remeshEdge2(CDT& cdt, const segment& edge, double edgeSize){

    double x0 = gx(edge.first), y0 = gy(edge.first);
    double xf = gx(edge.second), yf = gy(edge.second);
    double dx = xf-x0, dy =  yf-y0;
    double norm = std::sqrt(dx*dx + dy*dy);
    int div = std::ceil(norm/edgeSize);
    // std::cerr << geometry::wkt(edge.first) << "\n";
    // std::cerr << geometry::wkt(edge.second) << "\n\n";
    polygonCGAL polyEdge;
    polyEdge.push_back(CDT::Point(x0,y0));
    for (int k=1; k < div; k++){
        polyEdge.push_back( CDT::Point(x0 + (dx*k)/div, y0 + (dy*k)/div) );
    }
    polyEdge.push_back(CDT::Point(xf,yf));
    cdt.insert_constraint(polyEdge.vertices_begin(), polyEdge.vertices_end(), false);    
}

void segmentConstraint(const CDT& cdt, vector<size_t>& line_handle, const CDT::Vertex_handle v0,
    const CDT::Vertex_handle vf){

    
    double x0 = v0->point().x(), y0 = v0->point().y();
    double xf = vf->point().x(), yf = vf->point().y();
    double dx = xf - x0, dy = yf - y0;

    std::cerr << "\n\n" << v0->point() << " : " << v0->info() << " =========== ";
    std::cerr << vf->point() << " : " << vf->info() << "\n";
    vector<CDT::Vertex_handle> segment_handle;
    segment_handle.push_back(v0);
    size_t sizeSegment = 1;
    int C = 1;
    while (segment_handle.back()->info() != vf->info()){

        CDT::Edge_circulator incidents = cdt.incident_edges(segment_handle.back()), done(incidents);
        std::cerr << "COUN: " << C << "\n";

        if (incidents != 0) {
            do {
                if (cdt.is_constrained(*incidents)){
                    auto& f = *(incidents->first);
                    // CDT::Vertex_handle aa = f.vertex(f.cw(incidents->second));
                    CDT::Vertex_handle vt_h = f.vertex(f.ccw(incidents->second));
                    // std::cerr << aa->point() << " : " << aa->info() << "--\t--";
                    std::cerr << vt_h->point() << " : " << vt_h->info() << "\n";

                    if (vt_h->info() == vf->info()){
                        segment_handle.push_back(vt_h);
                        std::cerr << "final\n";
                        break;
                        // return;
                    }

                    if (vt_h->info() != line_handle.back()){
                        double xh = vt_h->point().x(), yh = vt_h->point().y();
                        double area = std::abs(0.5*(dx*(yh-y0) - dy*(xh-x0)));
                        std::cerr << "Area: " << area << "\n";
                        std::cerr << line_handle.back();
                        if (area < epsilon){
                            segment_handle.push_back(vt_h);
                            break;
                        }
                    }

                }
            } while(++incidents != done);
        }
        C++;
        if (sizeSegment == segment_handle.size() || C>100){
            throw GeomodelrException("PROBLEMSS...");
        }
        sizeSegment = segment_handle.size();
    }

    for (size_t k=1; k<segment_handle.size(); k++){
        line_handle.push_back(segment_handle[k]->info());
    }
   
}

std::map<wstring, vector<size_t>> getConstrainsEdges(const CDT& cdt, const rtree_s* constraints_tree,
    std::map<std::pair<size_t, size_t>, size_t>& edges_map){

    std::map<wstring, vector<size_t>> output;
    // Constraints
    for (auto eit = cdt.edges_begin(); eit != cdt.edges_end(); ++eit){
        if (cdt.is_constrained(*eit)){
            auto& f = *(eit->first);
            CDT::Vertex_handle pt_0 = f.vertex(f.cw(eit->second));
            CDT::Vertex_handle pt_f = f.vertex(f.ccw(eit->second));
            double x_min = std::min(pt_0->point().x(), pt_f->point().x());
            double x_max = std::max(pt_0->point().x(), pt_f->point().x());
            double y_min = std::min(pt_0->point().y(), pt_f->point().y());
            double y_max = std::max(pt_0->point().y(), pt_f->point().y());

            auto Bbox = box(point2(x_min - epsilon, y_min - epsilon), point2(x_max + epsilon, y_max + epsilon));
            auto X0 = point2(pt_0->point().x(), pt_0->point().y());
            auto Xf = point2(pt_f->point().x(), pt_f->point().y());

            auto test_segment = [&X0, &Xf](const value_s& v){
                return geometry::distance(X0,g0(v))<epsilon && geometry::distance(Xf,g0(v))<epsilon;
            };

            for (auto it = constraints_tree->qbegin( geometry::index::intersects(Bbox)&&
                geometry::index::satisfies(test_segment) ); it != constraints_tree->qend(); ++it){

                size_t idx_min = pt_0->info();
                size_t idx_max = pt_f->info();
                // output[g1(*it)].push_back( std::make_pair(pt_0->info(),pt_f->info()) );
                if (idx_min > idx_max){
                    std::swap(idx_min, idx_max);
                }
                output[g1(*it)].push_back(edges_map[std::make_pair(idx_min, idx_max)]);
                break;
            }
        }
    }
    // for (auto it = constraints_tree.)


    return output;
}

vectorLayers regularGrid_distance(const Model* geo_model, const CDT& cdt, int num_layers, double d){

    vectorLayers layers(num_layers + 1);
    int num_points = cdt.number_of_vertices();
    double bottom = g2(g0(geo_model->get_bbox())) + epsilon;
    
    // Topography
    auto& topography = layers.front();
    topography.reserve(num_points);
    for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit)    {
        topography.push_back (  geo_model->height(point2(vit->point().x(), vit->point().y())) );
    }

    // Layers
    for (int L = 1; L < num_layers; L++){
        vector<double>& layer = layers[L];
        layer.reserve(num_points); 
        for(double& topo_val: topography){
            double H = topo_val - bottom;
            // std::cerr << "d: " << d << "\n";
            // std::cerr << "other: " << 2.0*H/num_layers << "\n\n";
            double lthick = std::min(d, 2*H/num_layers - epsilon);
            double h = 2*(H - lthick*num_layers)/(num_layers*(num_layers - 1));
            layer.push_back(topo_val - lthick*L - L*(L-1)*h*0.5);
        }
    }
    layers[num_layers] = vector<double>(num_points, bottom);
    return layers;
}

vectorLayers regularGrid(const Model* geo_model, const CDT& cdt, int num_layers, double thickness){

    vectorLayers layers(num_layers + 1);
    size_t num_points = cdt.number_of_vertices();
    double bottom = g2(g0(geo_model->get_bbox())) + epsilon;

    // Define function fx and its derivative to use in the regular algorithm
    double a = 1.0 - 1.0/num_layers, log_a = std::log(a);
    double b = thickness/( g2(g1(geo_model->get_bbox())) - bottom), log_b = std::log(b);
    auto fx = [&a,&b](double p){
        return std::pow(a,p) + std::pow(b,p) - 1.0;
    };
    auto diff_fx = [&a, &b, &log_a, &log_b](double p){
        return  log_a*std::pow(a,p) + log_b*std::pow(b,p);
    };

    double p = Newton_Raphson(fx, diff_fx, 2.0);
    auto power = [&p](const double& t){
        return std::pow(1.0 - std::pow(1.0 - t, p), 1.0/p);
    };
    vector<double> vals_t = linspace(0.0, 1.0 , num_layers + 1);
    std::transform(vals_t.begin(), vals_t.end(), vals_t.begin(), power);
    
    // Topography
    vector<double>& topography = layers.front();
    topography.reserve(num_points);
    for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit)    {
        topography.push_back (  geo_model->height(point2(vit->point().x(), vit->point().y())) );
    }

    // Layers
    for (int L = 1; L < num_layers; L++){
        vector<double>& layer = layers[L];
        layer.reserve(num_points); 
        for(double& topo_val: topography){
            layer.push_back( topo_val + vals_t[L]*(bottom - topo_val) );
        }
    }
    layers[num_layers] = vector<double>(num_points, bottom);
    return layers;
}

vectorLayers adaptiveGrid(const Model* geo_model, const CDT& cdt, int num_layers, double rate){

    vectorLayers layers(num_layers + 1);
    int num_points = cdt.number_of_vertices();
    double bottom = g2(g0(geo_model->get_bbox())) + epsilon;

    // Resize
    for (auto& layer: layers){
        layer.resize(num_points);
    }
    
    // Topography
    int C = 0;
    for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){
        double xp = vit->point().x();
        double yp = vit->point().y();
        double topo_val = geo_model->height(point2(xp, yp));
        layers[0][C] = topo_val;
        double z_max = topo_val;

        vector<std::pair<double, int>> positions;
        positions.push_back( std::make_pair(topo_val, 0) );
        for (int L = 1; L < num_layers; L++){
          
            double z_min = (topo_val*(num_layers - L) + L*bottom)/num_layers;
            auto pair = find_unit_limits_cpp(geo_model, xp, yp, z_max, z_min, epsilon);
            if (pair.second){
                positions.push_back( std::make_pair(pair.first, L) );
            }
            z_max = z_min;
        }
        positions.push_back( std::make_pair(bottom, num_layers) );

        // std::cerr << "\nPoint: " << C << "\n";
        // std::cerr << "Size: " << positions.size() << "\n";
        for (size_t i = 0; i < positions.size()-1; i++){
            auto& A = positions[i];
            auto& B = positions[i+1];
            vector<double> vals = linspace(A.first, B.first, size_t(B.second - A.second + 1));
            for (int j = A.second + 1; j <= B.second; j++){
                layers[j][C] = vals[j - A.second - 1];
            }
        }
        C++;
    }

    return layers;
}

bool angle(double dx, double dy, double& dz, double Max_Tan){

    double norm = std::sqrt(dx*dx + dy*dy);
    double Tan_O = dz/norm;
    if (Tan_O>Max_Tan){
        dz = Max_Tan*norm;
        return true;
    }
    return false;

}

void layer_correction(const CDT& cdt, const vector<double>& prev_layer, vector<double>& layer, vector<bool>& fixed,
    vertexSet& fixed_vertex, double Max_Tan){

    for (auto it = fixed_vertex.begin(); it != fixed_vertex.end(); ++ it){

        CDT::Vertex_handle vertexH = it->second;
        CDT::Vertex_circulator incidents = cdt.incident_vertices(vertexH), done(incidents);
        double xp = vertexH->point().x(), yp = vertexH->point().y(), zp = -(it->first);

        if (incidents != 0) {
            do {
                if (not(cdt.is_infinite(incidents)) && not(fixed[incidents->info()])){
                    
                    double dx = xp - incidents->point().x();
                    double dy = yp - incidents->point().y();
                    size_t idx = incidents->info();
                    double dz = zp - layer[idx];
                    if (angle(dx, dy, dz, Max_Tan)){
                        layer[idx] = std::min(zp - dz, prev_layer[idx] -1.0 );
                        fixed[idx] = true;
                        fixed_vertex.insert( std::make_pair(-layer[idx], incidents ) );
                    }
                }
            } while(++incidents != done);
        }

    }


}


vectorLayers adaptiveGrid2(const Model* geo_model, const CDT& cdt, int num_layers, double Max_Tan){

    size_t idx = 0;
    vectorLayers Layers;
    size_t num_points = cdt.number_of_vertices();
    double bottom = g2(g0(geo_model->get_bbox())) + epsilon;

    // Topography    
    vector<double> topography;
    topography.reserve(num_points);
    for(CDT::Finite_vertices_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){
        topography.push_back (  geo_model->height(point2(vit->point().x(), vit->point().y())) );
    }
    Layers.push_back(topography);

    // Layer 1
    idx = 0;
    vector<bool> fixed(num_points);
    vector<double> layer(num_points);

    vertexSet fixed_vertex;
    vertexSet fixed_vertex_auxiliar;

    if (num_layers > 1){
        for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){
            double xp = vit->point().x();
            double yp = vit->point().y();
            double z_min = (topography[idx]*(num_layers - 1) + bottom)/num_layers;
            auto pair = find_unit_limits_cpp(geo_model, xp, yp, topography[idx], z_min, epsilon);
            layer[idx] = pair.first;
            fixed[idx] = pair.second;
            if (pair.second){
                fixed_vertex.insert( std::make_pair(-pair.first, vit) );
            }
            idx++;
        }
        Layers.push_back(layer);
    }

    layer_correction(cdt, Layers[0], Layers[1], fixed, fixed_vertex, Max_Tan);

    // Other layers
    vector<bool> fixed_test(num_points);
    for (int L = 2; L <= num_layers; L++){
        idx = 0;
        for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){
            double xp = vit->point().x();
            double yp = vit->point().y();
            double z_max = (topography[idx]*(num_layers - L + 1) + (L- 1)*bottom)/num_layers;
            double z_min = (topography[idx]*(num_layers - L) + L*bottom)/num_layers;
            if (L==num_layers){
                z_min = 0.5*(z_min + z_max);
            }
            auto pair = find_unit_limits_cpp(geo_model, xp, yp, z_max, z_min, epsilon);

            layer[idx] = pair.first;
            fixed_test[idx] = pair.second;

            if (not(fixed[idx])) {
                Layers.back()[idx] = pair.first;
                if (pair.second){
                    fixed_vertex.insert( std::make_pair(-pair.first, vit) );
                }
            } else{
                if (pair.second){
                    fixed_vertex_auxiliar.insert( std::make_pair(-pair.first, vit) );
                }
            }
            idx++;
        }

        layer_correction(cdt, Layers[Layers.size() - 2], Layers.back(), fixed, fixed_vertex, Max_Tan);
        vector<bool> result(num_points);
        std::transform(fixed.cbegin(),fixed.cend(),fixed_test.cbegin(),result.begin(), std::logical_and<bool>());

        if (std::any_of(result.cbegin(), result.cend(), [](bool elem){return elem;})){
            fixed = result;
            Layers.push_back(layer);
            fixed_vertex = fixed_vertex_auxiliar;
            fixed_vertex_auxiliar.clear();
        } else{

            // idx = -1;
            // auto change_values = [&layer, &idx](double val, bool elem){
            //     idx++;
            //     if (elem){
            //         return val;
            //     }
            //     return layer[idx];
            // };

            // std::transform(Layers.back().cbegin(),Layers.back().cend(),
            //     fixed.cbegin(),Layers.back().begin(), change_values);
            std::transform(fixed.cbegin(),fixed.cend(),fixed_test.cbegin(),fixed.begin(), std::logical_or<bool>());
        }
    }

    Layers.push_back( vector<double>(num_points, bottom) );

    return Layers;
}


feflowInfo prismaticMesh(const Model* geo_model, const polygon& domain,
    map<wstring, std::pair<point2, double> >& points, const vector<value_s>& constraints,
    const vector<point2>& riverCorners, double triSize, double edgeSize, int num_layers,
    double thickness, bool optimization, wstring algorithm, double Max_Tan){

    CDT cdt;
    // Domain
    setDomain(cdt, domain, riverCorners);

    // Constraints (rivers)
    for (auto& it: constraints){
        remeshEdge2(cdt, g0(it), edgeSize);
    }

    // Points (wells);
    remeshCorrector(points, domain, constraints);
    for (auto it = points.begin(); it != points.end(); it++){
        remeshPoint(cdt, it->second.first, it->second.second);
    }

    Mesher mesher(cdt);
    mesher.set_criteria(Criteria(0.125,triSize));
    mesher.init();
    mesher.refine_mesh();
    if (optimization){
        CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = 15);
        mesher.init();
        mesher.refine_mesh();
    }

    // Gabriel_Mesher gabriel(cdt);
    // std::cerr << "¿es Gabriel?  " << gabriel.is_conforming_Gabriel() << std::endl;
    // std::cerr << "¿es Delaunay?  " << gabriel.is_conforming_Delaunay();
    // gabriel.make_conforming_Gabriel();

    // for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); it++){
    //     const CGAL::Triangle_2<CGAL::Epick>& tri  = cdt.triangle(it);
    //     std::cout << " " << tri[0].x() << " " << tri[0].y() << " ";
    //     std::cout << " " << tri[1].x() << " " << tri[1].y() << " ";
    //     std::cout << " " << tri[2].x() << " " << tri[2].y() << " ";
    //     std::cout << std::endl << "-------------------" << std::endl;
    //     // std::cout << cdt.triangle(it->neighbor(0)) << std::endl;
    // }

    // std::cerr << "\nPoints:\n";
    // for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit)   
    // { 
    //     std::cerr << vit->point().x() << " " << vit->point().y() << std::endl;
    // }


    // // then make it conforming Gabriel
    // CGAL::make_conforming_Gabriel_2(cdt);
    // std::cout << "Number of vertices after make_conforming_Gabriel_2: "
    //         << cdt.number_of_vertices() << std::endl;

    // Set indices to the vertices
    size_t idx = 0;
    for(CDT::Finite_vertices_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){
        vit->info() = idx;
        idx++;
    }


    // for (auto it = output.begin(); it != output.end(); ++it){
    //     std::wcerr << it->first << "\n";
    //     for (auto& val: it->second){
    //         std::cerr << val.first << "\t " << val.second << "\n";
    //     }
    //     std::cerr << std::endl;
    // }

    rtree_s constraints_tree( constraints.begin(), constraints.end() );
    std::map<std::pair<size_t, size_t>, size_t> edges_map;
    triangMesh2D mesh_2D = cdtToMesh(cdt, edges_map);

    if (algorithm == L"adaptive"){
        return std::make_tuple(mesh_2D, adaptiveGrid2(geo_model, cdt, num_layers, Max_Tan),
               getConstrainsEdges(cdt, &constraints_tree, edges_map) );
    } else if (algorithm == L"regular"){
        return std::make_tuple(mesh_2D, regularGrid(geo_model, cdt, num_layers, thickness),
               getConstrainsEdges(cdt, &constraints_tree, edges_map) );
    } else{
        throw GeomodelrException("Algorithm must be equal to 'regular' or 'adaptive'.");
    }
}

