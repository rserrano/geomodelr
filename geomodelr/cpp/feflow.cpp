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

triangMesh2D cdtToMesh(const CDT& cdt){

    vector<CDT::Point> points;
    points.reserve(cdt.number_of_vertices());

    for(CDT::Finite_vertices_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit)    { 
        // vertex_index[vit] = idx;
        // idx++;
        points.push_back(vit->point());
    }

    vector<openvdb::Vec3I> triangles;
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++){

        if (fit->is_in_domain()){
            triangles.push_back(openvdb::Vec3I(fit->vertex(0)->info(), fit->vertex(1)->info(), fit->vertex(2)->info()));
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
 const map<wstring, multi_line>& constraints){


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
        for (auto it = constraints.begin(); it != constraints.end(); it++){
            dis =  std::min(dis, geometry::distance(pt, it->second));
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
    for (int k = 1; k < numPts; k++){

        point2 target = polyLine[k];        
        double xf = gx(target), yf = gy(target);
        double dx = xf-x0, dy =  yf-y0;
        double norm = std::sqrt(dx*dx + dy*dy);
        int div = std::ceil(norm/edgeSize);

        polygonCGAL edge;
        edge.push_back(CDT::Point(x0,y0));
        for (int k=1; k < div; k++){
            edge.push_back(CDT::Point(x0 + (dx*k)/div, y0 + (dy*k)/div) );
        }
        edge.push_back(CDT::Point(xf,yf));
        cdt.insert_constraint(edge.vertices_begin(), edge.vertices_end(), false);
        x0 = xf; y0 = yf;
    }
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

vectorLayers regularGrid_root(const Model* geo_model, const CDT& cdt, int num_layers, double thickness){

    vectorLayers layers(num_layers + 1);
    size_t num_points = cdt.number_of_vertices();
    double bottom = g2(g0(geo_model->get_bbox())) + epsilon;
    // thickness = std::exp(rate/2.0);

    double a = 1.0 - 1.0/num_layers, log_a = std::log(a);
    double b = thickness/( g2(g1(geo_model->get_bbox())) - bottom), log_b = std::log(b);

    auto fx = [&a,&b](double p){
        return std::pow(a,p) + std::pow(b,p) - 1.0;
    };
    auto diff_fx = [&a, &b, &log_a, &log_b](double p){
        return  log_a*std::pow(a,p) + log_b*std::pow(b,p);
    };

    // double p = Newton_Raphson(fx, diff_fx, 2.0);
    std::cerr << "a: " << a << "\n";
    std::cerr << "b: " << b << "\n";
    std::cerr << "Root: " << p << "\n";

    auto power = [&thickness](const double& t){
        // return std::pow(t, rate);
        return std::pow(1.0 - std::pow(1.0 - t, 1.0/2.1), 2.1);
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

vectorLayers adaptiveGrid2(const Model* geo_model, const CDT& cdt, int num_layers){

    // std::map<CDT::Vertex_handle, Int32> vertex_index;
    size_t idx = 0;

    // for(CDT::Finite_vertices_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){

    //     CDT::Vertex_circulator incidents = cdt.incident_vertices(vit), done(incidents);
    //     std::cerr << "\nPoint " << vit->info() << ": " << vit->point() << "\n";
    //     if (incidents != 0) {
    //         do {
    //             std::cerr << incidents->info() << " -> " << incidents->point() << std::endl;
    //         } while(++incidents != done);
    //     }
    // }

    vectorLayers Layers;
    size_t num_points = cdt.number_of_vertices();
    double bottom = g2(g0(geo_model->get_bbox())) + epsilon;

    // Topography    
    vector<double> topography;
    topography.reserve(num_points);
    for(CDT::Finite_vertices_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit)    {
        topography.push_back (  geo_model->height(point2(vit->point().x(), vit->point().y())) );
    }
    Layers.push_back(topography);

    // Layer 1
    idx = 0;
    boost::dynamic_bitset<> fixed(num_points);
    vector<double> layer(num_points);

    if (num_layers > 1){
        for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){
            double xp = vit->point().x();
            double yp = vit->point().y();
            double z_min = (topography[idx]*(num_layers - 1) + bottom)/num_layers;
            auto pair = find_unit_limits_cpp(geo_model, xp, yp, topography[idx], z_min, epsilon);
            layer[idx] = pair.first;
            fixed[idx] = pair.second;
            idx++;
        }
        Layers.push_back(layer);
    }
    boost::dynamic_bitset<> fixed_test(num_points);

    for (int L = 2; L < num_layers; L++){
        idx = 0;
        for(CDT::Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit){
            double xp = vit->point().x();
            double yp = vit->point().y();
            double z_max = (topography[idx]*(num_layers - L + 1) + (L- 1)*bottom)/num_layers;
            double z_min = (topography[idx]*(num_layers - L) + L*bottom)/num_layers;
            auto pair = find_unit_limits_cpp(geo_model, xp, yp, z_max, z_min, epsilon);

            layer[idx] = pair.first;
            fixed_test[idx] = pair.second;
            if (not(fixed.test(idx))){
                Layers.back()[idx] = layer[idx];
            }
            idx++;
        }

        if ((fixed & fixed_test).any()){
            fixed = fixed & fixed_test;
            Layers.push_back(layer);
        } else{
            fixed = fixed | fixed_test;
        }
    }

    Layers.push_back( vector<double>(num_points, bottom) );

    return Layers;
}


std::pair<triangMesh2D, vectorLayers > prismaticMesh(const Model* geo_model, const polygon& domain,
    map<wstring, std::pair<point2, double> >& points, const map<wstring, multi_line>& constraints,
    const vector<point2>& riverCorners, double triSize, double edgeSize, int num_layers, double rate,
    bool optimization, bool dist_alg){

    CDT cdt;
    // Domain
    setDomain(cdt, domain, riverCorners);

    // Constraints (rivers)
    for (auto it = constraints.begin(); it != constraints.end(); it++){
        for (auto& polyLine: it->second){
            remeshEdge(cdt, polyLine, edgeSize);
        }
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

    if (dist_alg){
        // return std::make_pair(cdtToMesh(cdt), regularGrid_distance(geo_model, cdt, num_layers, rate) );
        return std::make_pair(cdtToMesh(cdt), adaptiveGrid2(geo_model, cdt, num_layers) );
    } else{
        return std::make_pair(cdtToMesh(cdt), regularGrid_root(geo_model, cdt, num_layers, rate) );
    }
}

