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
#include <cmath>


pyobject Match::pytriangulate = python::object();

Match::Match( const Section * a, const Section * b )
:a(a), b(b), faultidx(nullptr)
{

}

Match::~Match( )
{
	if ( this->faultidx != nullptr ) 
	{
		delete this->faultidx;
	}
}

void Match::set( const vector<std::pair<int, int>>& match ){
	this->a_to_b.clear();
	this->b_to_a.clear();
	for ( size_t i = 0; i < match.size(); i++ ) {
		a_to_b[match[i].first].push_back(match[i].second);
		b_to_a[match[i].second].push_back(match[i].first);
	}
}

pylist Match::get() const {
	std::set<std::pair<int, int>> ret;
	for ( auto it = this->a_to_b.begin(); it != this->a_to_b.end(); it++ ) {
		const vector<int>& m = it->second;
		for ( size_t j = 0; j < m.size(); j++ ) {
			ret.insert(std::make_pair(it->first, m[j]));
		}
	}
	pylist tret;
	for ( auto it = ret.begin(); it != ret.end(); it++ ) {
		tret.append(python::make_tuple(it->first, it->second));
	}
	return tret;
}

void Match::match_polygons() {
	map<wstring, vector<int>> units_a;
	map<wstring, vector<int>> units_b;
	for ( size_t i = 0; i < this->a->polygons.size(); i++ ) {
		units_a[this->a->units[i]].push_back(i);
	}
	for ( size_t i = 0; i < this->b->polygons.size(); i++ ) {
		units_b[this->b->units[i]].push_back(i);
	}
	vector<std::pair<int, int>> m;
	for ( auto it = units_a.begin(); it != units_a.end(); it++ ) {
		if ( units_b.find(it->first) != units_b.end() ) {
			vector<int>& pols_a = it->second;
			vector<int>& pols_b = units_b[it->first];
			for ( size_t i = 0; i < pols_a.size(); i++ )
			{
				for ( size_t j = 0; j < pols_b.size(); j++ ) {
					if ( geometry::intersects(this->a->polygons[pols_a[i]], this->b->polygons[pols_b[j]]) ) {
						m.push_back(std::make_pair(pols_a[i], pols_b[j]));
					}
				}
			}
		}
	}
	this->set( m );
}

void Match::load_polygons_match( const pylist& match ) {
	vector<std::pair<int, int>> vmatch;
	size_t nmatch = python::len(match);
	for ( size_t i = 0; i < nmatch; i++ ) {
		int a = python::extract<int>(match[i][0]);
		int b = python::extract<int>(match[i][1]);
		vmatch.push_back(std::make_pair(a, b));
	}
	this->set(vmatch);
}

int Match::crosses_triangles(const point2& point, double cut) const {
	if ( this->faultidx == nullptr ) {
		return 0;
	}
	for ( auto it = this->faultidx->qbegin( geometry::index::contains(point) ); it != this->faultidx->qend(); it++ ) {
		const auto fl = this->faults.find(g1(*it));
		const AlignedTriangle& tr = fl->second[g2(*it)];
		int side = tr.crosses_triangle(point, cut);
		if ( side != 0 ) 
		{
			return side;
		}
	}
	return 0;
}

AlignedTriangle::AlignedTriangle(const std::tuple<point3, point3, point3>& triangle) {
	const point3& p0 = g0(triangle);
	const point3& p1 = g1(triangle);
	const point3& p2 = g2(triangle);
	point3 v1 = p1;
	point3 v2 = p2;
	// Find vectors and its cross product.
	geometry::subtract_point(v1, p0);
	geometry::subtract_point(v2, p0);
	this->normal = point3(gy(v1)*gz(v2)-gz(v1)*gy(v2), gz(v1)*gx(v2)-gx(v1)*gz(v2), gx(v1)*gy(v2)-gy(v1)*gx(v2));
	// Find norm of the vector.
	double norm = std::sqrt(gx(this->normal)*gx(this->normal) + gy(this->normal)*gy(this->normal) + gz(this->normal)*gz(this->normal));
	// Reverse if below 0.0
	if ( gz(normal) < 0 ) {
		norm *= -1.0;
	}
	// divide by norm.
	geometry::divide_value(this->normal, norm);
	
	this->point  = p0;
	ring& outer = this->triangle.outer();
	outer.push_back(point2(gx(p0), gy(p0)));
	outer.push_back(point2(gx(p1), gy(p1)));
	outer.push_back(point2(gx(p2), gy(p2)));
}

int AlignedTriangle::crosses_triangle(const point2& point, double cut ) const {
	if ( geometry::within(point, this->triangle) ) {
		point3 pt(gx(point), gy(point), cut);
		geometry::subtract_point(pt, this->point);
		double d = geometry::dot_product(pt, this->normal);
		if ( d < 0 ) {
			return -1;
		} else {
			return 1;
		}
	}
	return 0;
}

point3 angles( const point3& t ) 
{
	/*
	Given a triangle it finds the angles of the triangle.
	*/
	point3 st(gx(t)*gx(t), gy(t)*gy(t), gz(t)*gz(t));
	point3 cs(std::acos((gy(st) + gz(st) - gx(st))/(2.0*gy(t)*gz(t))),
		  std::acos((gx(st) + gy(st) - gz(st))/(2.0*gx(t)*gy(t))),
		  0.0);
	geometry::set<2>(cs, M_PI - (gx(cs) + gy(cs)));
	return cs;
}

std::tuple<edge, edge> next_and_check(const triangle& tri, const edge& edg, size_t na)
{
	/*
	Gets the next edge and also the edge that now has been finished.
	*/
	if ( g1(tri) < int(na) ) {
		// Case g0(tri), g1(tri) are in the same edge g2(tri) the contrary.
		// (g0(tri), g1(tri)) need to be removed.
		if ( g0(edg) == g0(tri) ) {
			// Case g0(edg) is g0(tri), then next should be (g1(tri), g2(tri))
			return std::make_tuple( std::make_tuple(g1(tri), g2(tri)), std::make_tuple(g0(tri), g1(tri)) );
		} else {
			// Case g0(edg) is g1(tri), then next should be (g0(tri), g2(tri))
			return std::make_tuple( std::make_tuple(g0(tri), g2(tri)), std::make_tuple(g0(tri), g1(tri)) );
		}
	} else {
		// Case g1(tri), g2(tri) are in the same edge, g0(tri) the contrary.
		// (g1(tri), g2(tri)) need to be removed.
		if ( g1(edg) == g1(tri) ) {
		    // Case g0(edg) is g0(tri), then next should be (g1(tri), g2(tri))
		    return std::make_tuple( std::make_tuple(g0(tri), g2(tri)), std::make_tuple(g1(tri), g2(tri)) );
		} else {
		    // Case g0(edg) is g1(tri), then next should be (g0(tri), g2(tri))
		    return std::make_tuple( std::make_tuple(g0(tri), g1(tri)), std::make_tuple(g1(tri), g2(tri)) );
		}
	}
}

vector<triangle> test_start(const vector<triangle>& tris, const vector<point3>& pa, const vector<point3>& pb, const edge& start)
{
	int na = pa.size();
	auto pt = [&] ( int i ) -> const point3& {
		if ( i >= na ) {
			return pb[i-na];
		}
		return pa[i];
	};
	auto weight = [&] ( const triangle& tri ) {
		/*
		Minimum angle in triangle denotes the weight.
		*/
		point3 dt(geometry::distance(pt(g0(tri)), pt(g1(tri))), 
			  geometry::distance(pt(g1(tri)), pt(g2(tri))), 
			  geometry::distance(pt(g2(tri)), pt(g0(tri))));
		point3 ang = angles(dt);
		return std::min(gx(ang), std::min(gy(ang), gz(ang)));
	};
	auto weight_comp = [weight] ( const triangle& a, const triangle& b ) {
		double wa = weight(a);
		double wb = weight(b);
		if ( wa < wb ) {
			return true;
		}
		if ( wa > wb ) {
			return false;
		}
		return a < b;
	};
	auto in_triangle = []( int i, const triangle& tri ) {
		return (i == g0(tri) or i == g1(tri) or i == g2(tri));
	};
	auto has_edge = [in_triangle](const triangle& tri, const edge& edg) {
		return in_triangle(g0(edg), tri) and in_triangle(g1(edg), tri);
	};
	
	std::list<triangle> sorted_triangles( tris.begin(), tris.end() );
	sorted_triangles.sort( weight_comp );
	
	edge nxtedg = start;
	bool has_next;
	vector<triangle> output_triangles;
	
	do {
		has_next = false;
		for ( auto it = sorted_triangles.begin(); it != sorted_triangles.end(); it++ ) {
			if ( has_edge(*it, nxtedg) ) {
				output_triangles.push_back(*it);
				sorted_triangles.erase(it);
				std::tuple<edge, edge> edgs = next_and_check(*it, nxtedg, na);
				edge burnedg = g1(edgs);
				auto to_burn = [&](const triangle& trr) {
					return (has_edge(trr, burnedg) or has_edge(trr, nxtedg));
				};
				sorted_triangles.remove_if(to_burn);
				nxtedg = g0(edgs);
				has_next = true;
				break;
			}
		}
	} while ( has_next );
	if ( sorted_triangles.size() != 0 ) {
		throw GeomodelrException("The faults were not completely triangulated");
	}
	return output_triangles;
}

 
vector<triangle> faultplane_for_lines(const vector<point3>& l_a, const vector<point3>& l_b)
{
	/*
	Get the faults plane between lines la, lb.
	*/
	int na = l_a.size();
	int nb = l_b.size();
	
	vector<triangle> tris = Match::triangulate(l_a, l_b);
	
	// pt is the point depending on the index.
	auto pt = [&](int i) -> point3 {
		if ( i >= na ) {
			return l_b[i-na];
		} else {
			return l_a[i];
		}
	};
	
	auto pt_dist = [pt]( const edge& e ) -> double {
		return geometry::distance(pt(g0(e)), pt(g1(e)));
	};
	
	auto edge_comp = [pt_dist]( const edge& a, const edge& b ) -> bool {
		double da = pt_dist(a);
		double db = pt_dist(b);
		if ( da < db ) {
			return true;
		}
		if ( da > db ) {
			return false;
		}
		return a < b;
	};
	// Order the possible edges by weight. 
	
	std::set<edge, decltype(edge_comp)> pos_start( edge_comp );
	
	// Obtain the edges count in point distance order.
	for ( const triangle& tri:  tris ) {
		vector<edge> edgs;
		edgs.push_back(std::make_tuple(g0(tri), g2(tri)));
		if ( g1(tri) < na ) {
			edgs.push_back(std::make_tuple(g1(tri), g2(tri)));
		} else {
			edgs.push_back(std::make_tuple(g0(tri), g1(tri)));
		}
		for ( const edge& edg: edgs ) {
			if ( ( g0(edg) == 0 or g0(edg) == na-1) and (g1(edg) == na or g1(edg) == na+nb-1) ) {
        			pos_start.insert(edg);
			}
		}
	}
	if ( pos_start.size() == 0 ) {
		throw GeomodelrException("Could not find a configuration to start the triangulation.");
	}
	// Obtain the minimum count for the edges, (hopefully is one, but, well, weird cases).
	for ( const edge& start: pos_start ) {
		try {
			return test_start(tris, l_a, l_b, start);
		} catch ( const GeomodelrException& e ) {
		
		}
	}
	throw GeomodelrException("The faults were not completely triangulated");
}

pylist test_faultplane_for_lines(const pylist& pyla, const pylist& pylb) {
	auto pypoint = []( const pyobject& pt ) {
		return point3(python::extract<double>(pt[0]), 
			      python::extract<double>(pt[1]), 
			      python::extract<double>(pt[2]));
	};
	auto pytovct = [pypoint]( const pylist& pyl ) {
		vector<point3> l;
		for ( int i = 0; i < python::len(pyl); i++ ) {
			l.push_back(pypoint(python::extract<pyobject>(pyl[i])));
		}
		return l;
	};
	auto vcttopy = []( const vector<triangle>& l ) {
		pylist pyl;
		for ( size_t i = 0; i < l.size(); i++ ) {
			pyl.append(python::make_tuple(g0(l[i]),g1(l[i]),g2(l[i])));
		}
		return pyl;
	};
	vector<point3> la = pytovct(pyla);
	vector<point3> lb = pytovct(pylb);
	
	vector<triangle> res = faultplane_for_lines(la, lb);
	return vcttopy(res);
}

map<wstring, vector<triangle_pt>> Match::match_lines()
{
	/*
	Creates the fault planes given the cross sections with faults with the same name.
	*/
	// Create map of related faults.
	map<wstring, std::tuple<int, int>> rel_faults;
	for ( size_t i = 0; i < this->a->lines.size(); i++ )
	{
		const wstring& name = this->a->lnames[i];
		if ( name != L"" ) {
			rel_faults[name] = std::make_tuple(i, -1);
		}
	}
	for ( size_t i = 0; i < this->b->lines.size(); i++ ) {
		const wstring& name = this->b->lnames[i];
		if ( name != L"" and rel_faults.find(name) != rel_faults.end() ) {
			g1(rel_faults[name]) = i;
		}
	}
	
	for( auto it = rel_faults.begin(); it != rel_faults.end(); ) {
		if( g1(it->second) == -1 ) {
			it = rel_faults.erase(it);
		} else {
			++it;
		}
	}
	map<wstring, vector<triangle_pt>> retfaults;
	for ( auto it = rel_faults.begin(); it != rel_faults.end(); it++ ) {
		const line& la = this->a->lines[g0(it->second)];
		const line& lb = this->b->lines[g1(it->second)];
		
		const wstring& name = it->first;
		vector<point3> pa;
		vector<point3> pb;
		std::transform(la.begin(), la.end(), std::back_inserter(pa), [&]( const point2& p ) { return point3(gx(p), gy(p), this->a->cut); });
		std::transform(lb.begin(), lb.end(), std::back_inserter(pb), [&]( const point2& p ) { return point3(gx(p), gy(p), this->b->cut); });
		try {
			vector<triangle> fplane = faultplane_for_lines(pa, pb);
			size_t na = la.size();
			auto pt = [&] ( size_t n ) {
				if ( n < na ) {
					return pa[n]; 
				} else {
					return pb[n-na]; 
				}
			};
			
			std::transform(fplane.begin(), fplane.end(), std::back_inserter(retfaults[name]),
				[&] ( const triangle& t ) {
					return std::make_tuple(pt(g0(t)), pt(g1(t)), pt(g2(t))); 
				} );
			
			std::transform(retfaults[name].begin(), retfaults[name].end(), std::back_inserter(this->faults[name]),
				[&] ( const triangle_pt& t ) -> AlignedTriangle {
					return AlignedTriangle(std::make_tuple(g0(t), g1(t), g2(t))); 
				} );
		
		} catch ( const GeomodelrException& e ) {
			string aname(this->a->name.begin(), this->a->name.end());
			string bname(this->b->name.begin(), this->b->name.end());
			string sname(name.begin(), name.end());
			std::cerr << "could not interpolate fault " << sname << " between " 
				  << aname << " and " << bname << " " << e.what() << "\n"; 
		}
	}
	vector<value_f> envelopes;
	for ( auto it = this->faults.begin(); it != this->faults.end(); it++ ) {
		const auto& triangles = it->second;
		for ( size_t i = 0; i < triangles.size(); i++ ) {
			const auto& tr = triangles[i];
			box trbox;
			geometry::envelope(tr.triangle, trbox);
			envelopes.push_back(std::make_tuple(trbox, it->first, i));
		}
	}
	this->faultidx = new rtree_f(envelopes.begin(), envelopes.end());
	return retfaults;
}


void Match::load_triangulate() {
	if ( Match::pytriangulate.is_none() ) {
		pyobject shared = python::import("geomodelr.shared");
		Match::pytriangulate = shared.attr("opposite_triangles");
	}
}

vector<triangle> Match::triangulate( const vector<point3>& l_a, const vector<point3>& l_b ) {
	pylist points;
	for ( const auto& p: l_a ) {
		points.append(python::make_tuple(gx(p), gy(p), gz(p)));
	}
	int na = python::len(points);
	for ( const auto& p: l_b ) {
		points.append(python::make_tuple(gx(p), gy(p), gz(p)));
	}
	pyobject pytris = Match::pytriangulate(points, na);
	vector<triangle> tris;
	int nt = python::len(pytris);
	for ( int i = 0; i < nt; i++ ) {
		const pyobject& t = python::extract<pyobject>(pytris[i]);
		int t0 = python::extract<int>(t[0]);
		int t1 = python::extract<int>(t[1]);
		int t2 = python::extract<int>(t[2]);
		
		tris.push_back(std::make_tuple(t0, t1, t2));
	}
	return tris;
}

/*
def find_fault_plane_intersection(fplane, x0, v1, v2, nv, intersect_line):
    """
    Finds the lines that result after intersecting a fault plane with a plane.
    """
    lines = []
    for tp in fplane:
        ints = []
        for i in range(3):
            # Find the ray parameters.
            ni = (i+1)%3
            p0 = np.array(tp[i])
            p1 = np.array(tp[ni])
            nt = p1-p0
            d  = la.norm(nt)
            nt = nt/d
            
            # Find if the ray and the plane are parallel. (plane normal perpendicular).
            div = np.dot(nv, nt)
            if math.fabs(div) < 1e-9:
                continue
            
            # Find the intersection between the ray and the plane.
            a = np.dot(x0 - p0, nv)/div
            
            # Find if it lies within the line.
            if a > (d + 1e-9) or (a < -1e-9):
                continue
            
            # Find if the point is already in the list.
            p = p0 + nt*a
            for pp in ints:
                if la.norm(pp-p) < 1e-9:
                    break
            else:
                ints.append(p)
        
        # Size of ints must be 2 to append it as lines.
        if len(ints) == 2:
            lines.append(ints)
    
    filt_lines = []
    # Cut the lines to the plane.
    for line in lines:
        res = intersect_line(line)
        filt_lines += res
    
    # Join the lines.
    j_lines = []
    for l in filt_lines:
        # Check that the resulting lines have two coordinates.
        assert len(l) == 2
        for i in range(2):
            p = l[i]
            pc = l[(i+1)%2]
            br = False
            for pl in j_lines:
                fp = pl[0]
                if dist(fp, p) < 1e-9:
                    pl.insert(0, pc)
                    br = True
                    break
                lp = pl[-1]
                if dist(lp, p) < 1e-9:
                    pl.append(pc)
                    br = True
                    break
            if br:
                break
        else:
            j_lines.append(deepcopy(l))
    return j_lines
    
def find_faults_plane_intersection(fplanes, plane):
    """
    Finds the lines that intersect a plane with the fault planes.
    fplanes:    the set of fault planes to intersect with the plane.
    plane:      plane to intersect. It's the four corners of the plane.
    """
    # Find the vectors in the plane that generate the space.
    pn = map(np.array, plane)
    # First is the vector with two corners in the top.
    v1 = pn[1] - pn[0]
    # Second is last point and first point.
    v2 = pn[-1] - pn[0]
    v1 = v1/la.norm(v1)
    
    # Gram-Schmidt is used to make them perpendicular.
    v2 = v2 - np.dot(v1, v2)*v1
    v2 = v2/la.norm(v2)
    
    # Normal orthogonal to both.
    nv = np.cross(v1, v2)
    nv = nv/la.norm(nv)
    x0 = pn[0]
    
    # planep is a function that returns the plane coordinate in terms of v1, v2.
    planep = lambda p: [np.dot(v1, p-pn[0]), np.dot(v2, p-pn[0])]
    
    # Get the polygon terms of v1, v2.
    poly = map(planep, pn)
    shpoly = Polygon(poly)
    
    # Corverts line to point coordinates.
    to_coords = lambda line: map(list, line.coords)
    # Converts returned multigeometry to coordinate lines.
    sep_lines = lambda lg: map( to_coords, shared.shape_list(lg, 'LineString') )
    def intersect_line(line):
        """
        Returns line after cutting it with a polygon.
        """
        p0 = planep(line[0])
        p1 = planep(line[1])
        l = LineString([p0, p1])
        inter = shpoly.intersection(l)
        return sep_lines(inter)
     
    joined_lines = {}
    for nam, fplane in fplanes.iteritems():
        joined_lines[nam] = find_fault_plane_intersection(fplane, x0, v1, v2, nv, intersect_line)
    
    return joined_lines

def find_faults_multiple_planes_intersection(fplanes, planes):
    # All joined lines.
    joined_lines = {}
    
    # Coordinate in x to add to the plane lines.
    startx = 0.0
    for plane in planes:
        lines = find_faults_plane_intersection(fplanes, plane)
        for nam, lines in lines.iteritems():
            for line in lines:
                for point in line:
                    point[0] += startx
            if nam in joined_lines:
                joined_lines[nam] += lines
            else:
                joined_lines[nam] = lines
        # Find the vectors in the plane that generate the space.
        pn = map(np.array, plane)
        # First is the vector with two corners in the top.
        v1 = pn[1] - pn[0]
        # Update the coordinate in x to start.
        startx += la.norm(v1)
    return joined_lines
*/
/*
def faultplane_to_feature(fplane, name):
    """
    Converts a fault plane to a Feature.
    """
    return { 'type': 'Feature',
             'geometry': { 'type': 'Surface', 
                           'coordinates': fplane },
             'properties': { 'name': name } }

def faultplanes_to_featurecollection(fplanes):
    """
    Converts a set of faultplanes to a FeatureCollection
    """
    features = [ faultplane_to_feature(fplane, name) for name, fplane in fplanes.iteritems() ]
    return  { 
        'type': 'FeatureCollection',
        'geology_type': 'faults',
        'features': features,
        'transform': {
                    'type': 'identity'
                }
            }
def align_fault_with(faults, model_point):
    """
    It returns a structure that can be used to know if the triangle is before of after the point.
    """
    ret = {}
    for feature in faults['features']:
        triangles = []
        for triangle in feature['geometry']['coordinates']:
            triangles.append(AlignedTriangle(map(model_point, triangle)))
        ret[feature['properties']['name']] = triangles 
    return ret
*/

