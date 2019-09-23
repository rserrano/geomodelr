#include "limiter.hpp"
#include "model.hpp"
#include <iostream>
// Compulsory destructors.
Limiter::~Limiter() {
}

// BBox Limiter
//   Non aligned
BBoxLimiter::BBoxLimiter(const bbox3& bbox, const Model * model): bbox(bbox), model(model) {
}

BBoxLimiter::~BBoxLimiter() {
}


double BBoxLimiter::limit_signed_distance(const point3& pt, double sdist) const {
	bool outside = false;
	double odist = 0.0;
	double idist = -std::numeric_limits<double>::infinity(); // In the future, fix distance below too.
	
	double x = gx(pt);
	double y = gy(pt);
	double z = gz(pt);
	
	double minx = g0(g0(bbox));
	double miny = g1(g0(bbox));
	double minz = g2(g0(bbox));
	
	double maxx = g0(g1(bbox));
	double maxy = g1(g1(bbox));
	double maxz = this->model->height( point2( x, y ) );
	
	double dists[6] = { minx - x, miny - y, minz - z, x - maxx, y - maxy, z - maxz };
	
	for ( size_t i = 0; i < 6; i++ ) {
		if ( dists[i] >= 0 ) {
			outside = true;
			odist += dists[i]*dists[i];
		} else {
			idist = std::max(idist, dists[i]);
		}
	}
	
	if ( outside ) {
		return std::max(sdist, std::sqrt(odist));
	}
	return std::max(sdist, idist);
}

// Polygon Limiter
//   Non Aligned
PolygonLimiter::PolygonLimiter(const polygon& poly, double bottom, const Model * model):bottom(bottom), model(model) {
  ring& outer = this->limit.outer();
  for ( auto pt: poly.outer() ) {

    if ( this->lpoly.size() && geometry::distance(this->lpoly.back(), pt) < boost_tol ) {
			continue;
		}
    geometry::append(outer, pt);
    geometry::append(this->lpoly, pt);
  }
  if ( outer.size() < 3 ) {
		throw GeomodelrException("Limit polygon has less than 3 points.");
	}
  geometry::validity_failure_type failure;
  if ( not geometry::is_valid(this->limit, failure) ) {
    geometry::correct(this->limit);
		string reason;
		if ( not geometry::is_valid(this->limit, reason) ) {
      string str = "Limit polygon has less than 3 points.";
      str += reason;
		  throw GeomodelrException(str.c_str());
    }
	}
	if ( not geometry::is_simple(this->limit) ) {
		  throw GeomodelrException("Polygon is not simple");
	}
  geometry::append(this->lpoly, poly.outer()[0]);
}

PolygonLimiter::~PolygonLimiter() {
}

double PolygonLimiter::limit_signed_distance(const point3& pt, double sdist) const {
  
  // Calculate the distance to the boundary.
  // Horizontal distance to the boundary.
  point2 pt2( gx(pt), gy(pt) );
  double dh = geometry::distance(pt2, limit);
  if ( dh == 0.0 ) {
    dh = -geometry::distance(pt2, lpoly);
  }

  double minz = g2(g0(this->model->bbox));
  // Check if there's a boolean there.
  if ( std::isfinite(this->bottom) ) {
    minz = this->bottom;
  }
	double maxz = this->model->height( point2( gx(pt), gy(pt) ) );
  double zdists[2] = { gz(pt)-maxz, minz-gz(pt) };
  // Vertical distance to the boundary.
  double dv;
  if ( zdists[0] > 0 ) {
    dv = zdists[0];
  } else if ( zdists[1] > 0 ) {
    dv = zdists[1];
  } else {
    dv = std::max(zdists[0], zdists[1]);
  }
  // Distance boundary.
  double db; 
  if ( dh > 0 ) {
    if ( dv > 0 ) {
      db = std::sqrt( dh*dh + dv*dv );
    } else {
      db = dh;
    }
  } else if ( dv > 0 ) {
    db = dv;
  } else {
    db = std::max( dv, dh );
  }
   
  if ( sdist < 0 ) {
    if ( db > 0 ) {
      // Inside boundary.
      return db;
    }
    // Inside object and boundary. (max as in least negative)
    return std::max( sdist, db );
  }
  if ( db > 0 ) {
    // Outside both.
    return std::max( db, sdist );
  }
  // Outside object but inside boundary.
  return sdist;
}

// Topography
Topography::Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims ):
point(point), sample(sample), heights(dims[0]*dims[1])
{
	this->dims[0] = dims[0];
	this->dims[1] = dims[1];
}


Topography::Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims, const vector<double>& heights ):
Topography( point, sample, dims ) {
	int rows = dims[0];
	int cols = dims[1];

	for ( int i = 0; i < rows; i++ ) {
		for ( int j = 0; j < cols; j++ ) {
			this->heights[ i*this->dims[1]+j ] = heights[i*this->dims[1]+j];
		}
	}
	this->max = *std::max_element(this->heights.begin(), this->heights.end());
	this->min = *std::min_element(this->heights.begin(), this->heights.end());
}

double Topography::height(const point2& pt) const {
	point2 pos(gx(pt), gy(pt));
	geometry::subtract_point(pos, this->point);
	geometry::divide_point(pos, this->sample);

	int i = std::floor(gx(pos));
	int j = std::floor(gy(pos));

	if (i < 0){i = 0;}
	if (j < 0){j = 0;}
	if (i >= this->dims[0]-1){i = this->dims[0]-2;}
	if (j >= this->dims[1]-1){j = this->dims[1]-2;}
	
	double A = this->heights[i*dims[1] + j];
	// return A;
	double B = this->heights[(i+1)*dims[1] + j];
	double C = this->heights[(i+1)*dims[1] + j+1];
	double D = this->heights[i*dims[1] + j+1];
	
	double x = 2.0*(gx(pos)- double(i)) - 1.0;
	double y = 2.0*(gy(pos)- double(j)) - 1.0;
	
	x = std::max( std::min( x, 1.0), -1.0 );
	y = std::max( std::min( y, 1.0), -1.0 );

	double v1 = B + A + x*(B - A);
	double v2 = C + D + x*(C - D);
	
	return 0.25*(v2 + v1 + y*(v2 - v1));
}

double Topography::height_new(const point2& pt) const {
	point2 pos(gx(pt), gy(pt));
	geometry::subtract_point(pos, this->point);
	geometry::divide_point(pos, this->sample);

	int i = std::floor(gx(pos));
	int j = std::floor(gy(pos));

	if (i < 0){i = 0;}
	if (j < 0){j = 0;}
	if (i >= this->dims[0]-1){i = this->dims[0]-2;}
	if (j >= this->dims[1]-1){j = this->dims[1]-2;}

	double A = this->heights[i*dims[1] + j];
	double B = this->heights[(i+1)*dims[1] + j];
	double C = this->heights[(i+1)*dims[1] + j+1];
	double D = this->heights[i*dims[1] + j+1];

	double alpha = gx(pos)- double(i);
	double betha = gy(pos)- double(j);
	double gamma = 1.0 - alpha - betha;
	if (gamma >= 0.0 ){
		assert(alpha>=-1e-5); assert(alpha<=1+1e-5);
		assert(betha>=-1e-5); assert(betha<=1+1e-5);
		assert(gamma>=-1e-5); assert(gamma<=1+1e-5);
		return gamma*A + alpha*B + betha*D;
	} else{
		alpha = 1.0 - alpha;
		betha = 1.0 - betha;
		gamma = 1.0 - alpha - betha;
		assert(alpha>=-1e-5); assert(alpha<=1+1e-5);
		assert(betha>=-1e-5); assert(betha<=1+1e-5);
		assert(gamma>=-1e-5); assert(gamma<=1+1e-5);
		return gamma*C + alpha*D + betha*B;
	}
}


bool check_intersection(const point3& pt, const point3& B, const point3& C){
	
	double det = gx(B)*gy(C) - gy(B)*gx(C);
	double alpha = (gy(C)*gx(pt) - gx(C)*gy(pt))/det;
	double beta  = (-gy(B)*gx(pt) + gx(B)*gy(pt))/det;
	double gamma = 1.0 - alpha - beta;

	return (std::abs(alpha-0.5)<=0.5) and (std::abs(beta-0.5)<=0.5) and (std::abs(gamma-0.5)<=0.5);

}


std::pair<bool, point3> projector(const point3& pt, const point3& projection, const point3& v1, const point3& v2){

	point3 n = geometry::cross_product(v1,v2);
	double dot = geometry::dot_product(projection, n);
	if (std::abs(dot) < epsilon){
		return std::make_pair(false, point3(0,0,0));
	}
	double t = -geometry::dot_product(n,pt)/dot;
	point3 output = point3( gx(pt) + t*gx(projection) , gy(pt) + t*gy(projection), gz(pt) + t*gz(projection));

	return std::make_pair(check_intersection( output, v1, v2 ), output);
}

std::pair<point3, double> Topography::intersection(const point3& pt, const point3& projection,
	const std::tuple<int, int, int, int>& limits ){
	
	int n = this->dims[1];
	const point2& x0 = this->point;
	double dx = gx(this->sample);
	double dy = gy(this->sample);

	std::pair<point3, double> output =  std::make_pair(point3(0,0,0), std::numeric_limits<double>::infinity() );
	if (std::get<0>(limits) == -1){
		return output;
	}

	for (int i = std::get<0>(limits); i <= std::get<2>(limits); i++){
		double x = gx(x0) + i*dx;
		for (int j = std::get<1>(limits); j <= std::get<3>(limits); j++){
			double y = gy(x0) + j*dy;

			std::wcerr << L"Size: " << this->heights.size() << L"\n";
			std::wcerr << L"(m, n): " << this->dims[0] << L"\t" << this->dims[1] << std::endl;
			double A =  this->heights[i*n + j];
			double B =  this->heights[(i+1)*n + j];
			double C =  this->heights[(i+1)*n + j+1];
			double D =  this->heights[i*n + j+1];

			point3 ft = point3(gx(pt)-x, gy(pt)-y, gz(pt)-A);
			auto pair = projector(ft, projection, point3(dx,0,B-A), point3(0,dy,D-A));
			if (pair.first){
				geometry::add_point(pair.second, point3(x,y,A));
				double distance = geometry::distance(pt, pair.second);
				if (distance < output.second){
					output = std::make_pair(pair.second, distance);
				}
			}

			ft = point3(gx(pt)-x-dx, gy(pt)-y-dy, gz(pt)-C);
			pair = projector(ft, projection, point3(-dx,0, D-C), point3(0,-dy,B-C));
			if (pair.first){
				geometry::add_point(pair.second, point3(x+dx, y+dy, C));
				double distance = geometry::distance(pt, pair.second);
				if (distance < output.second){
					output = std::make_pair(pair.second, distance);
				}
			}

		}
	}

	return output;
}


TopographyPython::TopographyPython( const pyobject& point, const pyobject& sample, 
                                    const pyobject& dims, const pylist& heights ):
                                    Topography(point2(python::extract<double>(point[0]),
                                                      python::extract<double>(point[1])),
	                                             point2(python::extract<double>(sample[0]), 
                                                      python::extract<double>(sample[1])),
		                                {python::extract<int>(dims[0]), python::extract<int>(dims[1])})
{
	int rows = python::len(heights);
	if ( rows != this->dims[0] ) {
		throw GeomodelrException("topography rows does not correspond with dims.");
	}
	for ( int i = 0; i < rows; i++ ) {
		int cols = python::len(heights[i]);
		if ( cols != this->dims[1] ) {
			throw GeomodelrException("topography columns does not correspond with dims.");
		}

		for ( int j = 0; j < cols; j++ ) {
			this->heights[i*this->dims[1]+j] = python::extract<double>(heights[i][j]);
		}
	}

	this->max = *std::max_element(this->heights.begin(), this->heights.end());
	this->min = *std::min_element(this->heights.begin(), this->heights.end());
	geometry::add_point(this->point, point2(0.5*gx(this->sample), 0.5*gy(this->sample)));
}

double TopographyPython::height( const pyobject& pypt ) const {
	return ((Topography *)this)->height( point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])) );
}

RestrictedFunction::RestrictedFunction( const pyobject& model, const wstring& restype, const pyobject& data ) {
  this->model = python::extract<const ModelPython *>(model);
  if ( restype == L"polygon" ) {

    polygon pbound;
    double bottom = -std::numeric_limits<double>::infinity();
    const pyobject& polygon = data[wstring(L"polygon")];
    size_t nnodes = python::len(polygon);
    ring& outer = pbound.outer();
    
    // Start filling the first ring.
    for ( size_t k = 0; k < nnodes; k++ ) {
    	pyobject pypt = polygon[k];
    	point2 aux = point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1]));
    	if ( outer.size() ) {
    		if ( geometry::distance(outer.back(), aux) < boost_tol ) {
    			continue;
    		}
    	}
    	outer.push_back(aux);
    }
    const pydict& bdict = python::extract<pydict>(data);
    if ( bdict.has_key(wstring(L"bottom")) ) {
      bottom = python::extract<double>(bdict[wstring(L"bottom")]);
    }
    this->limit.reset( new PolygonLimiter(pbound, bottom, (Model *)this->model) );
  } else if ( restype == L"bbox" ) {
    const pyobject& bbox = data[L"bbox"];
    double a = python::extract<double>(data[0]);
    double b = python::extract<double>(data[1]);
    double c = python::extract<double>(data[2]);
    auto min_bbox = std::make_tuple(a, b, c);
    a = python::extract<double>(data[3]);
    b = python::extract<double>(data[4]);
    c = python::extract<double>(data[5]);
    auto max_bbox = std::make_tuple(a, b, c);
    limit.reset( new BBoxLimiter(std::make_tuple(min_bbox, max_bbox), (Model *)this->model) );
  } else {
    throw GeomodelrException("Restriction not implemented.");
  }
}

double RestrictedFunction::signed_distance( const wstring& unit, const pyobject& point ) const {
  double a = python::extract<double>(point[0]), b = python::extract<double>(point[1]), c = python::extract<double>(point[2]);
  point3 pt(a, b, c);
  return this->model->signed_distance_bounded_restricted( unit, this->limit, pt );
}

