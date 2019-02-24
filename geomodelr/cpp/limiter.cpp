
#include "limiter.hpp"
#include "model.hpp"
// Compulsory destructors.
Limiter::~Limiter() {
}

AlignedLimiter::~AlignedLimiter() {
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

//   Aligned
BBoxAlignedLimiter::BBoxAlignedLimiter(const bbox3& bbox, const Model * model): abbox(bbox), model(model) {
}

//   Aligned
BBoxAlignedLimiter::~BBoxAlignedLimiter() {
}

double BBoxAlignedLimiter::limit_signed_distance(const point3& pt, double sdist) const {
  bool outside = false;
	double odist = 0.0;
	double idist = -std::numeric_limits<double>::infinity(); // In the future, fix distance below too.
	
	
	double minu = g0(g0(abbox));
	double minv = g1(g0(abbox));
	double minw = g2(g0(abbox));
	
	double maxu = g0(g1(abbox));
	double maxv = g1(g1(abbox));
	double maxw = g2(g1(abbox));
	
	double u = gx(pt);
	double v = gy(pt);
	double w = gz(pt);
	
	point3 in = this->model->inverse_point( point2( gx(pt), gy(pt) ), gz(pt) );
	double maxh = this->model->height( point2( gx(in), gy(in) ) );
	if ( this->model->horizontal ) {
		maxw = maxh;
	} else {
		maxv = maxh;
	}
	
	double dists[6] = { minu - u, minv - v, minw - w, u - maxu, v - maxv, w - maxw };
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
PolygonLimiter::PolygonLimiter(const polygon& poly, const Model * model): limit(poly), model(model) {
}

PolygonLimiter::~PolygonLimiter() {
}

double PolygonLimiter::limit_signed_distance(const point3& pt, double sds) const {
  return sds;
}

//   Aligned
PolygonAlignedLimiter::PolygonAlignedLimiter(const polygon& poly, const Model * model): alimit(poly), model(model) {
}

PolygonAlignedLimiter::~PolygonAlignedLimiter() {
}

double PolygonAlignedLimiter::limit_signed_distance(const point3& pt, double sds) const {
  return sds;
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
			this->heights[i*this->dims[1]+j] = heights[i*this->dims[1]+j];
		}
	}
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
}

double TopographyPython::height( const pyobject& pypt ) const {
	return ((Topography *)this)->height( point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])) );
}


