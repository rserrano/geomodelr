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

void Topography::set_tree(int i_min, int i_max, int j_min, int j_max, node * leaf){

	leaf->i_min = i_min;
	leaf->i_max = i_max;
	leaf->j_min = j_min;
	leaf->j_max = j_max;

	int delta_i = i_max - i_min;
	int delta_j = j_max - j_min;

	// Last leafs of the tree
	if (delta_i < 8 and delta_j < 8){

		leaf->z_min = std::numeric_limits<double>::infinity();
		leaf->z_max = -std::numeric_limits<double>::infinity();

		for ( int i = i_min; i <= i_max; i++ ) {
			for ( int j = j_min; j <= j_max; j++ ) {
				leaf->z_min = std::min(this->heights[ i*this->dims[1]+j ], leaf->z_min);
				leaf->z_max = std::max(this->heights[ i*this->dims[1]+j ], leaf->z_max);
			}
		}
		leaf->left = NULL;
		leaf->right = NULL;

	}else{
		leaf->left = new node;
		leaf->right = new node;

		// Split with vertical line
		if (delta_i > delta_j){
			int m = delta_i/2;
			set_tree(i_min , i_min + m, j_min, j_max, leaf->left);
			set_tree(i_min + m , i_max, j_min, j_max, leaf->right);

		} else{
			int m = delta_j/2;
			set_tree(i_min , i_max, j_min, j_min + m, leaf->left);
			set_tree(i_min , i_max, j_min + m, j_max, leaf->right);

		}
		leaf->z_min = std::min( leaf->left->z_min, leaf->right->z_min);
		leaf->z_max = std::max( leaf->left->z_max, leaf->right->z_max);

	}
}

// Topography
Topography::Topography( const point2& point, const point2& sample, const std::array<int, 2>& dims ):
point(point), sample(sample), heights(dims[0]*dims[1])
{
	this->dims[0] = dims[0];
	this->dims[1] = dims[1];
	// this->topography_tree = NULL;
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
	
	this->topography_tree = new node;
	this->set_tree(0, rows-1, 0, cols-1, this->topography_tree);

}


Topography::~Topography(){
	this->destroy_tree( this->topography_tree );
}

void Topography::destroy_tree( node * tree ){
	if(tree != NULL){
		this->destroy_tree(tree->left);
		this->destroy_tree(tree->right);
		delete tree;
	}
}

void Topography::print_level(int k, const node * leaf, int level){

	if (leaf != NULL){
		if (k == level){
			std::wcerr << L"i values:" << leaf->i_min << L", " << leaf->i_max << "\n";
			std::wcerr << L"j values:" << leaf->j_min << L", " << leaf->j_max << "\n";
			std::wcerr << L"z values:" << leaf->z_min << L", " << leaf->z_max << "\n\n";
		}else{
			this->print_level(k+1, leaf->left, level);
			this->print_level(k+1, leaf->right, level);
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

std::map<double, point2> Topography::get_points_section(const point2& p0, const point2& pf){

	double x0 = gx(this->point);
	double y0 = gy(this->point);

	double dx = gx(this->sample);
	double dy = gy(this->sample);

	double xf = x0 + (this->dims[0] - 1)*dx;
	double yf = y0 + (this->dims[1] - 1)*dy;

	double diff_x = gx(pf) - gx(p0);
	double diff_y = gy(pf) - gy(p0);

	std::map<double, point2> output;

	// Vertical checker
	if (std::abs(diff_x) > epsilon ){
		for (int i = 0; i < this->dims[0]; i++){
			
			double x = x0 + i*dx;
			double t = (x - gx(p0))/diff_x;
			double y = gy(p0) + t*diff_y;
			
			if (y < y0 or y > yf){
				continue;
			}
			if (t < 0.0 or t > 1.0){
				continue;
			}

			output[t] = point2(x, y);
		}	
	}

	// Horizontal checker
	if (std::abs(diff_y) > epsilon ){
		for (int j = 0; j < this->dims[1]; j++){
			
			double y = y0 + j*dy;
			double t = (y - gy(p0))/diff_y;
			double x = gx(p0) + t*diff_x;
			
			if (x < x0 or x > xf){
				continue;
			}
			if (t < 0.0 or t > 1.0){
				continue;
			}

			output[t] = point2(x, y);
		}	
	}

	// Diagonal
	double b = dy*(gx(p0) - x0) + dx*(gy(p0) - y0);
	double a = dy*diff_x + dx*diff_y;
	for (int k=1; k <= (this->dims[0] + this->dims[1] - 3); k++){

		double t = (k*dx*dy - b)/a;
		double x = gx(p0) + t*diff_x;
		double y = gy(p0) + t*diff_y;
		if (x < x0 or x > xf or y < y0 or y > yf){
			continue;
		}
		if (t < 0.0 or t > 1.0){
			continue;
		}		
		output[t] = point2(x, y);
	}

	// Add first and end points.

	output[0.0] = p0;
	output[1.0] = pf;
	double norm = std::sqrt( diff_x*diff_x + diff_y*diff_y );

	// Clean nearest points
	auto it_1 = output.begin();
	auto it_2 = std::next(it_1);
	do{
		if ( (it_2->first - it_1->first)*norm < epsilon ){
			output.erase(it_2);
		} else{
			++it_1;
		}
		it_2 = std::next(it_1);

	} while ( it_2 != output.end() );
	
	return output;
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

bool box_intersection(const point3& lower, const point3& upper, const point3& pt, const point3& vec, const point3& abs_vec){

	point3 pt0 = point3( gx(pt)-0.5*(gx(lower) + gx(upper)), gy(pt)-0.5*(gy(lower) + gy(upper)),
		gz(pt)-0.5*(gz(lower) + gz(upper)));

	double dx = gx(upper) - gx(lower);
	double dy = gy(upper) - gy(lower);
	double dz = gz(upper) - gz(lower);

    double D = -gz(vec)*gy(pt0) + gy(vec)*gz(pt0);
    double r = 0.5*(gz(abs_vec)*dy + gy(abs_vec)*dz);
    if (D>r or D<-r){
        return false;
    }

    D = gz(vec)*gx(pt0) - gx(vec)*gz(pt0);
    r = 0.5*(gz(abs_vec)*dx + gx(abs_vec)*dz);
    if (D>r or D<-r){
        return false;
    }

    D = -gy(vec)*gx(pt0) + gx(vec)*gy(pt0);
    r = 0.5*(gy(abs_vec)*dx + gx(abs_vec)*dy);
    if (D>r or D<-r){
        return false;
    }

    return true;
}

void Topography::tree_intersection(std::pair<point3, double>& output, const point3& pt, const point3& vec,
	const point3& abs_vec, const node * leaf){

	point3 lower = point3(gx(this->point)+(leaf->i_min)*gx(this->sample), gy(this->point)+(leaf->j_min)*gy(this->sample),
		leaf->z_min) ;
	point3 upper = point3(gx(this->point)+(leaf->i_max)*gx(this->sample), gy(this->point)+(leaf->j_max)*gy(this->sample),
		leaf->z_max) ;

	if (box_intersection(lower, upper, pt, vec, abs_vec )){
		if (leaf->right == NULL){

			int n = this->dims[1];
			const point2& x0 = this->point;
			double dx = gx(this->sample);
			double dy = gy(this->sample);

			for (int i = leaf->i_min; i < leaf->i_max; i++){
				double x = gx(x0) + i*dx;
				for (int j = leaf->j_min; j < leaf->j_max; j++){
					double y = gy(x0) + j*dy;

					double A =  this->heights[i*n + j];
					double B =  this->heights[(i+1)*n + j];
					double C =  this->heights[(i+1)*n + j+1];
					double D =  this->heights[i*n + j+1];

					point3 ft = point3(gx(pt)-x, gy(pt)-y, gz(pt)-A);
					auto pair = projector(ft, vec, point3(dx,0,B-A), point3(0,dy,D-A));
					if (pair.first){
						geometry::add_point(pair.second, point3(x,y,A));
						double distance = geometry::distance(pt, pair.second);
						if (distance < output.second){
							output = std::make_pair(pair.second, distance);
						}
					}

					ft = point3(gx(pt)-x-dx, gy(pt)-y-dy, gz(pt)-C);
					pair = projector(ft, vec, point3(-dx,0, D-C), point3(0,-dy,B-C));
					if (pair.first){
						geometry::add_point(pair.second, point3(x+dx, y+dy, C));
						double distance = geometry::distance(pt, pair.second);
						if (distance < output.second){
							output = std::make_pair(pair.second, distance);
						}
					}

				}
			}
		} else{
			this->tree_intersection(output, pt, vec, abs_vec, leaf->left);
			this->tree_intersection(output, pt, vec, abs_vec, leaf->right);
		}
	}
}


std::pair<point3, double> Topography::intersects(const point3& pt, const point3& projection){

	std::pair<point3, double> output = std::make_pair(point3(0,0,0), std::numeric_limits<double>::infinity() );
	if (std::abs((gz(projection)))<epsilon){
		return output;
	}

	double t_max = (this->max - gz(pt))/gz(projection);
	point3 vec = point3(gx(pt) + t_max*gx(projection), gy(pt) + t_max*gy(projection), this->max);

	geometry::subtract_point(vec, pt);
	point3 abs_vec = point3(std::abs(gx(vec)), std::abs(gy(vec)), std::abs(gz(vec)));

	this->tree_intersection(output, pt, vec, abs_vec, this->topography_tree);

	return output;
}

std::tuple<int, int, int, int> Topography::square_limits(const point3& pt, const point3& projection) const{

	if (std::abs((gz(projection)))<epsilon){
		return std::make_tuple(-1,-1,0,0);
	}

	double pro_x = gx(projection), pro_y = gy(projection), pro_z = gz(projection);
	double t_max = (this->max - gz(pt))/pro_z;
	double t_min = (this->min - gz(pt))/pro_z;

	point3 pt_0 = point3(gx(pt) + t_min*pro_x, gy(pt) + t_min*pro_y, this->min);
	point3 pt_1 = point3(gx(pt) + t_max*pro_x, gy(pt) + t_max*pro_y, this->max);

	// Minimum limits
	point2 ft = point2( std::min(gx(pt_0), gx(pt_1)), std::min(gy(pt_0), gy(pt_1)) );
	geometry::subtract_point(ft, this->point);
	geometry::divide_point(ft, this->sample);

	int i_min = std::floor(gx(ft)), j_min = std::floor(gy(ft));

	if (i_min < 0){i_min = 0;}
	if (j_min < 0){j_min = 0;}
	if (i_min >= this->dims[0]-1){i_min = this->dims[0]-2;}
	if (j_min >= this->dims[1]-1){j_min = this->dims[1]-2;}

	// Maximum limits
	ft = point2( std::max(gx(pt_0), gx(pt_1)), std::max(gy(pt_0), gy(pt_1)) );
	geometry::subtract_point(ft, this->point);
	geometry::divide_point(ft, this->sample);

	int i_max = std::floor(gx(ft)), j_max = std::floor(gy(ft));

	if (i_max < 0){i_max = 0;}
	if (j_max < 0){j_max = 0;}
	if (i_max >= this->dims[0]-1){i_max = this->dims[0]-2;}
	if (j_max >= this->dims[1]-1){j_max = this->dims[1]-2;}

	return std::make_tuple(i_min, j_min, i_max, j_max);
}

std::pair<point3, double> Topography::intersection(const point3& pt, const point3& projection){

	std::tuple<int, int, int, int> limits = this->square_limits(pt, projection);
	
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

			// std::wcerr << L"Size: " << this->heights.size() << L"\n";
			// std::wcerr << L"(m, n): " << this->dims[0] << L"\t" << this->dims[1] << std::endl;
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

	this->topography_tree = new node;
	this->set_tree(0, this->dims[0]-1, 0, this->dims[1]-1, this->topography_tree);
}

void TopographyPython::print_level(int level){
	((Topography *)this)->print_level(0, this->topography_tree, level);
}


double TopographyPython::height( const pyobject& pypt ) const {
	return ((Topography *)this)->height_new( point2(python::extract<double>(pypt[0]), python::extract<double>(pypt[1])) );
}

void TopographyPython::info(){
	std::wcerr << L"Base point: " << geometry::wkt(this->point) << "\n";
	std::wcerr << L"Sample: " << gx(this->sample) << ", " << gy(this->sample) << "\n";
	std::wcerr << L"Dims: " << this->dims[0] << ", " << this->dims[1] << "\n";
}

pylist TopographyPython::get_points_section(const pyobject& p0_py, const pyobject& pf_py) const{

	point2 p0 = point2(python::extract<double>(p0_py[0]), python::extract<double>(p0_py[1]));
	point2 pf = point2(python::extract<double>(pf_py[0]), python::extract<double>(pf_py[1]));
	auto output = ((Topography *)this)->get_points_section(p0, pf);

	pylist points;
	for (auto it = output.begin(); it != output.end(); ++it){
		points.append( python::make_tuple( gx(it->second), gy(it->second)) );
	}
	return points;

}

pytuple TopographyPython::inter_1(const pyobject& pt, const pyobject& pro){
	double x = python::extract<double>(pt[0]);
	double y = python::extract<double>(pt[1]);
	double z = python::extract<double>(pt[2]);

	double px = python::extract<double>(pro[0]);
	double py = python::extract<double>(pro[1]);
	double pz = python::extract<double>(pro[2]);

	auto aux = this->intersection(point3(x,y,z), point3(px,py,pz));

	const point3& p = aux.first;
	return python::make_tuple(python::make_tuple(gx(p), gy(p), gz(p)), aux.second);
}

pytuple TopographyPython::inter_2(const pyobject& pt, const pyobject& pro){
	double x = python::extract<double>(pt[0]);
	double y = python::extract<double>(pt[1]);
	double z = python::extract<double>(pt[2]);

	double px = python::extract<double>(pro[0]);
	double py = python::extract<double>(pro[1]);
	double pz = python::extract<double>(pro[2]);

	auto aux = this->intersects(point3(x,y,z), point3(px,py,pz));

	const point3& p = aux.first;
	return python::make_tuple(python::make_tuple(gx(p), gy(p), gz(p)), aux.second);
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

