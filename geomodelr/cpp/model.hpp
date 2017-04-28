
class Model {
	point2 base_point;
	point2 direction;
	vector<Section *> sections;
	vector<Match *> match;
	vector<double> cuts;
	Topography * topography;
	
	std::pair<int, double> closest_match(bool a, int a_idx, int pol_idx, const point2& pt) const;
	
	struct Possible {
		Possible(int a_match, int b_match, double a_dist, double b_dist);
		int a_match;
		int b_match;
		double a_dist;
		double b_dist;
		bool operator<( const Possible& other ) const;
		double distance( double c ) const;
	};
	
	void clear_matches();
	std::pair<double, double> cuts_range;
	// Returns all the possible matches of this 2d point, given the distance is unknown.
	std::pair<point2, double> model_point(const point3& pt) const;
	point3 to_inverse_point(const point2& p, double cut) const;
	vector<Possible> get_candidates(size_t a_idx, const point2& pt_a, const point2& pt_b ) const;
	vector<Possible> possible_closest(size_t a_idx, const point2& pt_a, const point2& pt_b) const;
	std::tuple<int, int, double> closest_to( size_t a_idx, const point2& pt_a, const point2& pt_b, double cut ) const;
public:
	static bool verbose;
	
	// Methods to create matches or load them from files.
	map<wstring, vector<triangle_pt>> make_matches(); // Returns the faults in global coordinates, (at least until moving plane-fault intersection to C++).
	void set_matches(const vector< std::tuple< std::tuple<wstring, wstring>, vector<int> > >& matching);
	
	// Methods to query matches.
	pylist possible_closest(const pyobject& pt) const;
	pytuple model_point(const pyobject& pt) const;
	pytuple inverse_point(const pyobject& pt) const;
	pytuple closest(const pyobject& pt) const;
	pytuple closest_topo(const pyobject& pt) const;
	pydict info() const;
	double height(const pyobject& pt) const;
};

class ModelPython: Model {
	Model(const pyobject& cuts_range,
	      const pyobject& basepoint, const pyobject& direction, 
	      const pyobject& map, const pyobject& topography,
	      const pylist& sections);
	virtual ~Model();
	
	// Methods to create matches or load them from files.
	pydict make_matches(); // Returns the faults in global coordinates, (at least until moving plane-fault intersection to C++).
	void set_matches(const pylist& matching);
	pylist get_matches() const;
	// Methods to query matches.
	pylist possible_closest(const pyobject& pt) const;
	pytuple model_point(const pyobject& pt) const;
	pytuple inverse_point(const pyobject& pt) const;
	pytuple closest(const pyobject& pt) const;
	pytuple closest_topo(const pyobject& pt) const;
	pydict info() const;
	double height(const pyobject& pt) const;
}

class Topography {
	point2 point;
	point2 sample;
	int dims[2];
	vector<double> heights;
public:
	double height(const point2&) const;
	Topography( const point2& point, const point2& sample, const & dims, const vector<double> heights );
};

class TopographyPython : Topography {
	Topography( const pyobject& point, const pyobject& sample, const pyobject& dims, const pylist& heights );
}


