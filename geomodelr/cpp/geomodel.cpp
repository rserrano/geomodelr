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
#include <boost/python/args.hpp>
#include "geomodel.hpp"

GeomodelrException::GeomodelrException(const string& what) 
:std::runtime_error(what)
{
}

PyObject *GeomodelrExceptionType = NULL;

PyObject* createExceptionClass(const char* name, PyObject* baseTypeObj = PyExc_Exception)
{
	using std::string;
	namespace bp = boost::python;
	
	string scopeName = bp::extract<string>(bp::scope().attr("__name__"));
	string qualifiedName0 = scopeName + "." + name;
	char* qualifiedName1 = const_cast<char*>(qualifiedName0.c_str());
	
	PyObject* typeObj = PyErr_NewException(qualifiedName1, baseTypeObj, 0);
	if(!typeObj) bp::throw_error_already_set();
	bp::scope().attr(name) = bp::handle<>(bp::borrowed(typeObj));
	return typeObj;
}

void translate(GeomodelrException const& e)
{
	// Use the Python 'C' API to set up an exception object
	boost::python::object pythonExceptionInstance(e);
	PyErr_SetObject(GeomodelrExceptionType, pythonExceptionInstance.ptr());
	PyErr_SetString(GeomodelrExceptionType, e.what());
}

void set_verbose( bool verbose ) {
	geomodelr_verbose = verbose;
}

BOOST_PYTHON_MODULE(cpp)
{
	const char* doc_possible_closest =	"Given a point, it finds all the possible geological units in the given line for\n"
						"the adjacent cross sections.\n\n"
						"It can be used to query a grid aligned with the model faster, for purposes of\n"
						"generating a triangulated mesh or a grid.\n\n"
						"Args:\n"
						"    (tuple) point:\n"
						"    The three coordinates of the point in the given coordinate system.\n"
						"Returns:\n"
						"    (list(tuple(), ...)):\n"
						"    a list with all the possible units, each unit with the distance to the previous\n"
						"    and next cross sections.\n";
	
	
	const char* doc_closest = 	"Given a point, it finds the geological unit that's defined as the closest to\n"
					"that point.\n\n"
					"The basic definition of the algorithm is that, given a match between geological units,\n"
					"the distance from the point to the unit is the sum of the in-section distance to the point\n"
					"averaged by the distance to the cross section.\n\n"
					"Args:\n"
					"    (tuple) point:\n"
					"    The three coordinates of the point in the given coordinate system.\n"
					"Returns:\n"
					"    (tuple):\n"
					"    A tuple with the geological unit and the defined distance to that unit.";
	
	const char* doc_closest_topo =	"Same as closest but it returns (AIR, inf) if the point is above the topography.\n\n"
					"It first looks if the point is above the topography and returns (AIR, inf) in that case.\n"
					"Otherwise it returns the same as closest.\n\n"
					"Args:\n"
					"    (tuple) point:\n"
					"    The three coordinates (easting, northing, altitude a.s.l) of the point in the given coordinate system.\n"
					"Returns:\n"
					"    (tuple):\n"
					"    A tuple with the geological unit and the defined distance to that unit or AIR if it's\n"
					"    above the topography.";
	
	const char* doc_height ="Returns: the height at the given point at the topography.\n\n"
				"It returns the height at the point stored in the topography. In case the point it's outside the bounds of\n"
				"the model, it returns the height of the closest point inside.\n\n"
				"Args:\n"
				"    (tuple)point:\n"
				"    The two coordinates (easting, northing) of the point in the given coordinate system.\n"
				"Returns:\n"
				"    (real)\n"
				"    The height as stored in the topography.\n";
	
	const char* doc_model_point =	"Translates the point to internal coordinates\n\n"
					"It returns in-section x coordinate, in-section y coordinate, cut coordinate\n\n"
					"Args:\n"
					"    (tuple) point:\n"
					"    The three coordinates (esting, norting, altitute a.s.l) of the point in the given coordinate system.\n"
					"Returns:\n"
					"    (tuple)\n"
					"    The point in the internal coordinate system.\n";
	
	const char* doc_inverse_point =	"From internal coordinates, it returns the point in the given coordinate system.\n\n"
					"It returns easting, northing and altitude from in-section x coordinate, in-section\n"
					"y coordinate, cut coordinate\n\n"
					"Args:\n"
					"    (tuple) point:\n"
					"    The three coordinates of the internal point.\n"
					"Returns:\n"
					"    (tuple)\n"
					"    The point in the given coordinate system\n";
	
	const char* doc_verb =	"Sets the operations as verbose.\n\n"
				"When creating the model, it will advice the user of problems with geometries or matchings.\n\n"
				"Args:\n"
				"    (boolean) verbose:\n"
				"    if geomodelr should be verbose when creating the model.\n";
	
	const char* doc_signed_distance = "Given unit U and a point P, it finds the geomodelr distance to U minus\n"
					  "the geomodelr distance to the closest unit different to U\n\n"
					  "It returns a signed distance that's zero at the boundary of the unit,\n"
					  "negative inside the unit and possitive outside the unit\n\n"
					  "Args:\n"
					  "    (string) unit:\n"
					  "    The unit to measure the signed distance to\n\n"
					  "    (tuple) point:\n"
					  "    The three coordinates (esting, norting, altitute a.s.l) of the point in the given coordinate system.\n"
					  "Returns:\n"
					  "    (double)\n"
					  "    The signed distance from the unit to the point.\n";
	
	const char* doc_signed_distance_bounded = "Given unit U and a point P, it finds the geomodelr distance to U minus\n"
					  	  "the geomodelr distance to the closest unit different to U\n\n"
					          "It returns a signed distance that's zero at the boundary of the unit,\n"
					          "negative inside the unit and possitive outside the unit\n\n"
                                                  "unlike signed_distance, when the point is outside the bounds of the model,\n"
						  "or above the topography, it returns a positive number (outside)\n\n"
					          "Args:\n"
					          "    (string) unit:\n"
					          "    The unit to measure the signed distance to\n\n"
					          "    (tuple) point:\n"
					          "    The three coordinates (esting, norting, altitute a.s.l) of the point in the given coordinate system.\n"
					          "Returns:\n"
					          "    (double)\n"
					          "    The signed distance from the unit to the point.\n";
	
	const char* doc_signed_distance_unbounded = "Given unit U and a point P, it finds the geomodelr distance to U minus\n"
					  	  "the geomodelr distance to the closest unit different to U\n\n"
					          "It returns a signed distance that's zero at the boundary of the unit,\n"
					          "negative inside the unit and possitive outside the unit\n\n"
                                                  "unlike signed_distance unbounded, it just returns a positive number\n"
						  "when the point is above the topography. It does not always produce solids\n\n"
					          "Args:\n"
					          "    (string) unit:\n"
					          "    The unit to measure the signed distance to\n\n"
					          "    (tuple) point:\n"
					          "    The three coordinates (esting, norting, altitute a.s.l) of the point in the given coordinate system.\n"
					          "Returns:\n"
					          "    (double)\n"
					          "    The signed distance from the unit to the point.\n";
	
	const char* doc_intersect_planes = "Intersects a set of planes with the faults of the Geological Model.\n"
        				  "Takes a set of plane represented with its four corners and returns the set\n"
        				  "of lines that intersect that plane with the faults. The coordinates start from\n"
        				  "the first plane lower corner, and increase by dist(plane[i][0], plane[i][1]) for the\n"
        				  "next plane.\n\n"
        				  "Args:\n"
        				  "    (list) plane: List with planes. Each plane has a list with four corners\n"
        				  "    that we want to intersect the fault with.\n"
        				  "Returns:\n"
        				  "    (dict): a dictionary with fault names as keys, and lines, (list of points)\n"
        				  "    as values.\n";
        const char* doc_intersect_plane = "Intersects a plane with the faults of the Geological Model.\n\n"
        				  "Takes a plane represented with its four corners and returns the set\n"
        				  "of lines that intersect that plane with the faults.\n\n"
        				  "Args:\n"
        				  "    (list) plane: list with the four corners of the plane that we \n"
        				  "    want to intersect the fault with.\n\n"
        				  "Returns:\n"
        				  "    (dict): a dictionary with fault names as keys, and lines,\n"
        				  "    (list of points) as values. The coordinates go from the\n"
        				  "    lower left corner, (0.0, 0.0).\n";
	// Register exception.
	python::class_<GeomodelrException> GeomodelrExceptionClass("GeomodelrException", boost::python::init<std::string>());
	
	GeomodelrExceptionType = createExceptionClass("GeomodelrException");
	
	python::register_exception_translator<GeomodelrException>(&translate);
	
	python::def("faultplane_for_lines", test_faultplane_for_lines);
	
	// Register triangle-plane intersection. 
	//python::def("find_faults_plane_intersection", find_faults_plane_intersection);
	python::def("find_faults_intersection", find_faults_multiple_planes_intersection_python);
	python::def("topography_intersection", find_faults_topography_intersection_python);
	python::def("join_lines_tree_test",join_lines_tree_test);
	
	// Register bbox calculation for section.
	python::def("calculate_section_bbox", calculate_section_bbox );
	python::def("extend_line", test_extend_line );
	
	// Set verbose shows the errors in stderr.
	python::def("set_verbose", set_verbose, python::args("verbose"), doc_verb);

	// Single section class. Mainly exported for testing purposes.
	python::class_<SectionPython>("Section", python::init<const wstring&, double, const pytuple&,
							      const pylist&, const pylist&, const pylist&,
							      const pylist&, const pylist&, const pylist&>())
							      .def("info", &SectionPython::info)
							      .def("closest", &SectionPython::closest);
	
	// Main exported class, Model.
	python::class_<ModelPython>("Model", python::init<const pylist&, const pyobject&, const pyobject&, const pyobject&,
					         const pyobject&, pylist&, pydict&>())
					    .def(python::init<const pylist&, const pyobject&,
					    	 const pyobject&, pylist&, pydict&>())
					    .def("make_matches", &ModelPython::make_matches)
					    .def("possible_closest", &ModelPython::possible_closest, python::args("point"), doc_possible_closest)
					    .def("model_point", &ModelPython::model_point, python::args("point"), doc_model_point)
					    .def("inverse_point", &ModelPython::inverse_point, python::args("internal_point"), doc_inverse_point)
					    .def("closest", &ModelPython::closest, python::args("point"), doc_closest)
					    .def("closest_topo", &ModelPython::closest_topo, python::args("point"), doc_closest_topo)
					    .def("signed_distance", &ModelPython::signed_distance, python::args("unit", "point"), doc_signed_distance)
					    .def("signed_distance_bounded", &ModelPython::signed_distance_bounded, python::args("unit", "point"), doc_signed_distance_bounded)
					    .def("signed_distance_unbounded", &ModelPython::signed_distance_unbounded, python::args("unit", "point"), doc_signed_distance_unbounded)
					    .def("height", &ModelPython::height, python::args("point"), doc_height)
					    .def("intersect_plane", &ModelPython::intersect_plane, doc_intersect_plane)
					    .def("intersect_planes", &ModelPython::intersect_planes, doc_intersect_planes)
					    .def("intersect_topography", &ModelPython::intersect_topography)
					    .def("find_unit_limits", &ModelPython::find_unit_limits)
					    .def("info", &ModelPython::info)
					    .add_property("bbox", &ModelPython::pybbox)
					    .add_property("matches", &ModelPython::get_matches, &ModelPython::set_matches)
					    .add_property("lines", &ModelPython::get_lines)
					    .add_property("not_extended_lines", &ModelPython::get_not_extended_lines)
					    .add_property("faults", &ModelPython::get_faults)
					    .add_property("not_extended_faults", &ModelPython::get_not_extended_faults)
					    .add_property("fracts", &ModelPython::get_fracts)
					    .add_property("not_extended_fracts", &ModelPython::get_not_extended_fracts)
					    .add_property("veins", &ModelPython::get_veins)
					    .add_property("not_extended_veins", &ModelPython::get_not_extended_veins);

}

wstring human_failure_type( const geometry::validity_failure_type& fail )
{
	switch ( fail ) {
		case geometry::validity_failure_type::no_failure:
			return L"no failure";
		case geometry::validity_failure_type::failure_few_points:
			return L"failure few points";
		case geometry::validity_failure_type::failure_wrong_topological_dimension:
			return L"failure wrong topological dimension";
		case geometry::validity_failure_type::failure_spikes:
			return L"failure spikes";
		case geometry::validity_failure_type::failure_duplicate_points:
			return L"failure duplicate points";
		case geometry::validity_failure_type::failure_not_closed:
			return L"failure not closed";
		case geometry::validity_failure_type::failure_self_intersections:
			return L"failure self intersections";
		case geometry::validity_failure_type::failure_wrong_orientation:
			return L"failure wrong orientation";
		case geometry::validity_failure_type::failure_interior_rings_outside:
			return L"failure interior rings outside";
		case geometry::validity_failure_type::failure_nested_interior_rings:
			return L"failure nested interior rings";
		case geometry::validity_failure_type::failure_disconnected_interior:
			return L"failure disconnected interior";
		case geometry::validity_failure_type::failure_intersecting_interiors:
			return L"failure intersecting interiors";
		case geometry::validity_failure_type::failure_wrong_corner_order:
			return L"failure wrong corner order";
		default:
			return L"unknown";
	}
}
