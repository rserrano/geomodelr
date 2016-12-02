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
	Model::verbose = verbose;
}

BOOST_PYTHON_MODULE(cpp)
{
	const char* doc_possible_closest =	"Given a point, it finds all the possible geological units in the given line for\n"
						"the adjacent cross sections.\n\n"
						"It can be used to query a grid aligned with the model faster, for purposes of\n"
						"generating a triangulated mesh or a grid.\n\n"
						"Parameters\n"
						"---------\n"
						"point : tuple\n"
						"    The three coordinates of the point in the given coordinate system.\n"
						"Returns\n"
						"-------\n"
						"list(tuple(), ...)\n"
						"    a list with all the possible units, each unit with the distance to the previous\n"
						"    and next cross sections.\n";
	
	const char* doc_closest = 	"Given a point, it finds the geological unit that's defined as the closest to\n"
					"that point.\n\n"
					"The basic definition of the algorithm is that, given a match between geological units,\n"
					"the distance from the point to the unit is the sum of the in-section distance to the point\n"
					"averaged by the distance to the cross section.\n"
					"Parameters\n"
					"---------\n"
					"point : tuple\n"
					"    The three coordinates of the point in the given coordinate system.\n"
					"Returns\n"
					"-------\n"
					"tuple\n"
					"    A tuple with the geological unit and the defined distance to that unit.";
	
	const char* doc_closest_topo =	"Same as closest but it returns (AIR, inf) if the point is above the topography.\n\n"
					"It first looks if the point is above the topography and returns (AIR, inf) in that case.\n"
					"Otherwise it returns the same as closest.\n\n"
					"Parameters\n"
					"---------\n"
					"point : tuple\n"
					"    The three coordinates (easting, northing, altitude a.s.l) of the point in the given coordinate system.\n"
					"Returns\n"
					"-------\n"
					"tuple\n"
					"    A tuple with the geological unit and the defined distance to that unit or AIR if it's\n"
					"    above the topography.";
	
	const char* doc_height ="Returns the height at the given point at the topography.\n\n"
				"It returns the height at the point stored in the topography. In case the point it's outside the bounds of\n"
				"the model, it returns the height of the closest point inside.\n\n"
				"Parameters\n"
				"----------\n"
				"point tuple\n\n"
				"    The two coordinates (easting, northing) of the point in the given coordinate system.\n\n"
				"Returns\n"
				"-------\n"
				"real\n"
				"    The height as stored in the topography.\n";
	
	const char* doc_model_point =	"Translates the point to internal coordinates\n\n"
					"It returns in-section x coordinate, in-section y coordinate, cut coordinate\n\n"
					"Parameters\n"
					"----------\n"
					"point : tuple\n"
					"    The three coordinates (esting, norting, altitute a.s.l) of the point in the given coordinate system.\n"
					"Returns\n"
					"----------\n"
					"tuple\n"
					"    The point in the internal coordinate system.\n";
	
	const char* doc_inverse_point =	"From internal coordinates, it returns the point in the given coordinate system.\n\n"
					"It returns easting, northing and altitude from in-section x coordinate, in-section\n"
					"y coordinate, cut coordinate\n\n"
					"Parameters\n"
					"----------\n"
					"point : tuple\n"
					"    The three coordinates of the internal point.\n"
					"Returns\n"
					"----------\n"
					"tuple\n"
					"    The point in the given coordinate system\n";
	
	const char* doc_verb =	"Sets the operations as verbose.\n\n"
				"When creating the model, it will advice the user of problems with geometries or matchings.\n\n"
				"Parameters\n"
				"----------\n"
				"verbose : boolean\n"
				"    if geomodelr should be verbose when creating the model.\n";
	
	// Register exception.
	python::class_<GeomodelrException> GeomodelrExceptionClass("GeomodelrException", boost::python::init<std::string>());
	GeomodelrExceptionType = createExceptionClass("GeomodelrException");
	python::register_exception_translator<GeomodelrException>(&translate);
	python::def("faultplane_for_lines", test_faultplane_for_lines);
	python::def("set_verbose", set_verbose, python::args("verbose"), doc_verb);
	python::class_<Section>("Section", python::init<const wstring&, double, const pylist&, 
							const pylist&, const pylist&, 
							const pylist&, const pylist&>())
							.def("info", &Section::info)
							.def("closest", &Section::closest);
	python::class_<Model>("Model", python::init<const pylist&, const pylist&, const pyobject&,
						    const pyobject&, const pylist&>())
						    .def("make_matches", &Model::make_matches)
						    .def("possible_closest", &Model::possible_closest, python::args("point"), doc_possible_closest)
						    .def("model_point", &Model::model_point, python::args("point"), doc_model_point)
						    .def("inverse_point", &Model::inverse_point, python::args("internal_point"), doc_inverse_point)
						    .def("closest", &Model::closest, python::args("point"), doc_closest)
						    .def("closest_topo", &Model::closest_topo, python::args("point"), doc_closest_topo)
						    .def("height", &Model::height, python::args("point"), doc_height)
						    .def("info", &Model::info)
						    .add_property("matches", &Model::get_matches, &Model::set_matches);
}

