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

	//PyErr_SetString(PyExc_RuntimeError, e.what());
}

BOOST_PYTHON_MODULE(cpp)
{
	// Register exception.
	python::class_<GeomodelrException> GeomodelrExceptionClass("GeomodelrException", boost::python::init<std::string>());
	GeomodelrExceptionType = createExceptionClass("GeomodelrException");
	python::register_exception_translator<GeomodelrException>(&translate);
	
	python::class_<Section>("Section", python::init<const wstring&, double, const pylist&, 
							const pylist&, const pylist&, 
							const pylist&, const pylist&>())
							.def("info", &Section::info)
							.def("closest", &Section::closest);
	python::class_<Model>("Model", python::init<const pylist&, const pylist&, 
						    const pylist&>())
						    .def("make_matches", &Model::make_matches)
						    .def("possible_closest", &Model::possible_closest)
						    .def("model_point", &Model::model_point)
						    .def("inverse_point", &Model::inverse_point)
						    .def("closest", &Model::closest)
						    .add_property("matches", &Model::get_matches, &Model::set_matches);

}

