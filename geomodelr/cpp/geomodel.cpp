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

BOOST_PYTHON_MODULE(cpp)
{
	python::class_<Section>("Section", python::init<double, const pylist&, 
							const pylist&, const pylist&, 
							const pylist&, const pylist&>())
							.def("info", &Section::info)
							.def("closest", &Section::closest);
	python::class_<Model>("Model", python::init<const pylist&, const pylist&, 
						    const pylist&>())
						    .def("make_matches", &Model::make_matches)
						    .def("possible_closest", &Model::possible_closest)
						    .def("model_point", &Model::to_model_point)
						    .def("closest", &Model::closest)
						    .add_property("matches", &Model::get_matches, &Model::set_matches);

}

