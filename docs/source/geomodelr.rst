
geomodelr package
*****************


Submodules
==========


geomodelr.model module
======================

**geomodelr.model.model_from_file(filename)**

   Entry point for the API. It creates the geological model  from the
   file path. The geological model is a model of  geomodelr.com,
   downloaded as a version.

   Args:
      (str) filename: The path to the Geological JSON file downloaded
      from  Geomodelr.com.

   Returns:
      (GeologicalModel): The output Geological model to query the
      geological units freely.

**class geomodelr.model.GeologicalModel(geolojson)**

   Bases: ``geomodelr.cpp.Model``

   Interface to query a Geological model from Geomodelr.com. The
   models in Geomodelr.com are saved in Geological JSON. A Geological
   JSON is a set of GeoJSON FeatureCollections  with a transformation.
   Go to Geomodelr.com, create a new model and use it with this  tool.

   **closest((object)point) -> tuple :**

      Given a point, it finds the geological unit that's defined as
      the closest to that point.

      The basic definition of the algorithm is that, given a match
      between geological units, the distance from the point to the
      unit is the sum of the in-section distance to the point averaged
      by the distance to the cross section.

      Args:
         (tuple) point: The three coordinates of the point in the
         given coordinate system.

      Returns:
         (tuple): A tuple with the geological unit and the defined
         distance to that unit.

      C++ signature :
         boost::python::tuple closest(Model
         {lvalue},boost::python::api::object)

   **closest_topo((object)point) -> tuple :**

      Same as closest but it returns (AIR, inf) if the point is above
      the topography.

      It first looks if the point is above the topography and returns
      (AIR, inf) in that case. Otherwise it returns the same as
      closest.

      Args:
         (tuple) point: The three coordinates (easting, northing,
         altitude a.s.l) of the point in the given coordinate system.

      Returns:
         (tuple): A tuple with the geological unit and the defined
         distance to that unit or AIR if it's above the topography.

      C++ signature :
         boost::python::tuple closest_topo(Model
         {lvalue},boost::python::api::object)

   **height((object)point) -> float :**

      Returns: the height at the given point at the topography.

      It returns the height at the point stored in the topography. In
      case the point it's outside the bounds of the model, it returns
      the height of the closest point inside.

      Args:
         (tuple)point: The two coordinates (easting, northing) of the
         point in the given coordinate system.

      Returns:
         (real) The height as stored in the topography.

      C++ signature :
         double height(Model {lvalue},boost::python::api::object)

   **info() -> dict :**

      C++ signature :
         boost::python::dict info(Model {lvalue})

   **intersect_plane(plane)**

      Intersects a plane with the faults of the Geological Model.

      Takes a plane represented with its four corners and returns the
      set  of lines that intersect that plane with the faults.

      Args:
         (list) plane: list with the four corners of the plane that we
         want to intersect the fault with.

      Returns:
         (dict): a dictionary with fault names as keys, and lines,
         (list of points) as values. The coordinates go from the
         lower left corner, (0.0, 0.0).

   **intersect_planes(planes)**

      Intersects a set of planes with the faults of the Geological
      Model.

      Takes a set of plane represented with its four corners and
      returns the set  of lines that intersect that plane with the
      faults. The coordinates start from the first plane lower corner,
      and increase by dist(plane[i][0], plane[i][1]) for the next
      plane.

      Args:
         (list) plane: List with planes. Each plane has a list with
         four corners  that we want to intersect the fault with.

      Returns:
         (dict): a dictionary with fault names as keys, and lines,
         (list of points)  as values.

   **inverse_point((object)internal_point) -> tuple :**

      From internal coordinates, it returns the point in the given
      coordinate system.

      It returns easting, northing and altitude from in-section x
      coordinate, in-section y coordinate, cut coordinate

      Args:
         (tuple) point: The three coordinates of the internal point.

      Returns:
         (tuple) The point in the given coordinate system

      C++ signature :
         boost::python::tuple inverse_point(Model
         {lvalue},boost::python::api::object)

   **make_matches()**

      Prepares the model to query by matching polygons and lines.

      It finds which polygons, when projected to the next cross
      section, intersect. After that, it tries to match faults with
      the same name by triangulating them and trying to find a
      continuous set of triangles between the two lines that go from
      the ends to the other side.

   **model_point((object)point) -> tuple :**

      Translates the point to internal coordinates

      It returns in-section x coordinate, in-section y coordinate, cut
      coordinate

      Args:
         (tuple) point: The three coordinates (esting, norting,
         altitute a.s.l) of the point in the given coordinate system.

      Returns:
         (tuple) The point in the internal coordinate system.

      C++ signature :
         boost::python::tuple model_point(Model
         {lvalue},boost::python::api::object)

   **possible_closest((object)point) -> list :**

      Given a point, it finds all the possible geological units in the
      given line for the adjacent cross sections.

      It can be used to query a grid aligned with the model faster,
      for purposes of generating a triangulated mesh or a grid.

      Args:
         (tuple) point: The three coordinates of the point in the
         given coordinate system.

      Returns:
         (list(tuple(), ...)): a list with all the possible units,
         each unit with the distance to the previous and next cross
         sections.

      C++ signature :
         boost::python::list possible_closest(Model
         {lvalue},boost::python::api::object)

   **print_information(verbose=False)**

      Prints the information of the geological model just loaded.

      Prints the version, coordinate system and valid coordinates
      that the geological model takes.

      Parameters:
         (boolean) verbose: You can print more information with
         verbose=True.

   **validate()**

      Validates that the Geological JSON has correct information.

geomodelr.cpp module
====================

**geomodelr.cpp.set_verbose((bool)verbose) -> None :**

   Sets the operations as verbose.

   When creating the model, it will advice the user of problems with
   geometries or matchings.

   Args:
      (boolean) verbose: if geomodelr should be verbose when creating
      the model.

   C++ signature :
      void set_verbose(bool)

