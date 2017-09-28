
geomodelr package
*****************


Submodules
==========


geomodelr.model module
======================

**class geomodelr.model.GeologicalModel(geolojson, delete=True)**

   Bases: ``geomodelr.cpp.Model``

   Interface to query a Geological model from Geomodelr.com. The
   models in Geomodelr.com are saved in Geological JSON. A Geological
   JSON is a set of GeoJSON FeatureCollections  with a transformation.
   Go to Geomodelr.com, create a new model and use it with this  tool.

   **closest((Model)arg1, (object)point) -> tuple :**

      Given a point, it finds the geological unit that’s defined as
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
         boost::python::tuple closest(ModelPython
         {lvalue},boost::python::api::object)

   **closest_topo((Model)arg1, (object)point) -> tuple :**

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
         distance to that unit or AIR if it’s above the topography.

      C++ signature :
         boost::python::tuple closest_topo(ModelPython
         {lvalue},boost::python::api::object)

   **height((Model)arg1, (object)point) -> float :**

      Returns: the height at the given point at the topography.

      It returns the height at the point stored in the topography. In
      case the point it’s outside the bounds of the model, it returns
      the height of the closest point inside.

      Args:
         (tuple)point: The two coordinates (easting, northing) of the
         point in the given coordinate system.

      Returns:
         (real) The height as stored in the topography.

      C++ signature :
         double height(ModelPython
         {lvalue},boost::python::api::object)

   **info((Model)arg1) -> dict :**

      C++ signature :
         boost::python::dict info(ModelPython {lvalue})

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

   **inverse_point((Model)arg1, (object)internal_point) -> tuple :**

      From internal coordinates, it returns the point in the given
      coordinate system.

      It returns easting, northing and altitude from in-section x
      coordinate, in-section y coordinate, cut coordinate

      Args:
         (tuple) point: The three coordinates of the internal point.

      Returns:
         (tuple) The point in the given coordinate system

      C++ signature :
         boost::python::tuple inverse_point(ModelPython
         {lvalue},boost::python::api::object)

   **make_matches()**

      Prepares the model to query by matching polygons and lines. It
      finds which polygons, when projected to the next cross section,
      intersect. After that, it tries to match faults with the same
      name by triangulating them and trying to find a continuous set
      of triangles between the two lines that go from the ends to the
      other side.

   **model_point((Model)arg1, (object)point) -> tuple :**

      Translates the point to internal coordinates

      It returns in-section x coordinate, in-section y coordinate, cut
      coordinate

      Args:
         (tuple) point: The three coordinates (esting, norting,
         altitute a.s.l) of the point in the given coordinate system.

      Returns:
         (tuple) The point in the internal coordinate system.

      C++ signature :
         boost::python::tuple model_point(ModelPython
         {lvalue},boost::python::api::object)

   **possible_closest((Model)arg1, (object)point) -> list :**

      Given a point, it finds all the possible geological units in the
      given line for the adjacent cross sections.

      It can be used to query a grid aligned with the model faster,
      for purposes of generating a triangulated mesh or a grid.

      Args:
         (tuple) point: The three coordinates of the point in the
         given coordinate system.

      Returns:
         (list(tuple(), …)): a list with all the possible units, each
         unit with the distance to the previous and next cross
         sections.

      C++ signature :
         boost::python::list possible_closest(ModelPython
         {lvalue},boost::python::api::object)

   **print_information(verbose=False)**

      Prints the information of the geological model just loaded.

      Prints the version, coordinate system and valid coordinates
      that the geological model takes.

      Args:
         (boolean) verbose: You can print more information with
         verbose=True.

   **signed_distance((Model)arg1, (unicode)unit, (object)point) ->
   float :**

      Given unit U and a point P, it finds the geomodelr distance to U
      minus the geomodelr distance to the closest unit different to U

      It returns a signed distance that’s zero at the boundary of the
      unit, negative inside the unit and possitive outside the unit

      Args:
         (string) unit: The unit to measure the signed distance to
         (tuple) point: The three coordinates (esting, norting,
         altitute a.s.l) of the point in the given coordinate system.

      Returns:
         (double) The signed distance from the unit to the point.

      C++ signature :
         double signed_distance(ModelPython
         {lvalue},std::__cxx11::basic_string<wchar_t,
         std::char_traits<wchar_t>, std::allocator<wchar_t>
         >,boost::python::api::object)

   **signed_distance_bounded((Model)arg1, (unicode)unit,
   (object)point) -> float :**

      Given unit U and a point P, it finds the geomodelr distance to U
      minus the geomodelr distance to the closest unit different to U

      It returns a signed distance that’s zero at the boundary of the
      unit, negative inside the unit and possitive outside the unit

      unlike signed_distance, when the point is outside the bounds of
      the model, or above the topography, it returns a positive number
      (outside)

      Args:
         (string) unit: The unit to measure the signed distance to
         (tuple) point: The three coordinates (esting, norting,
         altitute a.s.l) of the point in the given coordinate system.

      Returns:
         (double) The signed distance from the unit to the point.

      C++ signature :
         double signed_distance_bounded(ModelPython
         {lvalue},std::__cxx11::basic_string<wchar_t,
         std::char_traits<wchar_t>, std::allocator<wchar_t>
         >,boost::python::api::object)

   **signed_distance_unbounded((Model)arg1, (unicode)unit,
   (object)point) -> float :**

      Given unit U and a point P, it finds the geomodelr distance to U
      minus the geomodelr distance to the closest unit different to U

      It returns a signed distance that’s zero at the boundary of the
      unit, negative inside the unit and possitive outside the unit

      unlike signed_distance unbounded, it just returns a positive
      number when the point is above the topography. It does not
      always produce solids

      Args:
         (string) unit: The unit to measure the signed distance to
         (tuple) point: The three coordinates (esting, norting,
         altitute a.s.l) of the point in the given coordinate system.

      Returns:
         (double) The signed distance from the unit to the point.

      C++ signature :
         double signed_distance_unbounded(ModelPython
         {lvalue},std::__cxx11::basic_string<wchar_t,
         std::char_traits<wchar_t>, std::allocator<wchar_t>
         >,boost::python::api::object)

   **validate()**

      Validates that the Geological JSON has correct information.

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


geomodelr.cpp module
====================

**geomodelr.cpp.faultplane_for_lines((list)arg1, (list)arg2) -> list
:**

   C++ signature :
      boost::python::list
      faultplane_for_lines(boost::python::list,boost::python::list)

**geomodelr.cpp.set_verbose((bool)verbose) -> None :**

   Sets the operations as verbose.

   When creating the model, it will advice the user of problems with
   geometries or matchings.

   Args:
      (boolean) verbose: if geomodelr should be verbose when creating
      the model.

   C++ signature :
      void set_verbose(bool)


geomodelr.utils module
======================

**geomodelr.utils.generate_fdm_grid(query_func, bbox, grid_divisions,
max_refinements)**

   Generates a grid of points with a FDM like refinment method. It
   first generates a simple grid. then it checks if a cell needs
   refinement. If it does, it marks it as a cell to refine. Then it
   goes through every axis, creating planes where the cell needs
   refinements, plus marking the cells as not needing refinement.

   Args:
      (function) query_func: a function of the geological model that
      returns a unit.

      (list) bbox: the bounding box to search in.

      (int) grid_divisions: the number of points for the grid.

      (int) max_refinements: the number of refinements for this FDM
      scheme.

**geomodelr.utils.generate_octtree_grid(query_func, bbox,
grid_divisions, fdm_refine, oct_refine)**

   Generates an octree grid, starting with an FDM refined grid. The
   octtree grid divides each cell in 8 looking at the differences of
   material until reaching the number of refinements.

   Args:
      (function) query_func: a function of the geological model that
      returns a unit.

      (list) bbox: the bounding box to search in.

      (int) grid_divisions: the number of points for the grid.

      (int) fdm_refine: the number of refinements for the fdm scheme.

      (int) oct_refine: the number of refinements for the octree
      scheme

**geomodelr.utils.generate_simple_grid(query_func, bbox,
grid_divisions)**

   Returns a uniform grid of sizes grid_divisions x grid_divisions x
   grid_divisions  that covers the given bbox evaluated with the query
   function.

   Args:
      (function) query_func: a function of the geological model that
      returns a unit.

      (list) bbox: the bounding box to search in.

      (int) grid_divisions: the number of points for the grid.

**geomodelr.utils.octtree_volume_calculation(query_func, bbox,
grid_divisions, oct_refine)**

   An example of how to get the volumes of all units.

   Args:
      (function) query_func: a function of the geological model that
      returns a unit.

      (list) bbox: the bounding box to search in.

      (int) grid_divisions: the number of points for the grid.

      (int) oct_refine: the number of refinements for the octree
      scheme


geomodelr.isosurfaces module
============================

**geomodelr.isosurfaces.calculate_isosurface(model, unit,
grid_divisions, bounded=True, filter_by_normal=False,
normal_upwards=True)**

   Calculates an isosurface of a unit. It uses a signed distance and
   an isosurface algorithm present in skimage.measure.

   Args:
      (geomodelr.model.GeologicalModel) model: The geomodelr
      geological model.

      (unicode) unit: The unit to calculate the isosurface for.

      (int) grid_divisions: The number of divisions for all the axes.

      (bool) bounded: calculates the surface using the bounding box.
      This will result in a solid.

      (bool) filter_by_normal: filters by the normal of each triangle,
      depending on the normal_upwards argument.

      (bool) normal_upwards: if filter_by_normal is True, filters the
      triangles depending on its normal. It returns only the triangles
      pointing upwards if it’s True, otherwise it returns the
      triangles pointing downwards.

   Returns:
      (list) vertices: The list of vertices.

      (list) triangles: The list of triangles indexes to vertexes.

**geomodelr.isosurfaces.calculate_isovalues(model, unit,
grid_divisions, bbox, bounded=True)**

   Calculates a grid of isovalues in a bounding box to be used by
   method  skimage.measure to create a triangulation of the unit.

   Args:
      (geomodelr.model.GeologicalModel) geolojson: The Geological
      model created form a Geological JSON object.

      (int) grid_divisions: The divisions of the grid in the X, Y and
      Z directions.

      (list) bbox: the values of [minx, miny, minz, maxx, maxy, maxz].

      (bool) bounded: whether the object will be a solid, (the bounded
      signed distance will be used), or it will be an open surface.

**geomodelr.isosurfaces.calculate_normals(vertices, simplices)**

   Calculates the normals for all the simplices, (skimage returns
   normals per vertex).

   Args:
      (numpy.array) vertices: The vertices returned by the marching
      cubes algorithm.

      (numpy.array) simplices: The simplices (triangles) returned by
      the marching cubes algorithm.

**geomodelr.isosurfaces.check_bbox_surface(sd)**

   Checks if the bounding box of the object modeled is very small.
   When a geological unit covers a very small part of the model, it
   needs to be refined. The new bbox of the unit is returned to check
   that case.

   Args:
      (numpy.array) sd: The grid of signed distances.

**geomodelr.isosurfaces.plot_unit(model, unit, grid_divisions,
bounded=True, filter_by_normal=False, normal_upwards=False)**

   Plots a unit previously modeled with calculate_isosurface.

   Args:
      (geomodelr.model.GeologicalModel) model: The geomodelr
      geological model.

      (unicode) unit: The unit to calculate the isosurface for.

      (int) grid_divisions: The number of divisions for all the axes.

      (bool) bounded: calculates the surface using the bounding box.
      This will result in a solid.

      (bool) filter_by_normal: filters by the normal of each triangle,
      depending on the normal_upwards argument.

      (bool) normal_upwards: if filter_by_normal is True, filters the
      triangles depending on its normal. It returns only

      the triangles pointing upwards if it’s True, otherwise it
      returns the triangles pointing downwards.

**geomodelr.isosurfaces.save_unit(name, model, unit, grid_divisions,
bounded=True, filter_by_normal=False, normal_upwards=False)**

   Saves the wireframe of a geological unit to the disk. It uses a
   marching cubes and a signed distance from the model.

   Args:
      (str) name: the name to save the STL file.

      (geomodelr.model.GeologicalModel) model: the model to be
      queried.

      (unicode) unit: the unit to be meshed.

      (int) grid_divisions: the number of divisions of the grid to
      mesh the object.

      (bool) bounded: whether this surface is bounded by the bbox or
      only by the topography.

      (bool) filter_by_normal: whether to filter this mesh by normal
      to the surface. Useful  if you want to see the top or bottom of
      your formation.

      (bool) normal_upwards: if filter_by_normal is True, whether you
      want the triangles that look up or the triangles that look down.

**geomodelr.isosurfaces.stl_mesh(vertices, simplices)**

   Creates a numpy-stl mesh from a sets of vertices and triangles.

   Args:
      (list) vertices: vertices of the mesh.

      (list) simplices: triangles of the mesh.

   Returns:
      (stl.mesh.Mesh): a numpy-stl mesh.

**geomodelr.isosurfaces.triangulate_unit(model, unit, grid_divisions,
bounded=True, filter_by_normal=False, normal_upwards=False)**

   Triangulates a geological unit and returns it for further
   processing, (or saving it to the database). It uses a marching
   cubes and a signed distance from the model.

   Args:
      (geomodelr.model.GeologicalModel) model: the model to be
      queried.

      (unicode) unit: the unit to be meshed.

      (int) grid_divisions: the number of divisions of the grid to
      mesh the object.

      (bool) bounded: whether this surface is bounded by the bbox or
      only by the topography.

      (bool) filter_by_normal: whether to filter this mesh by normal
      to the surface. Useful  if you want to see the top or bottom of
      your formation.

      (bool) normal_upwards: if filter_by_normal is True, whether you
      want the triangles that look up or the triangles that look down.

   Returns:
      (dict): the triangulated unit with a few useful properties.


Module contents
===============
