# coding=utf

# Geomodelr query tool. Tools for using a geomodelr.com geological model.
# Copyright (C) 2016 Geomodelr, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.



import unittest
import os
import model
import json
import faults
import shared
import utils
import cpp
import numpy as np
from numpy import linalg as la
import cProfile, pstats, StringIO
import sys

class TestGeoModelR(unittest.TestCase):
    def setUp(self):
        pass
        # Profile GeoModelR
        # self.pr = cProfile.Profile()
        # self.pr.enable()
    def tearDown(self):
        pass
        # Print profiling of GeoModelR.
        # self.pr.disable()
        # s = StringIO.StringIO()
        # sortby = 'cumulative'
        # ps = pstats.Stats(self.pr, stream=s).sort_stats(sortby)
        # ps.print_stats()
        # print >> sys.stderr, s.getvalue()

    # Tests function fault plane for lines.
    def test_faultplane_for_lines(self):
        la = [(0, 0, 0), (1,1,0), (2,2,0), (3,3,0)]
        lb = [(0, 3, 1), (1,2,1), (2,1,1), (3,0,1)]
        # This is a weird testcase that should work. Star in two directions.
        
        self.assertEqual(cpp.faultplane_for_lines(la, lb), [(0, 1, 4), (1, 4, 5), (1, 2, 5), (2, 3, 5), (3, 5, 6), (3, 6, 7)])
        
        # This is a weird testcase that should not work.
        # Open circles in different directions.
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1,-2,1), (-2,-1,1), (-2,1,1),  (-1, 2,1), (1, 2,1), (2, 1,1), (2,-1,1), (1, -2, 1)]
        with self.assertRaises(cpp.GeomodelrException):
            cpp.faultplane_for_lines(la, lb)
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        
        # Open circles in different directions but one on top of the other.
        with self.assertRaises(cpp.GeomodelrException):
            cpp.faultplane_for_lines(la, lb)
        
        # One starts when the other finishes.
        la = [(-1,-2,0), (-2,-1,0), (-2,1,0),  (-1, 2,0), (1, 2,0), (2, 1,0), (2,-1,0), (1, -2, 0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        self.assertEqual(cpp.faultplane_for_lines(la, lb), [(0,1,  8), (1,  8,  9), (1,  2,  9), 
                                                                (2, 9,10), (2,  3, 10), (3, 10, 11), 
                                                                (3,11,12), (3,  4, 12), (4, 12, 13), 
                                                                (4, 5,13), (5, 13, 14), (5,  6, 14), 
                                                                (6,14,15), (6,  7, 15)])
        # Case from Aburra's Valley
        la = [(827675.59327569162, 1165500.0, 1628.0922045583072), 
              (827765.51647690905, 1165500.0, 1378.975507981105), 
              (827863.2279581792,  1165500.0, 1033.728274159945), 
              (827921.85484694131, 1165500.0, 740.5938303495259), 
              (828054.66909774847, 1165500.0, 402.1243423849806), 
              (828144.1721766293, 1165500.0, 132.04487628862262), 
              (828230.53500000108, 1165500.0, 0), 
              (842128.59300000034, 1165500.0, -30000)]
        lb = [(826965.29925951734, 1169500.0, 1712.9395692311227), 
              (826913.36802500673, 1169500.0, 1563.2193986860414), 
              (826903.19379854389, 1169500.0, 1319.0379635840654), 
              (826964.23915732, 1169500.0, 983.2884903208663), 
              (827025.28451609518, 1169500.0, 671.2788788030545), 
              (827055.80719548278, 1169500.0, 399.96617313536507), 
              (827160.72890587803, 1169500.0, 132.04487628862262), 
              (827315.46193332784, 1169500.0, -284.6906919181347), 
              (840876.00000000093, 1169500.0, -30000)]
        
        self.assertEqual(cpp.faultplane_for_lines(la, lb), [(0, 8, 9), (0, 9, 10), (0, 1, 10), 
                                                               (1, 10, 11), (1, 11, 12), (1, 12, 13), 
                                                               (1, 2, 13), (2, 3, 13), (3, 13, 14), 
                                                               (3, 4, 14), (4, 5, 14), (5, 14, 15), 
                                                               (5, 6, 15), (6, 15, 16), (6, 7, 16)])
        
        # Test that you can create cross sections and that bugs throw something in python (and not segfault).
    def test_sections(self):
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]]]
        units = ["unit1"]
        lines = []
        lnames = []
        
        # Try the simplest one.
        section = cpp.Section("A-A", 0, points, polygons, units, lines, lnames)
        
        # Try a bad index.
        polygons = [[[0, 1, 2, 3, 6]]]
        with self.assertRaises(IndexError):
            section = cpp.Section("B-B", 1, points, polygons, units, lines, lnames)
        
        # Try a bad number.
        points = [0, [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]]]
        with self.assertRaises(TypeError):
            section = cpp.Section("C-C", 2, points, polygons, units, lines, lnames)
        
        # Try a bad number of polygons.
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]], []]
        with self.assertRaises(IndexError):
            section = cpp.Section("D-D", 3, points, polygons, units, lines, lnames)
        
        # Try a bad polygon. The crossing lines goes from [1,1] to [0,0].
        points = [[0, 0], [1, 0], [0, 1], [1, 1]]
        polygons = [[[0, 1, 2, 3]]]
        section = cpp.Section("E-E", 4, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['polygons'], 1)
        
        # Test polygons with holes.
        points = [[0, 0], [1, 0], [0, 1], [1, 1]]
        polygons = [[[0, 1, 2, 3]]]
        section = cpp.Section("F-F", 5, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['polygons'], 1)
        
        # Test lines get created.
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]]
        lines = [[0, 2, 4], [1, 3]]
        lnames = ['fault1', 'fault2']
        section = cpp.Section("G-G", 6, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['lines'], 2)
        
        # Test polygons with holes.
        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]]]
        units = ['unit1', 'unit2']
        section = cpp.Section("H-H", 7, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['polygons'], 2)
    
    # Test the section closest function, without faults.
    def test_section_closest(self):
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], 
                  [0.75, 0.25], [0.75, 0.75], [0.25, 0.75], [2, 0], [2, 2], [0, 2]]

        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]], [[2, 1, 9]], [[4, 3, 2, 10, 11]]]
        units = ['unit1', 'unit2', 'unit3', 'unit4']
        section = cpp.Section("I-I", 8, points, polygons, units, [], [])
        self.assertEqual(section.info()['polygons'], 4)
        self.assertEqual(section.closest([0.5, 0.5]), (1, 'unit2'))
        self.assertEqual(section.closest([-0.5, -0.5]), (0, 'unit1'))
        self.assertEqual(section.closest([1.5, -0.5]), (2, 'unit3'))
        self.assertEqual(section.closest([1.5, 0.75]), (0, 'unit1'))
        self.assertEqual(section.closest([1.5, 0.25]), (2, 'unit3'))
        self.assertEqual(section.closest([-0.001, 1.001]), (3, 'unit4'))
         
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], 
                  [0.75, 0.25], [0.75, 0.75], [0.25, 0.75], [2, 0], [2, 2], [0, 2]]

        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]], [[2, 1, 9]], [[4, 3, 2, 10, 11]]]
        units = ['NONE', 'unit2', 'unit3', 'unit4']
        section = cpp.Section("J-J", 9, points, polygons, units, [], [])
        self.assertEqual(section.info()['polygons'], 3)
        self.assertEqual(section.closest([0.5, 0.5]), (0, 'unit2'))
        self.assertEqual(section.closest([-0.5, -0.5]), (0, 'unit2'))
        self.assertEqual(section.closest([1.5, -0.5]), (1, 'unit3'))
        self.assertEqual(section.closest([1.5, 0.75]), (1, 'unit3'))
        self.assertEqual(section.closest([1.5, 0.25]), (1, 'unit3'))
        self.assertEqual(section.closest([-0.001, 1.001]), (2, 'unit4'))
        
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0.9], [0.9, 0.9], [0.9, 0.1], [0, 0.1],
                  [0.25, 0.25], [0.9, 0.25], [0.9, 0.75], [0.25, 0.75] ]
        polygons = [[[0, 1, 2, 3, 4, 5, 6, 7]], [[8, 9, 10, 11]]]
        units = ["unit1", "unit2"]
        section = cpp.Section("K-K", 10, points, polygons, units, [], [])
        self.assertEqual(section.closest([0, 0.5]), (1, 'unit2'))
        self.assertEqual(section.closest([0, 0.3]), (0, 'unit1'))
    
    # Test the model matches correctly, test it launches GeomodelrException and that inherits from Exception.
    def test_model_matching(self):
        points_1   = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], 
                      [0.75, 0.25], [0.75, 0.75], [0.25, 0.75], [2, 0], [2, 2], [0, 2]]

        polygons_1 = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]], [[2, 1, 9]], [[4, 3, 2, 10, 11]]]
        units_1 = ["unit1", "unit2", "unit3", "unit4"]
        points_2   = [[0, 0], [1, 0], [1, 2], [0, 2], 
                      [0.25, 1.25], [0.75, 1.25], [0.75, 1.75], [0.25, 1.75], 
                      [1, 1], [2, 0], [2, 1], [2, 2]]
        polygons_2 = [[[0, 1, 2, 3], [7, 6, 5, 4]], [[4, 5, 6, 7]], [[1, 10, 8]], [[1, 9, 10]], [[8, 10, 11, 2]]]
        units_2 = ["unit1", "unit2", "unit3", "unit1", "unit4"]
        
        model = cpp.Model([0,0,2,2],[1, 0], [0, 1], {}, {}, [("A-A", 11, points_1, polygons_1, units_1, [], []), ("B-B", 12, points_2, polygons_2, units_2, [], [])])
        
        model.make_matches()
        self.assertEqual(model.matches, [((u'A-A', u'B-B'), [(0, 0), (0, 3), (2, 2), (3, 4)])])
        
        model.matches = [((u'A-A', u'B-B'), [(0, 0), (2, 2), (3, 4)])]
        self.assertEqual(model.matches, [((u'A-A', u'B-B'), [(0, 0), (2, 2), (3, 4)])])
        
        with self.assertRaises(cpp.GeomodelrException):
            model.matches = [((u'B-B', u'C-C'), [(0, 0), (2, 2), (3, 4)])]
        with self.assertRaises(Exception):
            model.matches = [((u'C-C', u'D-D'), [(0, 0), (2, 2), (3, 4)])]
    
    def test_closest_single(self):
        # Evaluate a model that has a hole in the middle.
        points_1 = [[0, 0], [3, 0], [3, 3], [0, 3], [1, 1], [2, 1], [2, 2], [1, 2]]
        points_2 = [[0, 0], [3, 0], [3, 3], [0, 3]]
        
        polygons_1 = [[[0, 1, 2, 3], [7, 6, 5, 4]], [[4, 5, 6, 7]]]
        polygons_2 = [[[0, 1, 2, 3]]] 
        
        units_1 = ["unit1", "unit2"]
        units_2 = ["unit1"]
        
        model = cpp.Model([0,0,3,3],[0, 0], [1, 0], {}, {}, [("A-A", 1, points_1, polygons_1, units_1, [], []), ("B-B", 2, points_2, polygons_2, units_2, [], [])])
        model.make_matches()
        self.assertEqual(model.model_point([1.5, 0.1, 1.5]), (1.5, 1.5, 0.1))
        pos_cls = model.possible_closest((1.5, 1.1, 1.5))
        self.assertEqual(pos_cls[0][0], 'unit1')
        cls = model.closest((1.5, 1.1, 1.5))
        self.assertEqual(cls[0], 'unit2')


    # Test the possible closest, all the units that can be a match given a line.
    def test_possible_closest(self):
        # Evaluate simple models.
        # Evaluate a model where a middle square changes between unit2 and unit3
        points_1   = [[0, 0], [3, 0], [0, 1], [2, 1], [3, 1], [0, 2], [2, 2], [3, 2]]
        points_2   = [[0, 0], [3, 0], [0, 1], [1, 1], [3, 1], [0, 2], [1, 2], [3, 2]]
        polygons_1 = [[[0, 1, 4, 3, 2]], [[2, 3, 6, 5]], [[3, 4, 7, 6]]]
        units_1 = ["unit1", "unit2", "unit3"]
        
        model = cpp.Model([0,0,3,3],[0, 0], [1, 0], {}, {}, [("A-A", 1, points_1, polygons_1, units_1, [], []), ("B-B", 2, points_2, polygons_1, units_1, [], [])])
        model.make_matches()
        self.assertEqual(model.matches, [((u'A-A', u'B-B'), [(0, 0), (1, 1), (2, 2)])])
        self.assertEqual(model.model_point([1.5, 1.5, 1.5]), (1.5, 1.5, 1.5))
        pos_cls = model.possible_closest([1.5, 1.5, 1.6])
        self.assertEqual(map(lambda v: v[0], pos_cls), ['unit2', 'unit3'])
        self.assertAlmostEqual(pos_cls[0][1], 0.0)
        self.assertAlmostEqual(pos_cls[0][2], 0.5)
        self.assertAlmostEqual(pos_cls[1][1], 0.5)
        self.assertAlmostEqual(pos_cls[1][2], 0.0)
        
        pos_cls = model.possible_closest([1.5, 1.5, 1.4])
        self.assertEqual(map(lambda v: v[0], pos_cls), ['unit2', 'unit3'])
        self.assertAlmostEqual(pos_cls[0][1], 0.0)
        self.assertAlmostEqual(pos_cls[0][2], 0.5)
        self.assertAlmostEqual(pos_cls[1][1], 0.5)
        self.assertAlmostEqual(pos_cls[1][2], 0.0)
        
        pos_cls = model.possible_closest([1.5, 1.5, 1.2])
        self.assertEqual(map(lambda v: v[0], pos_cls), ['unit1', 'unit2', 'unit3'])
        self.assertAlmostEqual(pos_cls[0][1], 0.2)
        self.assertAlmostEqual(pos_cls[0][2], 0.2)
        self.assertAlmostEqual(pos_cls[1][1], 0.0)
        self.assertAlmostEqual(pos_cls[1][2], 0.5)
        self.assertAlmostEqual(pos_cls[2][1], 0.5)
        self.assertAlmostEqual(pos_cls[2][2], 0.0)
        
        cls_1 = model.closest([1.5, 1.1, 1.2])
        cls_2 = model.closest([1.5, 1.5, 1.2])
        cls_3 = model.closest([1.5, 1.9, 1.2])
        
        self.assertEqual(cls_1[0], "unit2")
        self.assertEqual(cls_2[0], "unit1")
        self.assertEqual(cls_3[0], "unit3")
        
        self.assertAlmostEqual(cls_1[1], 0.05)
        self.assertAlmostEqual(cls_2[1], 0.2)
        self.assertAlmostEqual(cls_3[1], 0.05)
        
        # Evaluate all single units in both sides.
        points_1   = [[0, 0], [3, 0], [0, 1], [2, 1], [3, 1], [0, 2], [2, 2], [3, 2]]
        points_2   = [[0, 0], [3, 0], [0, 1], [1, 1], [3, 1], [0, 2], [1, 2], [3, 2]]
        polygons_1 = [[[0, 1, 4, 3, 2]], [[2, 3, 6, 5]], [[3, 4, 7, 6]]]
        units_1 = ["unit1", "unit2", "unit3"]
        units_2 = ["unit4", "unit5", "unit6"]
        
        model = cpp.Model([0,0,3,3],[0, 0], [1, 0], {}, {}, [("A-A", 1, points_1, polygons_1, units_1, [], []), ("B-B", 2, points_2, polygons_1, units_2, [], [])])
        model.make_matches()
        
        cls_1 = model.closest([1.5, 1.1, 1.2])
        cls_2 = model.closest([1.5, 1.499999999999, 1.2])
        cls_3 = model.closest([1.5, 1.9, 1.2])
        cls_4 = model.closest([1.5, 1.1, 1.00000000001])
        cls_5 = model.closest([1.5, 1.499999999999, 0.9])
        
        self.assertEqual(cls_1[0], "unit2")
        self.assertAlmostEqual(cls_1[1], 0.1)
        self.assertEqual(cls_2[0], "unit2")
        self.assertAlmostEqual(cls_2[1], 0.5)
        self.assertEqual(cls_3[0], "unit6")
        self.assertAlmostEqual(cls_3[1], 0.1)
        self.assertEqual(cls_4[0], "unit2")
        self.assertAlmostEqual(cls_3[1], 0.1)
        self.assertEqual(cls_5[0], "unit1")
        self.assertAlmostEqual(cls_3[1], 0.1)
        
        pos_cls = model.possible_closest([1.5, 1.5, 1.2])
        self.assertEqual(pos_cls[0][0], 'unit6')
        self.assertEqual(pos_cls[1][0], 'unit2')
        self.assertAlmostEqual(pos_cls[0][1], 1.0)
        self.assertAlmostEqual(pos_cls[1][1], 0.0)
        self.assertAlmostEqual(pos_cls[0][2], 0.0)
        self.assertAlmostEqual(pos_cls[1][2], 1.0)
        
        # Evaluate a shared unit, and the rest single units.
        units_1 = ["unit1", "unit2", "unit3"]
        units_2 = ["unit1", "unit5", "unit6"]
        
        model = cpp.Model([0,0,3,3],[0, 0], [1, 0], {}, {}, [("A-A", 1, points_1, polygons_1, units_1, [], []), ("B-B", 2, points_2, polygons_1, units_2, [], [])])
        model.make_matches()
        pos_cls = model.possible_closest([1.5, 1.5, 1.25])
        
        self.assertEqual(pos_cls[0][0], 'unit6')
        self.assertEqual(pos_cls[1][0], 'unit1')
        self.assertEqual(pos_cls[2][0], 'unit2')
        self.assertAlmostEqual(pos_cls[0][1], 1.0)
        self.assertAlmostEqual(pos_cls[1][1], 0.25)
        self.assertAlmostEqual(pos_cls[2][1], 0.0)
        
        self.assertAlmostEqual(pos_cls[0][2], 0.0)
        self.assertAlmostEqual(pos_cls[1][2], 0.25)
        self.assertAlmostEqual(pos_cls[2][2], 1.0)
    
    # Test the inverse point, plus other possible_closest and closest tests.
    def test_point_inverse(self):
        model = cpp.Model([0,0,2,2],[0, 0], [1, 0], {}, {}, [("A-A", 1, [], [], [], [], []), ("B-B", 2, [], [], [], [], [])])
        
        self.assertAlmostEqual(model.model_point([1,2,3]), (1,3,2))
        self.assertAlmostEqual(model.model_point([3,2,1]), (3,1,2))
        self.assertAlmostEqual(model.model_point([6,3,0]), (6,0,3))
        self.assertAlmostEqual(model.inverse_point([1,3,2]), (1,2,3))
        self.assertAlmostEqual(model.inverse_point([3,1,2]), (3,2,1))
        self.assertAlmostEqual(model.inverse_point([6,0,3]), (6,3,0))
        
        with self.assertRaises(cpp.GeomodelrException):
            model.closest([1.5,1.5,1.5])
        
        model.make_matches()
        self.assertEqual(model.closest([1.5,1.5,1.5]), ('NONE', float('inf')))
        self.assertEqual(model.possible_closest([1.5, 1.5, 1.5]), [])
        n = la.norm([-1,1])
        
        model = cpp.Model([0,0,2,2],[1, 1], list(np.array([-1, 1])/n), {}, {}, [])
        model.make_matches()
        
        self.assertEqual(model.closest([1,1,1]), ("NONE", float('inf')))
        self.assertEqual(model.possible_closest([-1,0,-1]), [])
        mp = model.model_point([1.5,0.5,0.5])
        self.assertAlmostEqual(mp[0], -n/2)
        self.assertAlmostEqual(mp[1], 0.5)
        self.assertAlmostEqual(mp[2], 0)
        mp = model.model_point([0.5,0.5,0.5])
        self.assertAlmostEqual(mp[0], 0.0)
        self.assertAlmostEqual(mp[1], 0.5)
        self.assertAlmostEqual(mp[2], n/2)
        ip = model.inverse_point([-n/2,0.5,0])
        self.assertAlmostEqual(ip[0], 1.5)
        self.assertAlmostEqual(ip[1], 0.5)
        self.assertAlmostEqual(ip[2], 0.5)
        ip = model.inverse_point([0.0,0.5,n/2])
        self.assertAlmostEqual(ip[0], 0.5)
        self.assertAlmostEqual(ip[1], 0.5)
        self.assertAlmostEqual(ip[2], 0.5)
    
    # Test that unicode can be input into model units.
    def test_unicode_models(self):
        points_1   = [[0, 0], [3, 0], [0, 1], [2, 1], [3, 1], [0, 2], [2, 2], [3, 2]]
        points_2   = [[0, 0], [3, 0], [0, 1], [1, 1], [3, 1], [0, 2], [1, 2], [3, 2]]
        polygons_1 = [[[0, 1, 4, 3, 2]], [[2, 3, 6, 5]], [[3, 4, 7, 6]]]
        units_1 = [u"unít1", u"unót2", u"unét3"]
        
        model = cpp.Model([0,0,3,3],[0, 0], [1, 0], {}, {}, [(u"A-Á", 1, points_1, polygons_1, units_1, [], []), (u"Ñ-Ñ", 2, points_2, polygons_1, units_1, [], [])])
        model.make_matches()
       	self.assertEqual(model.closest([1.5, 1.5, 1.5])[0], u"unót2")
    
    # Test that the model interpolates faults and that it takes the unit of the side it falls.
    def test_faults_model(self):
        points_1 = [[0, 0], [1, 0], [3, 0], [0, 0.5], [5.0/6.0, 0.5], [2.0/3.0, 1], [3, 1], [0.0, 1.5], [0.5, 1.5], [1.0/3.0, 2], [3, 2], [0, 2.5], [1.0/6.0, 2.5], [0, 3], [3,3]]
        points_2 = [[0, 0], [2, 0], [3, 0], [0, 0.5], [2, 0.5], [2, 1], [3, 1], [0, 1.5], [2, 1.5], [2, 2], [3, 2], [0, 2.5], [2, 2.5], [0, 3], [2, 3], [3, 3]]
        pols_1 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 14, 13, 12]], [[11, 12, 13]]]
        pols_2 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 15, 14, 12]], [[11, 12, 14, 13]]]
        units = ["unit1", "unit1", "unit2", "unit2", "unit3", "unit3", "unit4"]
        
        model = cpp.Model([0,0,3,3],[0, 0], [1, 0], {}, {}, [("s1", 1, points_1, pols_1, units, [], []), ("s2", 2, points_2, pols_2, units, [], [])])
        model.make_matches()
        self.assertEqual(model.matches, [((u's1', u's2'), [(0, 0), (1, 0), (1, 1), (2, 2), (3, 2), (3, 3), (4, 4), (5, 4), (5, 5), (6, 6)])])
        
        # Test what happens in the case of no faults.
        self.assertEqual(model.closest([7.0/6.0, 1.5, 2])[0], 'unit3')
        cls = model.closest([7.0/6.0, 1.4, 1.9])
        # Unit3 should dominate the closest.
        self.assertEqual(cls[0], 'unit3')
        self.assertAlmostEqual(cls[1], 0.06)
        
        # Instead, test what happens when there are faults.
        faults_1 = [[1, 4, 5, 8, 9, 12, 13]]
        faults_2 = [[1, 4, 5, 8, 9, 12, 14]]
        fname = ["F1"]
        model = cpp.Model([0,0,2,2],[0, 0], [1, 0], {}, {}, [("s1", 1, points_1, pols_1, units, faults_1, fname), ("s2", 2, points_2, pols_2, units, faults_2, fname)])
        model.make_matches()
        # It should have the unit of the side of the fault.
        cls = model.closest([7.0/6.0, 1.4, 1.9])
        self.assertEqual(cls[0], 'unit2')
        self.assertAlmostEqual(cls[1], 0.0)
        cls = model.closest([7.0/6.0, 1.6, 1.9])
        self.assertEqual(cls[0], 'unit3')
        self.assertAlmostEqual(cls[1], 0.0)
    
    # Test that you can load and test aburra Valley. Test grids and volumes.
    def test_aburra_valley(self):
        this_dir, this_filename = os.path.split(__file__)
        f = open(os.path.join(this_dir, 'test_files', 'aburra_version1.json'))
        m = model.GeologicalModel(json.loads(f.read()))
        m.height([813487.938015, 1164500.0])
        m.height([80000, 1100000])
        bbox = m.geojson['bbox']
        t = {}
        def query_func(p):
            q = m.closest(p)[0]
            if q in t:
                return t[q]
            else:
                l = len(t)
                t[q] = l
                return l
        # First test the simple grid.
        simple_grid = utils.generate_simple_grid(query_func, bbox, 6)
        units = simple_grid['units'].flatten()
        grid = simple_grid['grid'].flatten()
        self.assertEqual(grid.size, 1029)
        
        for i in range(3):
            self.assertEqual(bbox[i], grid[i])
            u = grid.size-3+i
            self.assertEqual(bbox[i+3], grid[u])
        
        usum = {}
        for u in units:
            if u in usum:
                usum[u] += 1
            else:
                usum[u] = 1
        # Now test performance.
        rt = { v: k for k, v in t.iteritems() }
        
        srt = sorted(usum.items(), key = lambda i: rt[i[0]])
        self.assertEqual(map( lambda x: x[1], srt ), [116, 43, 68, 4, 43, 5, 1, 52, 3, 8])
       
        # Then test the fdm refined grid.
        ref_grid = utils.generate_fdm_grid(query_func, bbox, 5, 5)
        
        self.assertEqual(ref_grid['grid'].shape[3], 3)
        self.assertGreaterEqual(ref_grid['grid'].shape[0], ref_grid['units'].shape[0])
        self.assertGreaterEqual(ref_grid['grid'].shape[1], ref_grid['units'].shape[1])
        self.assertGreaterEqual(ref_grid['grid'].shape[2], ref_grid['units'].shape[2])
        
        grid = utils.generate_octtree_grid(query_func, bbox, 3, 3, 3)
        vols, elems = utils.octtree_volume_calculation(query_func, bbox, 10, 2)

def main(args=None):
    unittest.main()

if __name__ == '__main__':
    main()
