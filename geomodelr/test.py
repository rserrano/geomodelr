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
import shared
import cpp
import numpy as np
from numpy import linalg as la
from datetime import datetime
import cProfile, pstats, StringIO
import sys
import matplotlib.pyplot as plt
from random import shuffle
import copy
import math
from shapely.geometry import Polygon, Point

class Tests(unittest.TestCase):
    def setUp(self):
        pass
        #Profile GeoModelR
        # self.pr = cProfile.Profile()
        # self.pr.enable()

    def tearDown(self):
        pass
        #profiling of GeoModelR.
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
        self.assertEqual(cpp.faultplane_for_lines(la, lb),[(0, 1, 4), (1, 4, 5), (1, 2, 5), (2, 5, 6), (2, 3, 6), (3, 6, 7)])
        
        # This is a weird testcase.
        # Open circles in different directions.
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1,-2,1), (-2,-1,1), (-2,1,1),  (-1, 2,1), (1, 2,1), (2, 1,1), (2,-1,1), (1, -2, 1)]
        cpp.faultplane_for_lines(la, lb)
        
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        
        # Open circles in different directions but one on top of the other.
        cpp.faultplane_for_lines(la, lb)
        
        # One starts when the other finishes.
        la = [(-1,-2,0), (-2,-1,0), (-2,1,0),  (-1, 2,0), (1, 2,0), (2, 1,0), (2,-1,0), (1, -2, 0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        
        self.assertEqual(cpp.faultplane_for_lines(la, lb), [(0, 1, 8),   (1, 8, 9),   (1, 2, 9), 
                                                            (2, 9, 10),  (2, 3, 10),  (3, 10, 11), 
                                                            (3, 4, 11),  (4, 11, 12), (4, 5, 12), 
                                                            (5, 12, 13), (5, 6, 13),  (6, 13, 14), 
                                                            (6, 7, 14),  (7, 14, 15)])
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
        
        self.assertEqual(cpp.faultplane_for_lines(la, lb), [(0, 8, 9), (0, 1, 9), (1, 9, 10), (1, 2, 10), 
                                                            (2, 10, 11), (2, 3, 11), (3, 11, 12), (3, 4, 12), 
                                                            (4, 12, 13), (4, 5, 13), (5, 13, 14), (5, 6, 14), 
                                                            (6, 14, 15), (6, 7, 15), (7, 15, 16)])
        
    def test_polygon_fault(self):    
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [0.25, -0.25], [0.75, -0.25] ]
        polygons = [[[0, 1, 2, 3]]]
        units = ["unit1"]
        lines = [ [4, 5] ]
        lnames = ["l1"]
        section = cpp.Section("A-A", 0, (0, 1, 0, 1), points, polygons, units, lines, lnames, [])
        section.params = {'faults': 'cover'}
        # Basic tests.
        self.assertEqual(section.closest(( 0.5, -0.25 ))[0], "NONE")
        self.assertEqual(section.closest(( 0.5, -0.3 ))[0], "NONE" )
        self.assertEqual(section.closest(( 0.5, -0.4 ))[0], "NONE" )
        self.assertEqual(section.closest(( 0.5, -0.49999 ))[0], "NONE")
        self.assertEqual(section.closest(( 0.5, -0.5 ))[0], "unit1")
        self.assertEqual(section.closest(( 0.5, -0.6 ))[0], "unit1")
        self.assertEqual(section.closest(( 0.5, -0.7 ))[0], "unit1")
        self.assertEqual(section.closest(( 0.5, -0.8 ))[0], "unit1")
        self.assertEqual(section.closest(( 0.5, -0.9 ))[0], "unit1")
        self.assertEqual(section.closest(( 0.5, -1.0 ))[0], "unit1")
        
        self.assertAlmostEqual(section.closest(( 0.5, -0.5 ))[1], 0.707106781187)
        self.assertAlmostEqual(section.closest(( 0.5, -0.6 ))[1], 0.737342165746)
        self.assertAlmostEqual(section.closest(( 0.5, -0.7 ))[1], 0.800771233187)
        self.assertAlmostEqual(section.closest(( 0.5, -0.8 ))[1], 0.878766979897)
        self.assertAlmostEqual(section.closest(( 0.5, -0.9 ))[1], 0.964273034574)
        self.assertAlmostEqual(section.closest(( 0.5, -1.0 ))[1], 1.054092553389)
        
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [-0.25, 0.25], [-0.25, 0.75] ]
        section = cpp.Section("A-A", 0, (0, 1, 0, 1), points, polygons, units, lines, lnames, [])
        section.params = {'faults': 'cover'}
        
        self.assertEqual(section.closest(( -0.25, 0.5 ))[0], "NONE")
        self.assertEqual(section.closest(( -0.3, 0.5 ))[0], "NONE" )
        self.assertEqual(section.closest(( -0.4, 0.5 ))[0], "NONE" )
        self.assertEqual(section.closest(( -0.49999, 0.5 ))[0], "NONE")
        self.assertEqual(section.closest(( -0.5, 0.5 ))[0], "unit1")
        self.assertEqual(section.closest(( -0.6, 0.5 ))[0], "unit1")
        self.assertEqual(section.closest(( -0.7, 0.5 ))[0], "unit1")
        self.assertEqual(section.closest(( -0.8, 0.5 ))[0], "unit1")
        self.assertEqual(section.closest(( -0.9, 0.5 ))[0], "unit1")
        self.assertEqual(section.closest(( -1.0, 0.5 ))[0], "unit1")
        
        self.assertAlmostEqual(section.closest(( -0.5, 0.5))[1], 0.707106781187)
        self.assertAlmostEqual(section.closest(( -0.6, 0.5))[1], 0.737342165746)
        self.assertAlmostEqual(section.closest(( -0.7, 0.5))[1], 0.800771233187)
        self.assertAlmostEqual(section.closest(( -0.8, 0.5))[1], 0.878766979897)
        self.assertAlmostEqual(section.closest(( -0.9, 0.5))[1], 0.964273034574)
        self.assertAlmostEqual(section.closest(( -1.0, 0.5))[1], 1.054092553389)
        
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [0, -0.1], [0, 1.1] ]
        section = cpp.Section("A-A", 0, (0, 1, 0, 1), points, polygons, units, lines, lnames, [])
        section.params = {'faults': 'cover'}
        self.assertEqual(section.closest(( -0.25, 0.5 ))[0], "NONE")
        
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0], [0, 1]]
        section = cpp.Section("A-A", 0, (0, 1, 0, 1), points, polygons, units, lines, lnames, [])
        section.params = {'faults': 'cover'}
        self.assertEqual(section.closest(( -0.25, 0.5 ))[0], "NONE")
    
    # Test that you can create cross sections and that bugs throw something in python (and not segfault).
    def test_sections(self):
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]]]
        units = ["unit1"]
        lines = []
        lnames = []
        
        # Try the simplest one.
        section = cpp.Section("A-A", 0, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        
        # Try a bad index.
        polygons = [[[0, 1, 2, 3, 6]]]
        with self.assertRaises(IndexError):
            section = cpp.Section("B-B", 1, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        
        # Try a bad number.
        points = [0, [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]]]
        with self.assertRaises(TypeError):
            section = cpp.Section("C-C", 2, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        
        # Try a bad number of polygons.
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]], []]
        with self.assertRaises(IndexError):
            section = cpp.Section("D-D", 3, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        
        # Try a bad polygon. The crossing lines goes from [1,1] to [0,0].
        points = [[0, 0], [1, 0], [0, 1], [1, 1]]
        polygons = [[[0, 1, 2, 3]]]
        section = cpp.Section("E-E", 4, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        self.assertEqual(section.info()['polygons'], 1)
        
        # Test polygons with holes.
        points = [[0, 0], [1, 0], [0, 1], [1, 1]]
        polygons = [[[0, 1, 2, 3]]]
        section = cpp.Section("F-F", 5, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        self.assertEqual(section.info()['polygons'], 1)
        
        # Test lines get created.
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]]
        lines = [[0, 2, 4], [1, 3]]
        lnames = ['fault1', 'fault2']
        section = cpp.Section("G-G", 6, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        self.assertEqual(section.info()['lines'], 2)
        
        # Test polygons with holes.
        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]]]
        units = ['unit1', 'unit2']
        section = cpp.Section("H-H", 7, (0, 0, 2, 1), points, polygons, units, lines, lnames, [])
        self.assertEqual(section.info()['polygons'], 2)
    
    # Test the section closest function, without faults.
    def test_section_closest(self):
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], 
                  [0.75, 0.25], [0.75, 0.75], [0.25, 0.75], [2, 0], [2, 2], [0, 2]]

        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]], [[2, 1, 9]], [[4, 3, 2, 10, 11]]]
        units = ['unit1', 'unit2', 'unit3', 'unit4']
        section = cpp.Section("I-I", 8, (0, 0, 2, 2), points, polygons, units, [], [], [])
        self.assertEqual(section.info()['polygons'], 4)
        self.assertEqual(section.closest([0.5, 0.5])[0], 'unit2')
        self.assertEqual(section.closest([-0.5, -0.5])[0], 'unit1')
        self.assertEqual(section.closest([1.5, -0.5])[0], 'unit3')
        self.assertEqual(section.closest([1.5, 0.75])[0], 'unit1')
        self.assertEqual(section.closest([1.5, 0.25])[0], 'unit3')
        self.assertEqual(section.closest([-0.001, 1.001])[0], 'unit4')
         
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], 
                  [0.75, 0.25], [0.75, 0.75], [0.25, 0.75], [2, 0], [2, 2], [0, 2]]

        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]], [[2, 1, 9]], [[4, 3, 2, 10, 11]]]
        units = ['NONE', 'unit2', 'unit3', 'unit4']
        section = cpp.Section("J-J", 9, (0, 0, 2, 2), points, polygons, units, [], [], [])
        self.assertEqual(section.info()['polygons'], 3)
        self.assertEqual(section.closest([0.5, 0.5])[0], 'unit2')
        self.assertEqual(section.closest([-0.5, -0.5])[0], 'unit2')
        self.assertEqual(section.closest([1.5, -0.5])[0], 'unit3')
        self.assertEqual(section.closest([1.5, 0.75])[0], 'unit3')
        self.assertEqual(section.closest([1.5, 0.25])[0], 'unit3')
        self.assertEqual(section.closest([-0.001, 1.001])[0], 'unit4')
        
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0.9], [0.9, 0.9], [0.9, 0.1], [0, 0.1],
                  [0.25, 0.25], [0.9, 0.25], [0.9, 0.75], [0.25, 0.75] ]
        polygons = [[[0, 1, 2, 3, 4, 5, 6, 7]], [[8, 9, 10, 11]]]
        units = ["unit1", "unit2"]
        section = cpp.Section("K-K", 10, (0, 0, 2, 2), points, polygons, units, [], [], [])
        self.assertEqual(section.closest([0, 0.5])[0], 'unit2')
        self.assertEqual(section.closest([0, 0.3])[0], 'unit1')
    
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
        
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[1, 0], [0, 1], [], {}, [["A-A", 11, points_1, polygons_1, units_1, [], [], []], ["B-B", 12, points_2, polygons_2, units_2, [], [], []]], {}, {})
        
        model.make_matches()
        self.assertEqual(model.matches, [((u'A-A', u'B-B'), [(0, 0), (3, 4)])])
        
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], [], {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_2, units_2, [], [], []]], {}, {})
        model.make_matches()
        self.assertEqual(model.model_point([1.5, 0.1, 1.5]), (1.5, 1.5, 0.1))
        cls = model.closest((1.5, 1.1, 1.5))
        self.assertEqual(cls[0], 'unit2')
        self.assertAlmostEqual(cls[1], 0.1)
        ip = model.inverse_point((1.5, 1.1, 1.5))
        clsa = model.closest_aligned(ip)
        self.assertEqual(clsa[0], 'unit2')
        self.assertAlmostEqual(clsa[1], 0.1)
        

    # Test the possible closest, all the units that can be a match given a line.
    def test_possible_closest(self):
        # Evaluate simple models.
        # Evaluate a model where a middle square changes between unit2 and unit3
        points_1   = [[0, 0], [3, 0], [0, 1], [2, 1], [3, 1], [0, 2], [2, 2], [3, 2]]
        points_2   = [[0, 0], [3, 0], [0, 1], [1, 1], [3, 1], [0, 2], [1, 2], [3, 2]]
        polygons_1 = [[[0, 1, 4, 3, 2]], [[2, 3, 6, 5]], [[3, 4, 7, 6]]]
        units_1 = ["unit1", "unit2", "unit3"]
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], [], {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_1, [], [], []]], {}, {})
        model.make_matches()
        self.assertEqual(model.matches, [((u'A-A', u'B-B'), [(0, 0), (1, 1), (2, 2)])])
        self.assertEqual(model.model_point([1.5, 1.5, 1.5]), (1.5, 1.5, 1.5))
        
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], [], {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {}, {})
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
        
        # Evaluate a shared unit, and the rest single units.
        units_1 = ["unit1", "unit2", "unit3"]
        units_2 = ["unit1", "unit5", "unit6"]
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], [], {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {}, {})
        model.make_matches()
        
    def test_closest_fault(self):

        xy = np.array([0,1,2])
        points=[]
        for i in range(3):
            for j in range(3):
                points.append([xy[j],xy[i]])

        polygons_1 = [[[3,4,7,6]],[[0,1,4,3]],[[1,2,5,8,7,4]]]
        polygons_2 = [[[3,4,7,6]],[[4,5,8,7]],[[0,1,4,3]],[[1,2,5,4]]]

        units_1 = ['unit1', 'unit2', 'unit3']
        units_2 = ['unit1', 'unit2', 'unit2', 'unit3']
        lines = [[1,4,7]]
        #lines=[]
        lnames = ['fault1']
        #lanmes = []
        model = cpp.Model([0,0,0,2,1,2],[0,0,0,2,1,2],[0,0],[1,0], [], {},
            [["A-A", 0, points, polygons_1, units_1, lines, lnames, []],
            ["B-B", 1., points, polygons_2, units_2, lines, lnames, []]], { "fault1": "FAULT" }, {'faults': 'cover'})
        model.make_matches()
        
        y_line = lambda u,v: np.sqrt((u-1)**2 + (v-1)**2)/( v-1 + np.sqrt((u-1)**2 + (v-1)**2))
        
        n = 10000
        epsilon = 1e-5;
        for k in range(n):
            p_xz = np.random.rand(2)*(1-2*epsilon) + epsilon +1.0
            py = np.random.rand()*(1-2*epsilon) + epsilon
            # y_val = y_line(p_xz[0],p_xz[1])
            y_val = 1.0/p_xz[1]
            unit = model.closest([p_xz[0],py,p_xz[1]])[0]
            if py>y_val:
                self.assertEqual(unit,u'unit2')
            else:
                self.assertEqual(unit,u'unit3')

        # crosses fault

        X = [0, 0.5, 1.5, 2.]
        Y = [0., 1., 2.]
        points=[]
        for i in range(3):
            for j in range(4):
                points.append([X[j],Y[i]])

        polygons_1 = [[[4,5,9,8]],[[0,1,5,4]],[[1,3,11,9]]]
        polygons_2 = [[[4,6,10,8]],[[6,7,11,10]],[[0,2,6,4]],[[2,3,7,6]]]

        lines_1 = [[1,5,9]]
        lines_2 = [[2,6,10]]

        model = cpp.Model([0,0,0,2,1,2],[0,0,0,2,1,2],[0,0],[1,0], [], {}, [["A-A", 0, points, polygons_1, units_1, lines_1, lnames, []],
            ["B-B", 1., points, polygons_2, units_2, lines_2, lnames, []]], { "fault1": "FAULT" }, {'faults': 'cover'})
        model.make_matches()

        for k in range(n):
            px = np.random.rand()*(2-2*epsilon) + epsilon
            py = np.random.rand()*(1-2*epsilon) + epsilon
            pz = np.random.rand()*(1-2*epsilon) + epsilon + 1.0
            y1 = px - 0.5
            y2 = 1/pz
            #print px, py, pz
            unit = model.closest([px,py,pz])[0]
            if py>y1:
                self.assertEqual(unit,u'unit1')
            elif py>y2:
                self.assertEqual(unit,u'unit2')
            else:
                self.assertEqual(unit,u'unit3')
    
    def test_point_inverse(self):
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[0, 0], [1, 0], [], {}, [["A-A", 1, [], [], [], [], [], []], ["B-B", 2, [], [], [], [], [], []]], {}, {})
        
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
        n = la.norm([-1,1])
        
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[1, 1], list(np.array([-1, 1])/n), [], {}, [], {}, {})
        model.make_matches()
        
        self.assertEqual(model.closest([1,1,1]), ("NONE", float('inf')))
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], [], {}, [[u"A-Á", 1, points_1, polygons_1, units_1, [], [], []], [u"Ñ-Ñ", 2, points_2, polygons_1, units_1, [], [], []]], {}, {})
        model.make_matches()
       	self.assertEqual(model.closest([1.5, 1.5, 1.5])[0], u"unót2")
    
    # Test that the model interpolates faults and that it takes the unit of the side it falls.
    def test_faults_model(self):
        points_1 = [[0, 0], [1, 0], [3, 0], [0, 0.5], [5.0/6.0, 0.5], [2.0/3.0, 1], [3, 1], [0.0, 1.5], [0.5, 1.5], [1.0/3.0, 2], [3, 2], [0, 2.5], [1.0/6.0, 2.5], [0, 3], [3,3]]
        points_2 = [[0, 0], [2, 0], [3, 0], [0, 0.5], [2, 0.5],       [2, 1],       [3, 1], [0, 1.5],   [2, 1.5],   [2, 2],       [3, 2], [0, 2.5], [2, 2.5],       [0, 3], [2, 3], [3, 3]]
        
        pols_1 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 14, 13, 12]], [[11, 12, 13]]]
        pols_2 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 15, 14, 12]], [[11, 12, 14, 13]]]
        
        units = ["unit1", "unit1", "unit2", "unit2", "unit3", "unit3", "unit4"]
        
        model = cpp.Model([0,0,0,3,3,0],[0,0,0,3,3,0],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, [], [], []], ["s2", 2, points_2, pols_2, units, [], [], []]], {}, {})
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
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[0, 0], [1, 0], [], {}, 
                          [["s1", 1, points_1, pols_1, units, faults_1, fname, []], 
                           ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], 
                          { "F1": "FAULT" }, {})
        
        model.make_matches()
        
        # print [points_1[n] for n in faults_1[0]]
        # print [points_2[n] for n in faults_2[0]]

        # It should have the unit of the side of the fault.
        cls = model.closest([7.0/6.0, 1.4, 1.9])
        self.assertEqual(cls[0], 'unit2')
        self.assertAlmostEqual(cls[1], 0.0)
        
        cls = model.closest([7.0/6.0, 1.6, 1.9])
        self.assertEqual(cls[0], 'unit3')
        self.assertAlmostEqual(cls[1], 0.0)
    
    def test_horizontal_models(self):
        # Evaluate simple models.
        # Evaluate a model where a middle square changes between unit2 and unit3
        points_1   = [[0, 0], [3, 0], [0, 1], [2, 1], [3, 1], [0, 2], [2, 2], [3, 2]]
        points_2   = [[0, 0], [3, 0], [0, 1], [1, 1], [3, 1], [0, 2], [1, 2], [3, 2]]
        polygons_1 = [[[0, 1, 4, 3, 2]], [[2, 3, 6, 5]], [[3, 4, 7, 6]]]
        units_1 = ["unit1", "unit2", "unit3"]
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3], [], {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_1, [], [], []]], {}, {})
        model.make_matches()
        self.assertEqual(model.matches, [((u'A-A', u'B-B'), [(0, 0), (1, 1), (2, 2)])])
        self.assertEqual(model.model_point([1.5, 1.5, 1.5]), (1.5, 1.5, 1.5))
        
        cls_1 = model.closest([1.5, 1.2, 1.1])
        cls_2 = model.closest([1.5, 1.2, 1.5])
        cls_3 = model.closest([1.5, 1.2, 1.9])
        
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[], {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {}, {})
        model.make_matches()
        
        cls_1 = model.closest([1.5, 1.2, 1.1])
        cls_2 = model.closest([1.5, 1.2, 1.499999999999])
        cls_3 = model.closest([1.5, 1.2, 1.9])
        cls_4 = model.closest([1.5, 1.00000000001, 1.1])
        cls_5 = model.closest([1.5, 0.9, 1.499999999999])
        
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
        
        # Evaluate a shared unit, and the rest single units.
        units_1 = ["unit1", "unit2", "unit3"]
        units_2 = ["unit1", "unit5", "unit6"]
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3], [], {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {}, {})
        model.make_matches()
        
    def test_faults_veins_fractures(self):
        this_dir, this_filename = os.path.split(__file__)
        m = model.model_from_file(os.path.join(this_dir, 'test_files', 'talud_faults_veins_fractures.json'))
        self.assertEqual(len(m.faults), 1)
        self.assertEqual(len(m.not_extended_faults), 1)
        self.assertEqual(len(m.fracts), 10)
        self.assertEqual(len(m.not_extended_fracts), 10)
        self.assertEqual(len(m.veins), 2)
        self.assertEqual(len(m.not_extended_veins), 2)

    def test_anchored_lines( self ):
        # cpp.set_verbose( True )
        # Test calculate_section_bbox
        # Direction Y
        self.assertAlmostEqual(cpp.calculate_section_bbox([0, 0, 0, 1, 1, 1], [0.1, 0.5], [0, 1], -0.2), (-0.5, 0.0, 0.5, 1.0))
        # Direction X
        self.assertAlmostEqual(cpp.calculate_section_bbox([0, 0, 0, 1, 1, 1], [0.1, 0.5], [1, 0], -0.2), (-0.1, 0.0, 0.9, 1.0))
        # Diagonal.
        self.assertAlmostEqual(cpp.calculate_section_bbox([0, 0, 0, 1, 1, 1], [0.5, 0.5], [0.7071067811865475, 0.7071067811865475], 0.0), (-0.7071067811865476, 0.0, 0.7071067811865476, 1.0))
        # Diagonal moved a bit.
        self.assertAlmostEqual(cpp.calculate_section_bbox([0, 0, 0, 1, 1, 1], [0.5, 0.5], [0.7071067811865475, 0.7071067811865475], 0.1), (-0.6071067811865476, 0.0, 0.6071067811865476, 1.0))
        # Pitagoric.
        res = cpp.calculate_section_bbox([0, 0, 0, 4, 3, 3], [0, 0], [4./5., -3./5.], 12./5.)
        self.assertAlmostEqual( res[2]-res[0], 5.0 )
        
        # Test extend line.

        line=[(0, 0), (.1, .1), (.5, .5), (1., 1.)] 
        with self.assertRaises( cpp.GeomodelrException ) as cm:
            cpp.extend_line( True, [0, 0, 1, 1], line )
        self.assertEqual( str( cm.exception ), "fault not extended: could not extend line." )
        
        with self.assertRaises( cpp.GeomodelrException ) as cm:
            cpp.extend_line( False, [0, 0, 1, 1], line )
        self.assertEqual( str( cm.exception ), "fault not extended: could not extend line." )
        
        with self.assertRaises( cpp.GeomodelrException ) as cm:
            cpp.extend_line( True, [-1e-20, -1e-20, 1, 1], line )
        self.assertEqual( str( cm.exception ), "fault not extended: the line actually goes to the bounds and does not need modification." )
        
        line=[(0, 0), (-1e-20, -1e-20), (.5, .5), (1., 1.)] 
        with self.assertRaises( cpp.GeomodelrException ) as cm:
            cpp.extend_line( True, [-1e-40, -1e-40, 1, 1], line )
        self.assertEqual( str( cm.exception ), "fault not extended: could not determine direction of line." )
        
        line=[(0, 0), (.3, .3), (.5, .5), (1., 1.)] 
        # Test ok extension.
        self.assertAlmostEqual(cpp.extend_line( False, [0, 0, 1.5, 1.5], line )[-1], (1.5, 1.5))
        self.assertAlmostEqual(cpp.extend_line( True, [-0.5, -0.5, 1, 1], line )[0], (-0.5, -0.5))
        
        line=[(0, .3), (.3, .3), (.5, .5), (1., 1.)] 
        self.assertAlmostEqual(cpp.extend_line( True, [-0.5, -0.5, 1, 1], line )[0], (-0.5, .3))
        line.reverse()
        self.assertAlmostEqual(cpp.extend_line( False, [-0.5, -0.5, 1, 1], line )[-1], (-0.5, .3))
        
        line=[(.3, 0), (.3, .3), (.5, .5), (1., 1.)]
        self.assertAlmostEqual(cpp.extend_line( True, [-0.5, -0.5, 1, 1], line )[0], (.3, -0.5))
        line.reverse()
        self.assertAlmostEqual(cpp.extend_line( False, [-0.5, -0.5, 1, 1], line )[-1], (.3, -0.5))
        
        points_1 = [[0, 0], [1, 0], [3, 0], [0, 0.5], [5.0/6.0, 0.5], [2.0/3.0, 1], [3, 1], [0.0, 1.5], [0.5, 1.5], [1.0/3.0, 2], [3, 2], [0, 2.5], [1.0/6.0, 2.5], [0, 3], [3,3]]
        points_2 = [[0, 0], [2, 0], [3, 0], [0, 0.5], [2, 0.5],       [2, 1],       [3, 1], [0, 1.5],   [2, 1.5],   [2, 2],       [3, 2], [0, 2.5], [2, 2.5],       [0, 3], [2, 3], [3, 3]]
        
        pols_1 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 14, 13, 12]], [[11, 12, 13]]]
        pols_2 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 15, 14, 12]], [[11, 12, 14, 13]]]
        
        units = ["unit1", "unit1", "unit2", "unit2", "unit3", "unit3", "unit4"]
        
        # Instead, test what happens when there are faults.
        faults_1 = [[1, 4, 5, 8, 9, 12, 13]]
        faults_2 = [[1, 4, 5, 8, 9, 12, 14]]
        fname = ["F1"]
        model = cpp.Model([0,0,0,4,4,4],[0,0,0,4,4,4],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False]]], 
                                                                 ["s2", 2, points_2, pols_2, units, faults_2, fname, [[0, False]]]], { "F1": "FAULT" }, {})
        model.make_matches()
        
        ext = model.not_extended_faults["F1"]
        fau = model.faults["F1"]
        
        self.assertEqual(len(ext), 35)
        self.assertEqual(len(fau), 39)
        sext = set(sum(map(list,ext), []))
        sfau = set(sum(map(list,fau), []))
        self.assertEqual(len( sfau - sext ), 2)
        self.assertEqual(len( sext ), 28)
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False]]], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, [[0, False]]]], { "F1": "FAULT" }, {})
        model.make_matches()
        ext2 = model.not_extended_faults["F1"]
        self.assertEqual(len(ext2), 36)
        self.assertEqual(len(model.faults["F1"]), 42)
        sext2 = set(sum(map(list,ext2), []))
        self.assertEqual(len(sext2), 28)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False]]], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" }, {})
        model.make_matches()
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 39)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, []], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" }, {})
        model.make_matches()
        
        f2 = model.not_extended_faults["F1"]
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 36)
        
        faults_1 = [[1, 4, 5, 8, 9, 12, 13]]
        faults_2 = [[1, 4, 5, 8, 9, 12, 14]]
        fname = ["F1"]
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], [], {}, 
                          [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False], [0, True]]], 
                           ["s2", 2, points_2, pols_2, units, faults_2, fname, [[0, False], [0, True]]]], { "F1": "FAULT" }, {})
        model.make_matches()
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 42)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False], [0, True]]], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" }, {})
        model.make_matches()
        
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 39)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, []], ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" }, {})
        model.make_matches()
        
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 36)
        
    def test_aligned(self):
        # Evaluate a model that has a hole in the middle.
        hsq = math.sqrt(2)/2
        points_1 = [[0, 0], [3*hsq, 0], [3*hsq, 3], [0, 3], [hsq, 1], [2*hsq, 1], [2*hsq, 2], [hsq, 2]]
        points_2 = [[0, 0], [3*hsq, 0], [3*hsq, 3], [0, 3]]
        
        polygons_1 = [[[0, 1, 2, 3], [7, 6, 5, 4]], [[4, 5, 6, 7]]]
        polygons_2 = [[[0, 1, 2, 3]]] 
        
        units_1 = ["unit1", "unit2"]
        units_2 = ["unit1"]
        
        model = cpp.Model([0,0,0,3,3,3],[-3*hsq/2,0,-3*hsq/2,3*hsq/2,3,3*hsq/2],[0.75, 2.25], [hsq, -hsq], [], {}, 
                          [["A-A", -3*hsq/2, points_1, polygons_1, units_1, [], [], []],
                           ["B-B",  3*hsq/2, points_2, polygons_2, units_2, [], [], []]], {}, {})
        
        model.make_matches()
        
        # Test model_point, inverse_point.
        self.assertAlmostEqual( model.model_point([1.5, 1.5, 0]), (3*hsq/2, 0, 0) )
        self.assertAlmostEqual( model.model_point([0, 1.5, 0]), (0, 0, -3*hsq/2) )
        self.assertAlmostEqual( model.inverse_point([0, 0, 0]), (0.75, 2.25, 0) )
        self.assertAlmostEqual( model.inverse_point((3*hsq, 0, 0)), (2.25,0.75,0) )
        
        # Test closest, and closest aligned.
        dp5 = math.sqrt( 2 * ( 0.05**2 ) )
        self.assertAlmostEqual(model.closest_aligned(( 3*hsq/2, 1.5,-3*hsq/2+0 )), (u'unit2', 0.0))
        self.assertAlmostEqual(model.closest_aligned(( 3*hsq/2, 1.5,-3*hsq/2+1*dp5)), (u'unit2', 0.07071067811865475))
        self.assertAlmostEqual(model.closest_aligned(( 3*hsq/2, 1.5,-3*hsq/2+2*dp5)), (u'unit2', 0.1414213562373095))
        self.assertAlmostEqual(model.closest_aligned(( 3*hsq/2, 1.5,-3*hsq/2+3*dp5)), (u'unit2', 0.21213203435596423))
        
        # Test signed_distance_aligned
        self.assertAlmostEqual(model.signed_distance_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2+0 )),    0.353553390593)
        self.assertAlmostEqual(model.signed_distance_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2+1*dp5)), 0.271057599455)
        self.assertAlmostEqual(model.signed_distance_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2+2*dp5)), 0.188561808316)
        self.assertAlmostEqual(model.signed_distance_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2+3*dp5)), 0.106066017178)
        self.assertAlmostEqual(model.signed_distance_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2+0 )),    -0.353553390593)
        self.assertAlmostEqual(model.signed_distance_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2+1*dp5)), -0.271057599455)
        self.assertAlmostEqual(model.signed_distance_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2+2*dp5)), -0.188561808316)
        self.assertAlmostEqual(model.signed_distance_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2+3*dp5)), -0.106066017178)
        # Test aligned bboxes unbounded.
        self.assertAlmostEqual(model.abbox, [-3*hsq/2,0,-3*hsq/2,3*hsq/2,3,3*hsq/2])
        
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-0 )), hsq/2)
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-1*dp5)), hsq/2)
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-2*dp5)), hsq/2)
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-3*dp5)), hsq/2)
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-0 )),-hsq/2)
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-1*dp5)),-hsq/2)
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-2*dp5)),-hsq/2)
        self.assertAlmostEqual(model.signed_distance_unbounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-3*dp5)),-hsq/2)
        # Test aligned bboxes bounded.
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-0 )),    hsq/2)
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-3*dp5)), hsq/2)
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-6*dp5)), 6*dp5)
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit1", ( 3*hsq/2, 1.5,-3*hsq/2-9*dp5)), 9*dp5)
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-0 )),    0)
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-1*dp5)), 1*dp5)
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-2*dp5)), 2*dp5)
        self.assertAlmostEqual(model.signed_distance_bounded_aligned("unit2", ( 3*hsq/2, 1.5,-3*hsq/2-3*dp5)), 3*dp5)
    
    def test_faults_matching( self ):
        points_1 = [[0, 0], [1, 0], [3, 0], [0, 0.5], [5.0/6.0, 0.5], [2.0/3.0, 1], [3, 1], [0.0, 1.5], [0.5, 1.5], [1.0/3.0, 2], [3, 2], [0, 2.5], [1.0/6.0, 2.5], [0, 3], [3,3]]
        points_2 = [[0, 0], [2, 0], [3, 0], [0, 0.5], [2, 0.5],       [2, 1],       [3, 1], [0, 1.5],   [2, 1.5],   [2, 2],       [3, 2], [0, 2.5], [2, 2.5],       [0, 3], [2, 3], [3, 3]]
        
        pols_1 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 14, 13, 12]], [[11, 12, 13]]]
        pols_2 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 15, 14, 12]], [[11, 12, 14, 13]]]
        
        units = ["unit1", "unit1", "unit2", "unit2", "unit3", "unit3", "unit4"]
        
        model = cpp.Model([0,0,0,3,3,0],[0,0,0,3,3,0],[0, 0], [1, 0], [], {}, [["s1", 1, points_1, pols_1, units, [], [], []], ["s2", 2, points_2, pols_2, units, [], [], []]], {}, {})
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
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[0, 0], [1, 0], [], {}, 
                          [["s1", 1, points_1, pols_1, units, faults_1, fname, []], 
                           ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], 
                          { "F1": "FAULT" }, {})
        
        model.make_matches()
        self.assertEqual(model.matches, [((u's1', u's2'), [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6)])])
        
        points_1 = [[0, 3], [0, 2], [0, 1], [0, 0], [2, 3], [3, 2], [3.5, 1.5], [4, 1], [3.5, 0.5], [3, 0], [8, 3], [8, 2], [8, 1], [8, 0]]
        points_2 = [[0, 3], [0, 2], [0, 1], [0, 0], [5, 3], [6, 2], [6.5, 1.5], [7, 1], [6.5, 0.5], [6, 0], [8, 3], [8, 2], [8, 1], [8, 0]]
        
        self.assertEqual( len(points_1), 14 )
        self.assertEqual( len(points_2), 14 )
        
        phole = lambda p: [p]
        pols = [[0, 1, 5, 4], [1, 2, 7, 6, 5], [2, 3, 9, 8, 7], [4, 5, 6, 11, 10], [6, 7, 8, 12, 11], [8, 9, 13, 12]]
        pols = map( phole, pols )
        units = ['A', 'B', 'C', 'A', 'B', 'C']
        faults_1 = [[4, 5, 6, 7, 8, 9]]
        faults_2 = [[4, 5, 6, 7, 8, 9]]
        fname = ["F1"]
        
        model = cpp.Model([0, 0, 0, 2, 2, 2], [ 0, 0, 0, 2, 2, 2], [0, 0], [1, 0], [], {},
                          [["s1", 1, points_1, pols, units, faults_1, fname, []], 
                           ["s2", 2, points_2, pols, units, faults_2, fname, []]], 
                           { "F1": "FAULT" }, {})
        
        model.make_matches()
        self.assertEqual(model.matches, [((u's1', u's2'), [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)])])
        points_2 = [[0, 3], [0, 2], [0, 1], [0, 0], [2, 3], [2, 2], [2, 1.5], [2, 1], [2, 0.5], [2, 0], [8, 3], [8, 2], [8, 1], [8, 0]]
        
        model = cpp.Model([0, 0, 0, 2, 2, 2], [ 0, 0, 0, 2, 2, 2], [0, 0], [1, 0], [], {},
                          [["s1", 1, points_1, pols, units, faults_1, fname, []], 
                           ["s2", 2, points_2, pols, units, faults_2, fname, []]], 
                           { "F1": "FAULT" }, {})
        
        model.make_matches()
        self.assertEqual( model.matches, [((u's1', u's2'), [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)])])
        
        points_2 = [[0, 3], [0, 2], [0, 1], [0, 0], [5, 3], [5, 2], [5, 1.5], [5, 1], [5, 0.5], [5, 0], [8, 3], [8, 2], [8, 1], [8, 0]]
        model = cpp.Model([0, 0, 0, 2, 2, 2], [ 0, 0, 0, 2, 2, 2], [0, 0], [1, 0], [], {},
                          [["s1", 1, points_1, pols, units, faults_1, fname, []], 
                           ["s2", 2, points_2, pols, units, faults_2, fname, []]], 
                           { "F1": "FAULT" }, {})
        
        model.make_matches()
        self.assertEqual( model.matches, [((u's1', u's2'), [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)])])
        
        points_2 = [[0, 3], [0, 2], [0, 1], [0, 0], [4, 3], [3, 2], [2.5, 1.5], [2, 1], [2.5, 0.5], [3, 0], [8, 3], [8, 2], [8, 1], [8, 0]]
        model = cpp.Model([0, 0, 0, 2, 2, 2], [ 0, 0, 0, 2, 2, 2], [0, 0], [1, 0], [], {},
                          [["s1", 1, points_1, pols, units, faults_1, fname, []], 
                           ["s2", 2, points_2, pols, units, faults_2, fname, []]], 
                           { "F1": "FAULT" }, {})
        
        model.make_matches()
        self.assertEqual( model.matches, [((u's1', u's2'), [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)])])

    def test_soils(self):
        # Evaluate a model that has a hole in the middle.
        points_map = [[0, 0], [3, 0], [3, 3], [0, 3]]
        points_1 = [[0, 0], [3, 0], [3, 3], [0, 3], [1, 1], [2, 1], [2, 2], [1, 2]]
        points_2 = [[0, 0], [3, 0], [3, 3], [0, 3]]
        
        polygons_map = [[[0, 1, 2, 3]]]
        polygons_1 = [[[0, 1, 2, 3], [7, 6, 5, 4]], [[4, 5, 6, 7]]]
        polygons_2 = [[[0, 1, 2, 3]]]
        
        units_map = ["unit2"]
        units_1 = ["unit1", "unit2"]
        units_2 = ["unit1"]
        
        topography = { 'point': [0, 0], 'dims': [1, 1], 'sample': [3, 3], 'heights': [[3]] }

        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], [points_map, polygons_map, units_map, [], [], []], topography, 
                          [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_2, units_2, [], [], []]], {}, {})
        model.make_matches()
        model.params = {'map': 'soils'}
        model.soil_depths = {'unit2': 0.5}
        
        self.assertEqual(model.height( ( 0.0, 0.0 ) ), 3.0)
        self.assertEqual(model.closest( ( 0.0, 0.0, 3.0 ) )[0], "unit2")
        self.assertEqual(model.closest( ( 0.0, 0.0, 2.5 ) )[0], "unit2")
        self.assertEqual(model.closest( ( 0.0, 0.0, 2.0 ) )[0], "unit1")
        self.assertEqual(model.closest( ( 0.0, 0.0, 1.5 ) )[0], "unit1")
        
        self.assertEqual(model.signed_distance( "unit2", ( 0.0, 0.0, 3.0 ) ), -0.5)
        self.assertEqual(model.signed_distance( "unit2", ( 0.0, 0.0, 2.5 ) ),  0.0)
        self.assertEqual(model.signed_distance( "unit2", ( 0.0, 0.0, 2.0 ) ),  0.5)
        self.assertEqual(model.signed_distance( "unit2", ( 0.0, 0.0, 1.5 ) ),  1.0)
        
        self.assertEqual(model.signed_distance( "unit1", ( 0.0, 0.0, 3.0 ) ),  0.5)
        self.assertEqual(model.signed_distance( "unit1", ( 0.0, 0.0, 2.5 ) ),  0.0)
        self.assertEqual(model.signed_distance( "unit1", ( 0.0, 0.0, 2.0 ) ), -0.5)
        self.assertEqual(model.signed_distance( "unit1", ( 0.0, 0.0, 1.5 ) ), -1.0)
        
        self.assertEqual(model.signed_distance( "unit2", ( 1.5, 0.0, 3.0 ) ),  -0.5 )
        self.assertEqual(model.signed_distance( "unit2", ( 1.5, 0.0, 2.5 ) ),   0.0 )
        self.assertEqual(model.signed_distance( "unit2", ( 1.5, 0.0, 2.25 ) ),  0.25)
        self.assertEqual(model.signed_distance( "unit2", ( 1.5, 0.0, 2.0 ) ),   0.0 )
        self.assertEqual(model.signed_distance( "unit2", ( 1.5, 0.0, 1.5 ) ),  -0.5 )
        
        self.assertEqual(model.signed_distance( "unit1", ( 1.5, 0.0, 3.0 ) ),   0.5 )
        self.assertEqual(model.signed_distance( "unit1", ( 1.5, 0.0, 2.5 ) ),   0.0 )
        self.assertEqual(model.signed_distance( "unit1", ( 1.5, 0.0, 2.25 ) ), -0.25)
        self.assertEqual(model.signed_distance( "unit1", ( 1.5, 0.0, 2.0 ) ),   0.0 )
        self.assertEqual(model.signed_distance( "unit1", ( 1.5, 0.0, 1.5 ) ),   0.5 )

    def test_v1(self):
        this_dir, this_filename = os.path.split(__file__)
        fn = os.path.join(this_dir, 'test_files', 'v1.json')
        geo_model = model.model_from_file(fn)
        geo_model.params = { "faults": "cover" }
        self.assertEqual(geo_model.matches, [((u'name2', u'name1'), [(0, 0), (1, 1), (2, 2)])])
        cls = geo_model.closest((834590.52,1174732.73,1450.80))
        self.assertEqual(cls[0], "U2")
        self.assertAlmostEqual(cls[1], 51.996853528295254)
        cls = geo_model.closest((834588,1174732.73,1450.80))
        self.assertEqual(cls[0], "U1")
        self.assertAlmostEqual(cls[1], 0)

def main(args=None):
    unittest.main()

if __name__ == '__main__':
    main()
