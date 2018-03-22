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
import utils
import isosurfaces
import cpp
import modflow
import numpy as np
from numpy import linalg as la
from datetime import datetime
import cProfile, pstats, StringIO
import sys
import matplotlib.pyplot as plt
from random import shuffle
import copy
import math


class TestGeoModelR(unittest.TestCase):
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
        section = cpp.Section("J-J", 9, (0, 0, 2, 2), points, polygons, units, [], [], [])
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
        section = cpp.Section("K-K", 10, (0, 0, 2, 2), points, polygons, units, [], [], [])
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
        
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[1, 0], [0, 1], {}, {}, [["A-A", 11, points_1, polygons_1, units_1, [], [], []], ["B-B", 12, points_2, polygons_2, units_2, [], [], []]], {})
        
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], {}, {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_2, units_2, [], [], []]], {})
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], {}, {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_1, [], [], []]], {})
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], {}, {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {})
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], {}, {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {})
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
        model = cpp.Model([0,0,0,2,1,2],[0,0,0,2,2,2],[0,0],[1,0], {}, {}, [["A-A", 0, points, polygons_1, units_1, lines, lnames, []],
            ["B-B", 1., points, polygons_2, units_2, lines, lnames, []]], { "fault1": "FAULT" })
        model.make_matches()

        y_line = lambda u,v: np.sqrt((u-1)**2 + (v-1)**2)/( v-1 + np.sqrt((u-1)**2 + (v-1)**2))

        n = 10000

        for k in range(n):
            p_xz = np.random.rand(2)+1.0
            py = np.random.rand()
            y_val = y_line(p_xz[0],p_xz[1])
            unit = model.closest([p_xz[0],py,p_xz[1]])[0]
            # if py>y_val:
            #     self.assertEqual(unit,u'unit2')
            # else:
            #     self.assertEqual(unit,u'unit3')
    
    # Test the inverse point, plus other possible_closest and closest tests.
    def test_point_inverse(self):
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[0, 0], [1, 0], {}, {}, [["A-A", 1, [], [], [], [], [], []], ["B-B", 2, [], [], [], [], [], []]], {})
        
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
        
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[1, 1], list(np.array([-1, 1])/n), {}, {}, [], {})
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 0], [1, 0], {}, {}, [[u"A-Á", 1, points_1, polygons_1, units_1, [], [], []], [u"Ñ-Ñ", 2, points_2, polygons_1, units_1, [], [], []]], {})
        model.make_matches()
       	self.assertEqual(model.closest([1.5, 1.5, 1.5])[0], u"unót2")
    
    # Test that the model interpolates faults and that it takes the unit of the side it falls.
    def test_faults_model(self):
        points_1 = [[0, 0], [1, 0], [3, 0], [0, 0.5], [5.0/6.0, 0.5], [2.0/3.0, 1], [3, 1], [0.0, 1.5], [0.5, 1.5], [1.0/3.0, 2], [3, 2], [0, 2.5], [1.0/6.0, 2.5], [0, 3], [3,3]]
        points_2 = [[0, 0], [2, 0], [3, 0], [0, 0.5], [2, 0.5],       [2, 1],       [3, 1], [0, 1.5],   [2, 1.5],   [2, 2],       [3, 2], [0, 2.5], [2, 2.5],       [0, 3], [2, 3], [3, 3]]
        
        pols_1 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 14, 13, 12]], [[11, 12, 13]]]
        pols_2 = [[[0, 1, 4, 3]], [[1, 2, 6, 5, 4]], [[3, 4, 5, 8, 7]], [[5, 6, 10, 9, 8]], [[7, 8, 9, 12, 11]], [[9, 10, 15, 14, 12]], [[11, 12, 14, 13]]]
        
        units = ["unit1", "unit1", "unit2", "unit2", "unit3", "unit3", "unit4"]
        
        model = cpp.Model([0,0,0,3,3,0],[0,0,0,3,3,0],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, [], [], []], ["s2", 2, points_2, pols_2, units, [], [], []]], {})
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
        model = cpp.Model([0,0,0,2,2,2],[0,0,0,2,2,2],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, []], ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" })
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3], {}, {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_1, [], [], []]], {})
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],{}, {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {})
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
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3], {}, {}, [["A-A", 1, points_1, polygons_1, units_1, [], [], []], ["B-B", 2, points_2, polygons_1, units_2, [], [], []]], {})
        model.make_matches()
        
    def test_faults_above(self):
        this_dir, this_filename = os.path.split(__file__)
        m = model.model_from_file(os.path.join(this_dir, 'test_files', 'Modelo_Hidro.json'))
    
    def test_faults_veins_fractures(self):
        this_dir, this_filename = os.path.split(__file__)
        m = model.model_from_file(os.path.join(this_dir, 'test_files', 'talud_faults_veins_fractures.json'))
        self.assertEqual(len(m.faults), 1)
        self.assertEqual(len(m.not_extended_faults), 1)
        self.assertEqual(len(m.fracts), 10)
        self.assertEqual(len(m.not_extended_fracts), 10)
        self.assertEqual(len(m.veins), 2)
        self.assertEqual(len(m.not_extended_veins), 2)

    # Test that you can load and test aburra Valley. Test grids and volumes.
    def test_aburra_valley(self):
        this_dir, this_filename = os.path.split(__file__)
        f = open(os.path.join(this_dir, 'test_files', 'aburra_version1.json'))
        m = model.GeologicalModel(json.loads(f.read()))
        m.height([813487.938015, 1164500.0])
        m.height([80000, 1100000])
        bbox = m.bbox
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
        self.assertEqual(map( lambda x: x[1], srt ),[117, 43, 68, 7, 43, 5, 1, 52, 1, 6])
       
        # Then test the fdm refined grid.
        ref_grid = utils.generate_fdm_grid(query_func, bbox, 5, 5)
        
        self.assertEqual(ref_grid['grid'].shape[3], 3)
        self.assertGreaterEqual(ref_grid['grid'].shape[0], ref_grid['units'].shape[0])
        self.assertGreaterEqual(ref_grid['grid'].shape[1], ref_grid['units'].shape[1])
        self.assertGreaterEqual(ref_grid['grid'].shape[2], ref_grid['units'].shape[2])
        
        grid = utils.generate_octtree_grid(query_func, bbox, 3, 3, 3)
        vols, elems = utils.octtree_volume_calculation(query_func, bbox, 10, 2)
        self.assertAlmostEqual(vols[t["Anfibolitas"]]/1e13, 2.71, 2)
        
        # Calculate bounded.
        verts, triangs = isosurfaces.calculate_isosurface(m, "Anfibolitas", 50 )
        self.assertEqual(len(verts), 11054)
        self.assertEqual(len(triangs), 22104)
        
        # Calculate unbounded
        verts, triangs = isosurfaces.calculate_isosurface(m, "Anfibolitas", 50, False )
        self.assertEqual(len(verts), 8892)
        self.assertEqual(len(triangs), 17386)
        
        # Filter by normal.
        verts, triangs = isosurfaces.calculate_isosurface(m, "Anfibolitas", 50, False, True, True )
        self.assertEqual(len(verts), 5864)
        self.assertEqual(len(triangs), 10529)
        
        # Filter by normal, negative.
        verts, triangs = isosurfaces.calculate_isosurface(m, "Anfibolitas", 50, False, True, False )
        self.assertEqual(len(verts), 4612)
        self.assertEqual(len(triangs), 8546)

    def test_modflow(self):
        
        this_dir, this_filename = os.path.split(__file__)
        fn = os.path.join(this_dir, 'test_files', 'Modelo_Hidro.json')
        geo_model = model.model_from_file(fn)

        Geo_Data = True
        Graph = True

        Rows = 50;  Cols = 25; Layers = 10; Angle = 30; DZ = 1.0
        Units = geo_model.units
        Kh = np.arange(len(Units))*0
        ani = np.ones(len(Units))
        Kz = Kh/10.0
        Act_Uni = np.ones(len(Units))

        Faults = geo_model.faults.keys()
        Kh_f = np.arange(len(Faults))+1.0
        ani_f = np.arange(len(Faults))+2.0
        Kz_f = Kh_f/10.0
        Act_faults = np.ones(len(Faults))
        Act_faults[0]=0
        Act_faults[2]=0
        Alg = 'adaptive'

        units_data = { Units[k]: (Kh[k],ani[k],Kz[k],Act_Uni[k]) for k in range(len(Units)) }
        faults_data = { Faults[k]: (Kh_f[k],ani_f[k],Kz_f[k],Act_faults[k]) for k in range(len(Faults)) }
        Bounding_Box = geo_model.bbox
        file_name = 'Files_Test'

        if False:
            modflow.create_modflow_inputs(file_name,geo_model,units_data,
            length_units=2, rows=Rows, cols=Cols,layers=Layers,
            bbox  = Bounding_Box, angle = Angle, dz_min=DZ, time_units=4,
            algorithm='adaptive',faults_data=faults_data)

        os.system('rm ' + file_name + '*')

        # Test LAYER CORRECTION
        Pos_List = [(0,0)]; Rows=Cols=2
        Mat_Order = np.zeros((Rows,Cols))
        Z_Bool_Top = np.array([[True,False],[False,False]])

        Z_top = np.array([[1.0, -0.6],[1.0,-2.45]])
        Max_Tan = np.tan(45 * np.pi/180.0)
        dz_min=0.25; dX=1.5; dY=2.0

        Z_Bottoms = [np.array([[1.0-dz_min, -1.0],
            [1.0-dz_min-dY*Max_Tan-0.1, -2.76]])]
        Layer=0

        #print Z_Bottoms
        #print Z_Bool_Top, '\n'
        modflow.Layer_Correction(Pos_List,Mat_Order,Z_Bool_Top,Z_top,Z_Bottoms,Layer,
            Max_Tan,Rows,Cols,dX,dY,dz_min)
        #print Z_Bottoms
        #print Z_Bool_Top
        #print Pos_List
        #print Mat_Order

        self.assertGreaterEqual(1E-5,abs(Z_Bottoms[0][0,0]-0.75))
        self.assertGreaterEqual(1E-5,abs(Z_Bottoms[0][0,1]+0.85))
        self.assertGreaterEqual(1E-5,abs(Z_Bottoms[0][1,0]+1.25))
        self.assertGreaterEqual(1E-5,abs(Z_Bottoms[0][1,1]+2.75))

        self.assertEqual(Mat_Order[0,1],1)
        self.assertEqual(Mat_Order[1,0],1)
        self.assertEqual(Mat_Order[1,1],2)

        # find_unit_limits test
        N = 0
        MIN = geo_model.bbox[0]; DELTA = geo_model.bbox[3]-MIN
        px = DELTA*np.random.rand(N) + MIN

        MIN = geo_model.bbox[1]; DELTA = geo_model.bbox[4]-MIN
        py = DELTA*np.random.rand(N) + MIN

        MIN = geo_model.bbox[2]

        for k in range(N):
            DELTA = geo_model.height((px[k], px[k])) - MIN
            Z = DELTA*np.random.rand(2) + MIN
            z_max = np.max(Z); z_min = np.min(Z)
            z1,b1 = modflow.find_unit_limits(geo_model, px[k], px[k], z_max, z_min, 1E-3)
            z2,b2 = geo_model.find_unit_limits(px[k], px[k], z_max, z_min, 1E-3)

            self.assertEqual(b1,b2)
            self.assertAlmostEqual(z1-z2,0)


    def test_faults(self):

        this_dir, this_filename = os.path.split(__file__)
        # Modelo de Argos
        fn = os.path.join(this_dir, 'test_files', 'Modelo_Argos.json')
        mfile = open(fn)
        geo_model = model.GeologicalModel(json.loads(mfile.read()),delete=False)

        Geo_Data = True
        Graph = True

        plane_1 = [[1066800,894300,960],[1067100,894300,960],[1067050.,894496.31385359,1069.93575801],[1066850, 894518.12650399,1082.15084223]]
        plane_2 = [[1066800,894300,1150],[1067200,894300,1150],[ 1.06727500e+06,8.94584037e+05,9.45493346e+02],[1.06683000e+06,8.94523172e+05,9.89316200e+02]]
        plane_3 = [[1066800,894300,1050],[1067100,894300,1050],[ 1.06727500e+06,8.94584037e+05,9.1050],[1.06683000e+06,8.94523172e+05,1050]]

        planes = [plane_1,plane_2,plane_3]

        start_time = datetime.now()
        Fault_int = geo_model.intersect_planes(planes)
        end_time = datetime.now()
        total_time = end_time - start_time
        #print 'Total plane time ARGOS: ', total_time.total_seconds(), ' s'

        # planes
        faults_size=[2,5,3,3]
        lines_size=[[6,9],[8,5,8,6,10],[15,19,16],[7,12,9]]
        
        c_fp=-1
        for name,fp in Fault_int.iteritems():
            c_fp+=1
            self.assertEqual(len(fp),faults_size[c_fp])
            #print '\n', len(fp)
            c_ls=0
            for ls in fp:
                self.assertEqual(len(ls),lines_size[c_fp][c_ls])
                c_ls+=1
                #print len(ls),


        # Modelo Hidrogeológico
        fn = os.path.join(this_dir, 'test_files', 'Modelo_Hidro.json')
        mfile = open(fn)
        geo_model = model.GeologicalModel(json.loads(mfile.read()),delete=False)

        plane_1 = [[1062218,1063707,2100],[1074654,1063707,2100],[ 1072218.,1073707.,2100.],[1062218,1075106,2100]]
        plane_2 = [[1062218,1063707,2350],[1074654,1063707,2350],[ 1070218.,1071707.,2350.],[1062818,1075106,2350]]
        plane_3 = [[1062218,1063707,2000],[1070000,1063707,2000],[ 1072218.,1073707.,2000.],[1068400,1075106,2000]]
        planes = [plane_1,plane_2,plane_3]

        start_time = datetime.now()
        Fault_int = geo_model.intersect_planes(planes)
        end_time = datetime.now()
        total_time = end_time - start_time
        #print 'Total plane time HIDRO: ', total_time.total_seconds(), ' s'        
        faults_size=[5,3,4,3]
        #lines_size=[[9,5,23,61,19],[16,10,12],[4,13,18,30],[58,39,60]]
        lines_size=[[9,5,23,61,19],[16,10,12],[3,13,15,30],[58,39,60]]
        c_fp=-1
        for name,fp in Fault_int.iteritems():
            c_fp+=1
            #print '\n', len(fp)
            self.assertEqual(len(fp),faults_size[c_fp])            
            c_ls=0
            for ls in fp:
                #print len(ls),
                self.assertEqual(len(ls),lines_size[c_fp][c_ls])
                c_ls+=1
 
        print '\n'
        # topography
        topo = geo_model.geojson[u'features'][0][u'transform']
        start_time = datetime.now()
        Fault_int = geo_model.intersect_topography( topo)
        end_time = datetime.now()
        total_time = end_time - start_time
        #print 'Total topography time HIDRO: ', total_time.total_seconds(), ' s', '\n'
        faults_size=[7,8,1,1]
        lines_size=[[898,35,44,3,153,11,47],[12,39,40,9,27,7,126,10],[565],[1077]]
        c_fp=-1
        for name,fp in Fault_int.iteritems():
            c_fp+=1
            #print '\n', len(fp)
            self.assertEqual(len(fp),faults_size[c_fp])            
            c_ls=0
            for ls in fp:
                #print len(ls),
                self.assertEqual(len(ls),lines_size[c_fp][c_ls])
                c_ls+=1

        pt_a = np.array(Fault_int.values()[0][0][187])
        pt_b = np.array([1063853.65210747,1065178.47682209])
        self.assertAlmostEqual(np.linalg.norm(pt_a-pt_b),0.0)
        pt_a = np.array(Fault_int.values()[1][6][87])
        pt_b = np.array([1072653.34031817,1072775.52557088])
        self.assertAlmostEqual(np.linalg.norm(pt_a-pt_b),0.0)
        pt_a = np.array(Fault_int.values()[2][0][246])
        pt_b = np.array([1064408.43639531,1069913.85073445])
        self.assertAlmostEqual(np.linalg.norm(pt_a-pt_b),0.0)
        pt_a = np.array(Fault_int.values()[3][0][1076])
        pt_b = np.array([1068678.98780714,1074722.21063005])
        self.assertAlmostEqual(np.linalg.norm(pt_a-pt_b),0.0)

        # Modelo Aburra
        fn = os.path.join(this_dir, 'test_files', 'aburra_version2.json')
        mfile = open(fn)
        geo_model = model.GeologicalModel(json.loads(mfile.read()),delete=False)

        val_z = -20000
        plane_1 = [[813449,1164500 ,val_z],[863000,1164500 ,val_z],[863000,1205000 ,val_z],[813449,1205000 ,val_z]]
        val_z = -2000
        plane_2 = [[813449,1164500 ,val_z],[863000,1164500 ,val_z],[863000,1205000 ,val_z],[813449,1205000 ,val_z]]
        val_z = 000
        plane_3 = [[813449,1164500 ,val_z],[863000,1164500 ,val_z],[863000,1205000 ,val_z],[813449,1205000 ,val_z]]
        planes = [plane_1,plane_2,plane_3]

        start_time = datetime.now()
        Fault_int = geo_model.intersect_planes(planes)
        end_time = datetime.now()
        total_time = end_time - start_time
        #print 'Total plane time ABURRA: ', total_time.total_seconds(), ' s'

        faults_size=[3,2,1,2,2,0,0,1,2]
        lines_size=[[24,24,25],[27,27],[23],[9,9],[7,7],[],[],[9],[8,10]]

        c_fp=-1
        for name,fp in Fault_int.iteritems():
            c_fp+=1
            self.assertEqual(len(fp),faults_size[c_fp])
            c_ls=0
            for ls in fp:
                self.assertEqual(len(ls),lines_size[c_fp][c_ls])
                c_ls+=1

        eps_z = 0.0
        eps_z2 = 1e-25
        planes=[[[1,1,eps_z],[-1,1,eps_z],[-1,-1,eps_z],[1,-1,eps_z]]]
        Faults={'Triangle_1':[[[1,0,0],[0,1,eps_z2],[0,0,1]]],'Triangle_2':[[[1,0,0],[0,0,-1],[0,1,eps_z2]]]}
        Fault_int = cpp.find_faults_intersection(Faults,planes)
        self.assertEqual(len(Fault_int.values()[0][0]),2)
        self.assertEqual(len(Fault_int.values()[1]),0)


        lines = [[0,0],[0.5,0],[1,0],[0.75,0.5],[0.25,0.5],[0.25,0.75],[0.35,1.],[0.75,1.],[1.,1.],[0.75,0.75],[1.,0.75]]
        segment=[]
        for k in range(len(lines)-1):
            segment.append([lines[k],lines[k+1]])

        pos = [[i] for i in range(len(segment))]
        shuffle(pos)
        segment_rand=copy.deepcopy(segment)
        for i in range(len(segment)):
            segment_rand[i] = segment[pos[i][0]]

        output = map(np.array, (cpp.join_lines_tree_test(segment_rand))[0])
        lines = map(np.array,lines)

        for k in range(len(lines)):
            self.assertAlmostEqual(np.linalg.norm(output[k]-lines[k]),0.0)

   
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
        self.assertEqual( str( cm.exception ), "Could not extend line." )
        
        with self.assertRaises( cpp.GeomodelrException ) as cm:
            cpp.extend_line( False, [0, 0, 1, 1], line )
        self.assertEqual( str( cm.exception ), "Could not extend line." )
        
        with self.assertRaises( cpp.GeomodelrException ) as cm:
            cpp.extend_line( True, [-1e-20, -1e-20, 1, 1], line )
        self.assertEqual( str( cm.exception ), "The line actually goes to the bounds and does not need modification." )
        
        line=[(0, 0), (-1e-20, -1e-20), (.5, .5), (1., 1.)] 
        with self.assertRaises( cpp.GeomodelrException ) as cm:
            cpp.extend_line( True, [-1e-40, -1e-40, 1, 1], line )
        self.assertEqual( str( cm.exception ), "Could not determine direction of line." )
        
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
        model = cpp.Model([0,0,0,4,4,4],[0,0,0,4,4,4],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False]]], 
                                                                 ["s2", 2, points_2, pols_2, units, faults_2, fname, [[0, False]]]], { "F1": "FAULT" })
        model.make_matches()
        
        ext = model.not_extended_faults["F1"]
        fau = model.faults["F1"]
        
        self.assertEqual(len(ext), 35)
        self.assertEqual(len(fau), 39)
        sext = set(sum(map(list,ext), []))
        sfau = set(sum(map(list,fau), []))
        self.assertEqual(len( sfau - sext ), 2)
        self.assertEqual(len( sext ), 28)
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False]]], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, [[0, False]]]], { "F1": "FAULT" })
        model.make_matches()
        ext2 = model.not_extended_faults["F1"]
        self.assertEqual(len(ext2), 36)
        self.assertEqual(len(model.faults["F1"]), 42)
        sext2 = set(sum(map(list,ext2), []))
        self.assertEqual(len(sext2), 28)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False]]], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" })
        model.make_matches()
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 39)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, []], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" })
        model.make_matches()
        
        f2 = model.not_extended_faults["F1"]
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 36)
        
        faults_1 = [[1, 4, 5, 8, 9, 12, 13]]
        faults_2 = [[1, 4, 5, 8, 9, 12, 14]]
        fname = ["F1"]
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], {}, {}, 
                          [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False], [0, True]]], 
                           ["s2", 2, points_2, pols_2, units, faults_2, fname, [[0, False], [0, True]]]], { "F1": "FAULT" })
        model.make_matches()
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 48)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, [[0, False], [0, True]]], 
                                                                    ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" })
        model.make_matches()
        
        self.assertEqual(len(model.not_extended_faults["F1"]), 35)
        self.assertEqual(len(model.faults["F1"]), 42)
        
        model = cpp.Model([-1,-1,-1,4,4,4],[-1,-1,-1,4,4,4],[0, 0], [1, 0], {}, {}, [["s1", 1, points_1, pols_1, units, faults_1, fname, []], ["s2", 2, points_2, pols_2, units, faults_2, fname, []]], { "F1": "FAULT" } )
        model.make_matches()
        
        self.assertEqual(len(model.not_extended_faults["F1"]), 36)
        self.assertEqual(len(model.faults["F1"]), 36)
        
    def test_aligned(self):
        # Evaluate a model that has a hole in the middle.
        points_1 = [[0, 0], [3, 0], [3, 3], [0, 3], [1, 1], [2, 1], [2, 2], [1, 2]]
        points_2 = [[0, 0], [3, 0], [3, 3], [0, 3]]
        
        polygons_1 = [[[0, 1, 2, 3], [7, 6, 5, 4]], [[4, 5, 6, 7]]]
        polygons_2 = [[[0, 1, 2, 3]]] 
        
        units_1 = ["unit1", "unit2"]
        units_2 = ["unit1"]
        
        model = cpp.Model([0,0,0,3,3,3],[0,0,0,3,3,3],[0, 3], [math.sqrt(2)/2.0, math.sqrt(2)/2.0], {}, {}, [["A-A", -1, points_1, polygons_1, units_1, [], [], []], ["B-B", 1, points_2, polygons_2, units_2, [], [], []]], {})
        model.make_matches()
        self.assertEqual(model.model_point([1.5, 0.1, 1.5]), (1.5, 1.5, 0.1))
        cls = model.closest((1.5, 1.1, 1.5))
        self.assertEqual(cls[0], 'unit2')
        self.assertAlmostEqual(cls[1], 0.1)
        ip = model.inverse_point((1.5, 1.1, 1.5))
        clsa = model.closest_aligned(ip)
        self.assertEqual(clsa[0], 'unit2')
        self.assertAlmostEqual(clsa[1], 0.1)

def main(args=None):
    unittest.main()

if __name__ == '__main__':
    main()
