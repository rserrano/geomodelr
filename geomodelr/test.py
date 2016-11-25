"""
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
"""


import unittest
import os
import model
import json
import faults
import shared
import cpp
import numpy as np
from numpy import linalg as la

class TestGeoModelR(unittest.TestCase):
    def test_create_model(self):
        pass
    def test_faultplane_for_lines(self):
        """
        Tests function fault plane for lines.
        """
        la = [(0, 0, 0), (1,1,0), (2,2,0), (3,3,0)]
        lb = [(0, 3, 1), (1,2,1), (2,1,1), (3,0,1)]
        # This is a weird testcase that should work. Star in two directions.
        self.assertEqual(faults.faultplane_for_lines(la, lb), [(0, 1, 4), (1, 4, 5), (1, 2, 5), (2, 3, 5), (3, 5, 6), (3, 6, 7)])
        
        # This is a weird testcase that should not work.
        # Open circles in different directions.
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1,-2,1), (-2,-1,1), (-2,1,1),  (-1, 2,1), (1, 2,1), (2, 1,1), (2,-1,1), (1, -2, 1)]
        with self.assertRaises(shared.GeometryException):
            faults.faultplane_for_lines(la, lb)
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        
        # Open circles in different directions but one on top of the other.
        with self.assertRaises(shared.GeometryException):
            faults.faultplane_for_lines(la, lb)
        
        # One starts when the other finishes.
        la = [(-1,-2,0), (-2,-1,0), (-2,1,0),  (-1, 2,0), (1, 2,0), (2, 1,0), (2,-1,0), (1, -2, 0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        self.assertEqual(faults.faultplane_for_lines(la, lb), [(0,1,  8), (1,  8,  9), (1,  2,  9), 
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
        
        self.assertEqual(faults.faultplane_for_lines(la, lb), [(0, 8, 9), (0, 9, 10), (0, 1, 10), 
                                                               (1, 10, 11), (1, 11, 12), (1, 12, 13), 
                                                               (1, 2, 13), (2, 3, 13), (3, 13, 14), 
                                                               (3, 4, 14), (4, 5, 14), (5, 14, 15), 
                                                               (5, 6, 15), (6, 7, 15), (7, 15, 16)])
        
    
    def test_next_and_check(self):
        """
        Test helper function next_and_check.
        """
        self.assertEqual( faults.next_and_check( (1,2,6), (2,6), 4), ( (1,6), (1,2) ) )
        self.assertEqual( faults.next_and_check( (1,2,6), (1,6), 4), ( (2,6), (1,2) ) )
        self.assertEqual( faults.next_and_check( (1,5,6), (1,6), 4), ( (1,5), (5,6) ) )
        self.assertEqual( faults.next_and_check( (1,5,6), (1,5), 4), ( (1,6), (5,6) ) )

    def test_aburra_valley(self):
        """
        Test that you can load and test aburra Valley.
        """
        this_dir, this_filename = os.path.split(__file__)
        f = open(os.path.join(this_dir, 'test_files', 'aburra_version1.json'))
        m = model.GeologicalModel(json.loads(f.read()))
        m.calc_cache()
    
    def test_sections(self):
        """
        Test that you can create cross sections and that bugs are catched.
        """
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]]]
        units = ["unit1"]
        lines = []
        lnames = []
        
        # Try the simplest one.
        section = cpp.Section(0, points, polygons, units, lines, lnames)
        
        # Try a bad index.
        polygons = [[[0, 1, 2, 3, 6]]]
        with self.assertRaises(IndexError):
            section = cpp.Section(1, points, polygons, units, lines, lnames)
        
        # Try a bad number.
        points = [0, [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]]]
        with self.assertRaises(TypeError):
            section = cpp.Section(2, points, polygons, units, lines, lnames)
        
        # Try a bad number of polygons.
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1]]
        polygons = [[[0, 1, 2, 3, 4]], []]
        with self.assertRaises(IndexError):
            section = cpp.Section(3, points, polygons, units, lines, lnames)
        
        # Try a bad polygon. The crossing lines goes from [1,1] to [0,0].
        points = [[0, 0], [1, 0], [0, 1], [1, 1]]
        polygons = [[[0, 1, 2, 3]]]
        section = cpp.Section(4, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['polygons'], 0)
        
        # Test polygons with holes.
        points = [[0, 0], [1, 0], [0, 1], [1, 1]]
        polygons = [[[0, 1, 2, 3]]]
        section = cpp.Section(5, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['polygons'], 0)
        
        # Test lines get created.
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]]
        lines = [[0, 2, 4], [1, 3]]
        lnames = ['fault1', 'fault2']
        section = cpp.Section(6, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['lines'], 2)
        
        # Test polygons with holes.
        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]]]
        units = ['unit1', 'unit2']
        section = cpp.Section(7, points, polygons, units, lines, lnames)
        self.assertEqual(section.info()['polygons'], 2)
        
    def test_section_closest(self):
        points = [[0, 0], [1, 0], [2, 1], [1, 1], [0, 1], [0.25, 0.25], 
                  [0.75, 0.25], [0.75, 0.75], [0.25, 0.75], [2, 0], [2, 2], [0, 2]]

        polygons = [[[0, 1, 2, 3, 4], [5, 8, 7, 6]], [[5, 6, 7, 8]], [[2, 1, 9]], [[4, 3, 2, 10, 11]]]
        units = ['unit1', 'unit2', 'unit3', 'unit4']
        section = cpp.Section(8, points, polygons, units, [], [])
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
        section = cpp.Section(9, points, polygons, units, [], [])
        self.assertEqual(section.info()['polygons'], 4)
        self.assertEqual(section.closest([0.5, 0.5]), (1, 'unit2'))
        self.assertEqual(section.closest([-0.5, -0.5]), (1, 'unit2'))
        self.assertEqual(section.closest([1.5, -0.5]), (2, 'unit3'))
        self.assertEqual(section.closest([1.5, 0.75]), (2, 'unit3'))
        self.assertEqual(section.closest([1.5, 0.25]), (2, 'unit3'))
        self.assertEqual(section.closest([-0.001, 1.001]), (3, 'unit4'))
        
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0.9], [0.9, 0.9], [0.9, 0.1], [0, 0.1],
                  [0.25, 0.25], [0.9, 0.25], [0.9, 0.75], [0.25, 0.75] ]
        polygons = [[[0, 1, 2, 3, 4, 5, 6, 7]], [[8, 9, 10, 11]]]
        units = ["unit1", "unit2"]
        section = cpp.Section(10, points, polygons, units, [], [])
        self.assertEqual(section.closest([0, 0.5]), (1, 'unit2'))
        self.assertEqual(section.closest([0, 0.3]), (0, 'unit1'))

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
        
        model = cpp.Model([1, 0], [0, 1], [(11, points_1, polygons_1, units_1, [], []), (12, points_2, polygons_2, units_2, [], [])])
        
        model.make_matches()
        self.assertEqual(model.matches, [[(0, 0), (0, 3), (2, 2), (3, 4)]])
        model.matches = [[(0, 0), (2, 2), (3, 4)]]
        self.assertEqual(model.matches, [[(0, 0), (2, 2), (3, 4)]])
    
    def test_possible_closest(self):
        # Evaluate simple models.
        points_1   = [[0, 0], [3, 0], [0, 1], [2, 1], [3, 1], [0, 2], [2, 2], [3, 2]]
        points_2   = [[0, 0], [3, 0], [0, 1], [1, 1], [3, 1], [0, 2], [1, 2], [3, 2]]
        polygons_1 = [[[0, 1, 4, 3, 2]], [[2, 3, 6, 5]], [[3, 4, 7, 6]]]
        units_1 = ["unit1", "unit2", "unit3"]
        
        model = cpp.Model([0, 0], [1, 0], [(1, points_1, polygons_1, units_1, [], []), (2, points_2, polygons_1, units_1, [], [])])
        model.make_matches()
        self.assertEqual(model.matches, [[(0, 0), (1, 1), (2, 2)]])
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
        
        # Evaluate all single.
        points_1   = [[0, 0], [3, 0], [0, 1], [2, 1], [3, 1], [0, 2], [2, 2], [3, 2]]
        points_2   = [[0, 0], [3, 0], [0, 1], [1, 1], [3, 1], [0, 2], [1, 2], [3, 2]]
        polygons_1 = [[[0, 1, 4, 3, 2]], [[2, 3, 6, 5]], [[3, 4, 7, 6]]]
        units_1 = ["unit1", "unit2", "unit3"]
        units_2 = ["unit4", "unit5", "unit6"]
        
        model = cpp.Model([0, 0], [1, 0], [(1, points_1, polygons_1, units_1, [], []), (2, points_2, polygons_1, units_2, [], [])])
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
        
        units_1 = ["unit1", "unit2", "unit3"]
        units_2 = ["unit1", "unit5", "unit6"]
        
        model = cpp.Model([0, 0], [1, 0], [(1, points_1, polygons_1, units_1, [], []), (2, points_2, polygons_1, units_2, [], [])])
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

if __name__ == '__main__':
    unittest.main()
