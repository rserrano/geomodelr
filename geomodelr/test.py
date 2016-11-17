import unittest
import os
import model
import json
import faults
import shared
import cpp

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
        

if __name__ == '__main__':
    unittest.main()
