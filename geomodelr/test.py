import unittest
import faults
import shared
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

if __name__ == '__main__':
    unittest.main()
