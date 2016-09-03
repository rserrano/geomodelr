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
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1,-2,1), (-2,-1,1), (-2,1,1),  (-1, 2,1), (1, 2,1), (2, 1,1), (2,-1,1), (1, -2, 1)]
        # Open circles in different directions.
        with self.assertRaises(shared.GeometryException):
            faults.faultplane_for_lines(la, lb)
        la = [(-1, 2,0), (-2,1, 0), (-2,-1,0), (-1,-2,0), (1,-2,0), (2,-1,0), (2, 1,0), (1,2,0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        # Open circles in different directions but one on top of the other.
        with self.assertRaises(shared.GeometryException):
            faults.faultplane_for_lines(la, lb)
        la = [(-1,-2,0), (-2,-1,0), (-2,1,0),  (-1, 2,0), (1, 2,0), (2, 1,0), (2,-1,0), (1, -2, 0)]
        lb = [(-1, 2, 1), (-2, 3, 1), (-2, 5, 1), (-1, 6, 1), (1, 6, 1), (2, 5, 1), (2, 3, 1), (1, 2, 1)]
        self.assertEqual(faults.faultplane_for_lines(la, lb), [(0,1,  8), (1,  8,  9), (1,  2,  9), 
                                                                (2, 9,10), (2,  3, 10), (3, 10, 11), 
                                                                (3,11,12), (3,  4, 12), (4, 12, 13), 
                                                                (4, 5,13), (5, 13, 14), (5,  6, 14), 
                                                                (6,14,15), (6,  7, 15)])
    
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
