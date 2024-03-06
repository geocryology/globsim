import unittest
from globsim.interp import interpolate_level_data, calculate_weights, ele_interpolate
import numpy as np


class TestEleInterpolate(unittest.TestCase):

    def setUp(self):
        self.ele1 = np.array([[300,200,100],
                              [310, 210,110]], dtype='float')
        
    def test_simple(self):
        elev_diff, va, vb = ele_interpolate(self.ele1, 150, 3)
        self.assertEqual(1, va[0])
        self.assertEqual(4, va[1])
        self.assertEqual(2, vb[0])
        self.assertEqual(5, vb[1])
    

class TestLevelInterpolation(unittest.TestCase):
    
    def test_on_level(self):
        h = 200
        data = np.array([[4,8,16],
                         [31, 21, 11]], dtype='float')
        ele = np.array([[300,200,100],
                        [310, 210, 110]], dtype='float')
        ipol = interpolate_level_data(data, ele, h)
        self.assertEqual(8, ipol[0])

    def test_on_lowest_level(self):
        h = 100
        data = np.array([[4,8,16],
                         [31, 21, 11]], dtype='float')
        ele = np.array([[300,200,100],
                        [310, 210, 110]], dtype='float')
        ipol = interpolate_level_data(data, ele, h)
        self.assertEqual(16, ipol[0])

    def test_below_lowest_level(self):
        h = 50
        data = np.array([[4,8,16],
                         [31, 21, 11]], dtype='float')
        ele = np.array([[300,200,100],
                        [310, 210, 110]], dtype='float')
        ipol = interpolate_level_data(data, ele, h)
        self.assertEqual(16, ipol[0])

    def test_on_highest_level(self):
        h = 300
        data = np.array([[4,8,16],
                         [31, 21, 11]], dtype='float')
        ele = np.array([[300,200,100],
                        [310, 210, 110]], dtype='float')
        ipol = interpolate_level_data(data, ele, h)
        self.assertEqual(4, ipol[0])
    
    def test_backwards_elevations(self): 
        h = 125
        data = np.array([[4,8,16],[31, 21, 11]], dtype='float')
        ele = np.array([[100,200,300],[110, 210, 310]], dtype='float')
        ipol = interpolate_level_data(data, ele, h)
        self.assertEqual(ipol[0], 5)

    def test_two_stations(self): 
        h = np.array([125,250])
        data = np.array([[[4,8,16],
                          [31, 21, 11]],
                         [[4,8,16],
                          [31, 21, 11]]], dtype='float')
        ele = np.array([[[100,200,300],
                         [110, 210, 310]],
                        [[100,200,300],
                         [110, 210, 310]]], dtype='float')
        data = data.transpose(1,2,0)
        ele = ele.transpose(1,2,0)
        ipol = interpolate_level_data(data, ele, h)
        self.assertEqual(ipol[0,0], 5)
        self.assertEqual(ipol[0,1], 12)


if __name__ == "__main__":
    unittest.main()