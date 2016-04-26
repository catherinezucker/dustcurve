
from unittest import TestCase
import unittest
from dustcurve import pixclass
from dustcurve import model
import numpy as np

# Test PixStars Class

class PixStarsTestCase(unittest.TestCase): 
   def TestPixStars(self):
      """Tests PixClass class"""
      #Create PixClass object
      pixTest = pixclass.PixStars('simulated_data.h5')
      # Is the stellar posterior array the right shape? 
      self.assertTrue(pixTest.get_p()[0,:,:].shape==(700,120))
      
class ModelLineIntegralTestCase(unittest.TestCase): 
   def TestLineIntegral(self):
      """Tests Line Integral Function Critical to Model"""
      pixTest = pixclass.PixStars('simulated_data.h5')
      stack_array=np.array([0,0,0,0,0,0,2,0,0,0,3,0])
      for i in range(0,pixTest.get_n_stars()-1)
         co_array=np.vstack((co_array, stack_array))
      post_array=pixTest.get_p()
      coeff_array=np.ones((12))
      dist_array=np.array([0,0,0,0,0,0,7.75,0,0,0,14,0])
      #check that the line integral you're getting is above 215, the approximate probability you would get if you summed
      #over the "true" reddening profile given by the above distance array
      self.assertTrue(model.get_line_integral(co_array,post_array,dist_array,coeff_array)>215.00)
      
if __name__ == '__main__':
    unittest.main()

