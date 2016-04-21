
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
      pixTest = pixclass.PixStars('testdata.h5', 'pixel 1024-5753839')
      # Is the stellar posterior array the right shape? 
      self.assertTrue(pixTest.get_p()[0,:,:].shape==(120,700))
      
class ModelLineIntegralTestCase(unittest.TestCase): 
   def TestLineIntegral(self):
      """Tests Line Integral Function Critical to Model"""
      #create artificial arrays we fill with known values
      post_array=np.zeros((120,700))
      dist_array=np.array([4.3,5.5,6.2,7,8.3,9.1,10,10,11.45,12.9,13.7,14.01])
      coeff_array=np.ones((12))
      co_array=np.ones((12))*0.25
      dbins,redbins=model.convert_to_bins(co_array,dist_array,coeff_array)
       #manipulate the posterior array so we know exactly what the line integral should be, given the dist_array, co_array, coeff_array   
      post_array[0:2,0]=1
      post_array[2:12,25]=1
      post_array[12:17,50]=1
      post_array[17:24,75]=1
      post_array[24:34,100]=1
      post_array[34:40,125]=1
      post_array[40:48,150]=1
      post_array[48:59,200]=1
      post_array[59:71,225]=1
      post_array[71:77,250]=1
      post_array[77:80,275]=1
      post_array[80:-1,300]=1
      #test that the value of the integral over the stellar posterior array is 1 for this specific test case
      self.assertTrue(model.get_line_integral(get_line_integral(co_array,post_array,dist_array,coeff_array)==119.0))
      
if __name__ == '__main__':
    unittest.main()

