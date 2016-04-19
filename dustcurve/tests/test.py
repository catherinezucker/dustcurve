
from unittest import TestCase
import unittest
from dustcurve import pixclass 

# Test PixStars Class

class PixStarsTestCase(unittest.TestCase): 
   def TestPixStars(self):
      """Tests PixClass class"""
      #Create PixClass object
      pixObj = PixClass('testdata.h5', 'pixel 1024-5753839')
      # Is the stellar posterior array the right shape? 
