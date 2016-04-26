
from unittest import TestCase
import unittest
from dustcurve import pixclass
from dustcurve import model
import numpy as np
import emcee

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
      ndim=12
      nsteps=1
      nwalkers=1
      pos_array=[0,0,0,0,0,0,7.75,0,0,0,14,0]
      std_array=[1. for i in range(ndim)]
      starting_positions = emcee.utils.sample_ball((pos_array),(std_array),nwalkers) #set up the initial 
      file='simulated_data.h5'
      sampler = emcee.MHSampler(np.diagflat(np.ones(ndim)), ndim, model.log_posterior, args=[file])
      (pos,lnprob)=sampler.run_mcmc(starting_positions[0],nsteps)
      #check that the line integral you're getting is above 215, the approximate probability you would get if you summed
      #over the "true" reddening profile given by the above distance array
      self.assertTrue(lnprob>215.00)
      
if __name__ == '__main__':
    unittest.main()

