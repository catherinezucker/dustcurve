
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
      ncoeff = 12
      ndist= 12
      ndim=24
      nwalkers = 50
      nsteps = 500
      
      fname='simulated_data.h5'
      pixObj=pixclass.PixStars(fname)
      co_array=np.asarray(pixObj.get_co()[:,:])
      post_array=np.asarray(pixObj.get_p()[:,:,:])
      nstars=pixObj.get_n_stars()
      data=(co_array,post_array,nstars) 

      allsamples = np.empty((1,ndim))
      pos_array=[4,4.5,4.75,5,5,5,7.75,8,8,8,14,15] #set off walkers at random distances within acceptable range
      pos_array.extend([1 for i in range(ncoeff)]) #set of walkers near expected gas to dust coefficient in literatre (0.03 mags of reddening per CO intensity unit K)
      std_array=[0.01 for i in range(ndim)]
      starting_positions = emcee.utils.sample_ball((pos_array),(std_array),nwalkers) #set up the initial position vectors for our walkers

      sampler = emcee.MHSampler(np.diagflat(np.ones(ndim)), ndim, model.log_posterior, args=(data))
      
      (pos,lnprob,rstate)=sampler.run_mcmc(starting_positions[0],nsteps)
      #check that the line integral you're getting is above 215, the approximate probability you would get if you summed
      #over the "true" reddening profile given by the above distance array
      self.assertTrue(lnprob>215.00)
      
if __name__ == '__main__':
    unittest.main()

