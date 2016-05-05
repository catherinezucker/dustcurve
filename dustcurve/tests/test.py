from unittest import TestCase
import unittest
from dustcurve import pixclass
from dustcurve import model
import numpy as np
import emcee
from dustcurve import io

#test pixstars class
class PixStarsTestCase(unittest.TestCase): 
   def TestPixStars(self):
      """Tests PixClass class"""
      #Create PixClass object
      pixTest = pixclass.PixStars('simulated_data.h5')
      # Is the stellar posterior array the right shape? 
      self.assertTrue(pixTest.get_p()[0,:,:].shape==(700,120))
      #Is the CO data array the right shape?
      self.assertTrue(pixTest.get_co()[0,:].shape==(12))
      
      
#test io module
class IOTestCase(unittest.TestCase): 
   def TestIO(self):
      """Tests PixClass class"""
      #Call io.fetch_args
      ldata,pdata = io.fetch_args('simulated_data.h5', [4,19], 1.0)
      # Are the arguments returned in the right format? 
      self.assertTrue(pdata.shape==(2,1))
      #Is the CO data array the right shape?
      self.assertTrue(pixTest.get_co()[0,:].shape==(12))
      

#test model module      
class ModelLineIntegralTestCase(unittest.TestCase): 
   def TestLineIntegral(self):
        ndim=12
        nwalkers = 50
        nsteps = 1

        #fetch the required likelihood and prior arguments for PTSampler
        ldata,pdata=io.fetch_args('simulated_data.h5',[4,19],1.0)

        #set up the walkers at a location where we know what probability should be returned
        result=[5,5,5,5,5,5,7.75,8,8,8,14,15]

        #set up starting positions, all at same distance 
        starting_positions = [[result for i in range(nwalkers)] for j in range(ntemps)]

        #set up the sampler object
        sampler = emcee.PTSampler(ntemps, nwalkers, ndim, model.log_likelihood, model.log_prior, loglargs=(ldata), logpargs=[pdata])
     
        (pos,lnprob,rstate)=sampler.run_mcmc(starting_positions,nsteps)
        #check that the line integral you're getting is above 215, the approximate probability you would get if you summed
        #over the "true" reddening profile given by the above distance array
        self.assertTrue(lnprob>215.00)
      
if __name__ == '__main__':
    unittest.main()

