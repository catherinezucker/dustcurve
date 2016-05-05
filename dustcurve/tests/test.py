from unittest import TestCase
import unittest
from dustcurve import pixclass
from dustcurve import model
import numpy as np
import emcee
from dustcurve import io

#test pixstars class pt 1
class PixStarsTestPost(unittest.TestCase): 
   def TestPixStarsPost(self):
      """Tests PixClass class"""
      #Create PixClass object
      pixTest = pixclass.PixStars('simulated_data.h5')
      # Is the stellar posterior array the right shape? 
      self.assertTrue(pixTest.get_p()[0,:,:].shape==(700,120))

      
#test pixstars class pt 2
class PixStarsTestCO(unittest.TestCase): 
   def TestPixStarsCO(self):
      """Tests PixClass class"""
      #Create PixClass object
      pixTest = pixclass.PixStars('simulated_data.h5')
      #Is the CO data array the right shape?
      self.assertTrue(pixTest.get_co()[0,:].shape==(12,))
      

#test io module pt 2
class IOTestPrior(unittest.TestCase): 
   def TestIOPrior(self):
      """Tests PixClass class"""
      #Call io.fetch_args
      ldata,pdata = io.fetch_args('simulated_data.h5', [4,19], 1.0)
      # Are the arguments returned in the right format? 
      self.assertTrue(pdata.shape==(2,))

#test io module
class IOTestLikelihood(unittest.TestCase): 
   def TestIOLikelihood(self):
      """Tests PixClass class"""
      #Call io.fetch_args
      ldata,pdata = io.fetch_args('simulated_data.h5', [4,19], 1.0)
      #Is the ldata the right type?
      self.assertTrue(type(ldata)==tuple)

#test model module pt 1 
class ModelTestCase(unittest.TestCase): 
   def TestLineIntegral(self):
        ndim=12
        nwalkers = 50
        nsteps = 1
        ntemps=5

        #fetch the required likelihood and prior arguments for PTSampler
        ldata,pdata=io.fetch_args('simulated_data.h5',[4,19],1.0)

        #set up the walkers at a location where we know what probability should be returned
        result=[5,5,5,5,5,5,7.75,8,8,8,14,15]

        #set up starting positions, all at same distance with known probability 
        starting_positions = [[result for i in range(nwalkers)] for j in range(ntemps)]

        #set up the sampler object
        sampler = emcee.PTSampler(ntemps, nwalkers, ndim, model.log_likelihood, model.log_prior, loglargs=(ldata), logpargs=[pdata])
     
        sampler.run_mcmc(starting_positions,nsteps)
        lnprob=sampler.lnprobability
        
        #check that the line integral you're getting is above XX, the minimum reasonable probability you would get if you summed
        #over the "true" reddening profile given by the above distance array
        
        self.assertTrue(lnprob[0,0,0]>5.00)

#test model module pt 2 
class ModelTestCasePrior(unittest.TestCase): 
   def TestLineIntegral(self):
        ndim=12
        nwalkers = 50
        nsteps = 1
        ntemps=5

        #fetch the required likelihood and prior arguments for PTSampler
        ldata,pdata=io.fetch_args('simulated_data.h5',[4,19],1.0)

        #set up the walkers outside the prior bounds, so we should hopefully return -np.inf
        result=[0,5,5,5,5,3,7.75,8,8,8,14,20]

        #set up starting positions, all at same distance with known probability (-np.inf)
        starting_positions = [[result for i in range(nwalkers)] for j in range(ntemps)]

        #set up the sampler object
        sampler = emcee.PTSampler(ntemps, nwalkers, ndim, model.log_likelihood, model.log_prior, loglargs=(ldata), logpargs=[pdata])
     
        sampler.run_mcmc(starting_positions,nsteps)
        lnprob=sampler.lnprobability
        
        #check that the line integral you're getting is -np.inf since we instantiated our walkers outside the prior bounds         
        self.assertTrue(lnprob[0,0,0]==-np.inf)
    
      
if __name__ == '__main__':
    unittest.main()

