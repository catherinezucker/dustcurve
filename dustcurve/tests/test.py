from unittest import TestCase
import unittest
from dustcurve import pixclass
from dustcurve import model
import numpy as np
import emcee
from dustcurve import io
import os

def get_example_data_file_path(filename, data_dir='example_data'):
    # __file__ is the location of the source file currently in use (so
    # in this case io.py). We can use it as base path to construct
    # other paths from that should end up correct on other machines or
    # when the package is installed
    start = os.path.abspath(__file__)
    start_dir = os.path.dirname(start)
    # If you need to go up another directory (for example if you have
    # this function in your tests directory and your data is in the
    # package directory one level up) you can use
    up_dir = os.path.split(start_dir)[0]
    data_dir = os.path.join(up_dir, data_dir)
    return os.path.join(start_dir, data_dir, filename)

#test pixstars class pt 1
class PixStarsTestPost(unittest.TestCase):
   def TestPixStarsPost(self):
      """Tests PixClass class"""
      #Create PixClass object
      pixTest = pixclass.PixStars(get_example_data_file_path('simulated_data.h5'))
      # Is the stellar posterior array the right shape?
      self.assertTrue(pixTest.get_p()[0,:,:].shape==(700,120))


#test pixstars class pt 2
class PixStarsTestCO(unittest.TestCase):
   def TestPixStarsCO(self):
      """Tests PixClass class"""
      #Create PixClass object
      pixTest = pixclass.PixStars(get_example_data_file_path('simulated_data.h5'))
      #Is the CO data array the right shape?
      self.assertTrue(pixTest.get_co()[0,:].shape==(12,))


#test io module pt 2
class IOTestPrior(unittest.TestCase):
   def TestIOPrior(self):
      """Tests PixClass class"""
      #Call io.fetch_args
      ldata,pdata = io.fetch_args(get_example_data_file_path('simulated_data.h5'), [4,19], 1.0)
      # Are the arguments returned in the right format?
      self.assertTrue(pdata.shape==(2,))

#test io module
class IOTestLikelihood(unittest.TestCase):
   def TestIOLikelihood(self):
      """Tests PixClass class"""
      #Call io.fetch_args
      ldata,pdata = io.fetch_args(get_example_data_file_path('simulated_data.h5'), [4,19], 1.0)
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
        ldata,pdata=io.fetch_args(get_example_data_file_path('simulated_data.h5'),[4,19],1.0)

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

        self.assertTrue(lnprob[0,0,0]>100.0)

#test model module pt 2
class ModelTestCasePrior(unittest.TestCase):
   def TestLineIntegral(self):
        ndim=12
        nwalkers = 50
        nsteps = 1
        ntemps=5

        #fetch the required likelihood and prior arguments for PTSampler
        ldata,pdata=io.fetch_args(get_example_data_file_path('simulated_data.h5'),[4,19],1.0)

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
