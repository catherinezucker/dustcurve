from dustcurve import pixclass
import h5py
import healpy
import emcee 
from dustcurve import model 
import numpy as np

# the model has 12 parameters; we'll use 50 walkers and 500 steps each
ndim = 12
nwalkers = 50
nsteps = 50

# set up the walkers in a "Gaussian ball" around the literature estimate for distance to Cepheus cloud (distance mod=10)
ls_result = np.arange(4,16,1)
starting_positions = [ls_result + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

#pixObj=pixclass.PixStars('testdata.h5', 'pixel 1024-5753839')
#co_array=np.asarray(pixObj.get_co()[0,:])
#post_array=np.asarray(pixObj.get_p()[0,:,:])
#nstars=np.asarray(pixObj.get_n_stars())



# set up the sampler object
sampler = emcee.EnsembleSampler(nwalkers, ndim, model.log_posterior, args=['pixel 1024-5753839'])
                                
# run the sampler. We use iPython's %time directive to tell us 
# how long it took (in a script, you would leave out "%time")
sampler.run_mcmc(starting_positions, nsteps)

print('Done')