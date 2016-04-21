from dustcurve import pixclass
import h5py
import healpy
import emcee 
from dustcurve import model 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# the model has 12 parameters; we'll use 50 walkers and 500 steps each
ndim = 12
nwalkers = 100
nsteps = 1000

# set up the walkers in a "Gaussian ball" around the literature estimate for distance to Cepheus cloud (distance mod=10)
#ls_result = np.arange(4,16,1)
ls_result=np.arange(4,16,1)
starting_positions = [ls_result + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# set up the sampler object
sampler = emcee.EnsembleSampler(nwalkers, ndim, model.log_posterior, args=['pixel 1024-5753839'])
                                
# run the sampler. We use iPython's %time directive to tell us 
# how long it took (in a script, you would leave out "%time")
sampler.run_mcmc(starting_positions, nsteps)

print('Done')

# Burn off initial steps
samples = sampler.chain[:, 500:, :].reshape((-1, ndim))

# Plot traces to see how much you need to burn off, and to check for convergence 
plt.plot(sampler.chain[:, :, 0].T, color="k", alpha=0.3)
plt.axhline(F_true, color='#4682b4')
plt.ylabel('Parameter 1')
plt.xlabel('Number of Steps')




fig, (ax_d1, ax_d2) = plt.subplots(2)
ax_d1.set(ylabel='d1')
ax_d2.set(ylabel='d2')
for i in range(10):
    sns.tsplot(sampler.chain[i,:,0], ax=ax_d1)
    sns.tsplot(sampler.chain[i,:,1], ax=ax_d2)
    
plt.show()
