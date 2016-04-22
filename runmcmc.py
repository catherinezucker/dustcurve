from dustcurve import pixclass
import h5py
import healpy
import emcee 
from dustcurve import model 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# the model has 12 parameters; we'll use 50 walkers and 500 steps each
ndim = 12
nwalkers = 50
nsteps = 500

# set up the walkers in a "Gaussian ball" around the literature estimate for distance to Cepheus cloud (distance mod=10)
ls_result=[4.5,4.6,5,5.1,6.2,7.0,7.75,8.0,9,12,14,16]
starting_positions = [ls_result + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# set up the sampler object
sampler = emcee.EnsembleSampler(nwalkers, ndim, model.log_posterior, args=('simulated_data.h5', 'pixel0000'))

# run the sampler. We use iPython's %time directive to tell us 
# how long it took (in a script, you would leave out "%time")
testprob=sampler.run_mcmc(starting_positions, nsteps)

print('Done')

fig, (ax_d7, ax_d11) = plt.subplots(2)

ax_d7.set(ylabel='d7')
ax_d11.set(ylabel='d11')


for i in range(0,nwalkers):
    sns.tsplot(sampler.chain[i,:,6], ax=ax_d7)
    sns.tsplot(sampler.chain[i,:,10], ax=ax_d11)
    
# Burn off initial steps
samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
    
traces = samples.reshape(-1, ndim).T

parameter_samples = pd.DataFrame({'d1': traces[0], 'd2': traces[1], 'd3': traces[2], 'd4': traces[3], 'd5': traces[4], 'd6': traces[5], 'd7': traces[6], 'd8': traces[7], 'd9': traces[8], 'd10': traces[9], 'd11': traces[10], 'd12': traces[11]})


q = parameter_samples.quantile([0.16,0.50,0.84], axis=0)

                                            	
print("d7 = {:.2f} + {:.2f} - {:.2f}".format(q['d7'][0.50], 
                                            q['d7'][0.84]-q['d7'][0.50],
                                            q['d7'][0.50]-q['d7'][0.16]))
print("d11 = {:.2f} + {:.2f} - {:.2f}".format(q['d11'][0.50], 
                                            q['d11'][0.84]-q['d11'][0.50],
                                            q['d11'][0.50]-q['d11'][0.16]))
                                                           
plt.show()
