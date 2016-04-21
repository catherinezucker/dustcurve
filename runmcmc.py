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
nwalkers = 24
nsteps = 500

# set up the walkers in a "Gaussian ball" around the literature estimate for distance to Cepheus cloud (distance mod=10)

ls_result=[4.5,4.6,5,5.1,6.2,7.0,7.75,8.0,9,12,14,16]
starting_positions = [ls_result + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

#starting_positions = [ls_result + 1e-4*np.random.randn(ndim)]

# set up the sampler object
sampler = emcee.EnsembleSampler(nwalkers, ndim, model.log_posterior, args=('simulated_data.h5', 'pixel0000'))

#cov=np.eye(12)
#sampler = emcee.MHSampler(cov,ndim, model.log_posterior, args=('simulated_data.h5', 'pixel0000'))

                                
# run the sampler. We use iPython's %time directive to tell us 
# how long it took (in a script, you would leave out "%time")
sampler.run_mcmc(starting_positions, nsteps)

print('Done')

# Plot traces to see how much you need to burn off, and to check for convergence 
fig, (ax_d1, ax_d2, ax_d3, ax_d4, ax_d5, ax_d6, ax_d7, ax_d8, ax_d9, ax_d10, ax_d11, ax_d12) = plt.subplots(12)

ax_d1.set(ylabel='d1')
ax_d2.set(ylabel='d2')
ax_d3.set(ylabel='d3')
ax_d4.set(ylabel='d4')
ax_d5.set(ylabel='d5')
ax_d6.set(ylabel='d6')
ax_d7.set(ylabel='d7')
ax_d8.set(ylabel='d8')
ax_d9.set(ylabel='d8')
ax_d10.set(ylabel='d10')
ax_d11.set(ylabel='d11')
ax_d12.set(ylabel='d12')

for i in range(10):
    sns.tsplot(sampler.chain[i,:,0], ax=ax_d1)
    sns.tsplot(sampler.chain[i,:,1], ax=ax_d2)
    sns.tsplot(sampler.chain[i,:,2], ax=ax_d3)
    sns.tsplot(sampler.chain[i,:,3], ax=ax_d4)
    sns.tsplot(sampler.chain[i,:,4], ax=ax_d5)
    sns.tsplot(sampler.chain[i,:,5], ax=ax_d6)
    sns.tsplot(sampler.chain[i,:,6], ax=ax_d7)
    sns.tsplot(sampler.chain[i,:,7], ax=ax_d8)
    sns.tsplot(sampler.chain[i,:,8], ax=ax_d9)
    sns.tsplot(sampler.chain[i,:,9], ax=ax_d10)
    sns.tsplot(sampler.chain[i,:,10], ax=ax_d11)
    sns.tsplot(sampler.chain[i,:,11], ax=ax_d12)
    
# Burn off initial steps
samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
    
traces = samples.reshape(-1, ndim).T

parameter_samples = pd.DataFrame({'d1': traces[0], 'd2': traces[1], 'd3': traces[2], 'd4': traces[3], 'd5': traces[4], 'd6': traces[5], 'd7': traces[6], 'd8': traces[7], 'd9': traces[8], 'd10': traces[9], 'd11': traces[10], 'd12': traces[11]})

q = parameter_samples.quantile([0.16,0.50,0.84], axis=0)

#for i in range(1,13):
#	print("d%i = {:.2f} + {:.2f} - {:.2f}".format(q['d%i'][0.50], 
#                                            	q['d%i'][0.84]-q['d%i'][0.50],
#                                            	q['d%i'][0.50]-q['d%i'][0.16]) % (i,i,i,i,i))
#
#print("d%i = {:.2f} + {:.2f} - {:.2f}".format(q['d7'][0.50], 
#                                            	q['d7'][0.84]-q['d7'][0.50],
#                                            	q['d7'][0.50]-q['d7'][0.16]))
                                                           
plt.show()
