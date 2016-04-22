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

#This MH Code has been adapted from code snippets found in the PHYS 201 week 9 MCMC notebook (written by Vinny and Tom) and the PHYS 201 week 9 homework solutions (written by Tom and Kevin Shane)
sampler = emcee.MHSampler(np.diagflat(np.ones(ndim)), ndim, model.log_posterior, args=('simulated_data.h5', 'pixel0000'))

allsamples = np.empty((0,ndim))
pos_array=[np.random.randint(4,19) for i in range(ndim)]
std_array=[0.1 for i in range(ndim)]
starting_positions = emcee.utils.sample_ball((pos_array),(std_array),nwalkers)

# run the sampler. 
for i in range(nwalkers):
    sampler.run_mcmc(starting_positions[i], nsteps)
    allsamples = np.vstack((allsamples,sampler.chain[100:,:]))
print('Done')

fig, (ax_d7, ax_d11) = plt.subplots(2)

ax_d7.set(ylabel='d7')
ax_d11.set(ylabel='d11')

sns.tsplot(allsamples[:,6], ax=ax_d7)
sns.tsplot(allsamples[:,10], ax=ax_d11)

parameter_samples = pd.DataFrame({'d7': allsamples[:,6], 'd11': allsamples[:,10]})


print("d7 = {:.2f} + {:.2f} - {:.2f}".format(q['d7'][0.50], 
                                            q['d7'][0.84]-q['d7'][0.50],
                                            q['d7'][0.50]-q['d7'][0.16]))
print("d11 = {:.2f} + {:.2f} - {:.2f}".format(q['d11'][0.50], 
                                            q['d11'][0.84]-q['d11'][0.50],
                                            q['d11'][0.50]-q['d11'][0.16]))

plt.show()



