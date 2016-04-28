import emcee
from dustcurve import model
import seaborn as sns
import numpy as np
from dustcurve import pixclass

# the model has 12 parameters; we'll use 50 walkers and 500 steps each
ndim = 24
nwalkers = 24
nsteps = 500

fname='simulated_data.h5'
pixObj=pixclass.PixStars(fname)
co_array=np.asarray(pixObj.get_co()[:,:])
post_array=np.asarray(pixObj.get_p()[:,:,:])
nstars=pixObj.get_n_stars()

data=(co_array,post_array,nstars)

#This MH Code has been adapted from code snippets found in the PHYS 201 week 9 MCMC notebook (written by Vinny and Tom) 
#and the PHYS 201 week 9 homework solutions (written by Tom and Kevin Shane)
sampler = emcee.MHSampler(np.diagflat(np.ones(ndim)), ndim, model.log_posterior, args=(data))

allsamples = np.empty((1,ndim))
pos_array=[np.random.randint(4,19) for i in range(ndim)]
std_array=[1. for i in range(ndim)]
starting_positions = emcee.utils.sample_ball((pos_array),(std_array),nwalkers) #set up the initial position vectors for our walkers

# set up and run the sampler 50 different times, and create array of chains
sampler.run_mcmc(starting_positions[0],nsteps)
for i in range(nwalkers):
        sampler.run_mcmc(sampler.chain[-1,:],nsteps)
print('Done')

fig, (ax_d7, ax_d11) = plt.subplots(2)
ax_d7.set(ylabel='d7')
ax_d11.set(ylabel='d11')

sns.tsplot(sampler.chain[1000:,6], ax=ax_d7)
sns.tsplot(sampler.chain[1000:,10], ax=ax_d11)
    
parameter_samples = pd.DataFrame({'d7': sampler.chain[1000:,6], 'd11': sampler.chain[1000:,10]})

q = parameter_samples.quantile([0.16,0.50,0.84], axis=0)

#what values do we get?
print("d7 = {:.2f} + {:.2f} - {:.2f}".format(q['d7'][0.50], 
                                            q['d7'][0.84]-q['d7'][0.50],
                                            q['d7'][0.50]-q['d7'][0.16]))
print("d11 = {:.2f} + {:.2f} - {:.2f}".format(q['d11'][0.50], 
                                            q['d11'][0.84]-q['d11'][0.50],
                                            q['d11'][0.50]-q['d11'][0.16]))
    