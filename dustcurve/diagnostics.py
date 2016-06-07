import numpy as np
import emcee
from dustcurve import model
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from dustcurve import io
import warnings
import line_profiler
from dustcurve import hputils
from dustcurve import kdist

def gelman_rubin(chain_ensemble):
    """
    returns: Gelman-Rubin scale reduction factor for set of chains
    Parameters:
    chain_ensemble: the ensemble of chains to use; this paramater expects that the first nsteps/2 of each chain have already been discarded;
    chain_ensemble should be an array of the format nchains x nsteps x ndim
    
    for more info on GR diagnostic see: http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf

    code adapted from function github.com/hmc/convergence.py available under GNU open source license
    """

    nruns,nsteps,ndim=chain_ensemble.shape

    #calculate the mean of each chain
    mean=np.mean(chain_ensemble,axis=2)
        
    #calculate the variance of each chain
    var=np.var(chain_ensemble,axis=2)
    
    #calculate the mean of the variances of each chain
    W=np.mean(var, axis=0)

    #calculate the variance of the chain means multiplied by n
    B=nsteps*np.var(mean, axis=0)

    #calculate estimated variance:
    sigmasq=(1-1/nsteps) * W + (1/nsteps)*B

    #sampling variability
    V=sigmasq+B/(nruns*nsteps)

    #calculate the potential scale reduction factor
    R=V/W
    
    return R

def run_chains(fnames, nwalkers=100, nsteps=1000, ntemps=5, bounds=[4,19,0,10], runs=5, ratio=0.06, simulated=False):

    from dustcurve import kdist

    #set up sampler
    nwalkers=nwalkers
    nsteps=nsteps
    ntemps=ntemps
    bounds=bounds
    ratio=ratio

    ndim=24
    nslices=12

    #fetch the required likelihood and prior arguments for PTSampler
    ldata,pdata=io.fetch_args(fnames,bounds,ratio)

    #setting off the walkers at the kinematic distance given by the literature, assuming a flat rotation curve, theta=220 km/s, R=8.5 kpc
    #Details on rotation curve given in Rosolowsky and Leroy 2006
    if simulated==False:
        vslices=np.linspace(-15.6,-1.3,12)
        klong=np.ones(12)*hputils.pix2lb_scalar(256,int(fnames[:-3]))[0]
        klat=np.ones(12)*hputils.pix2lb_scalar(256,int(fnames[:-3]))[1]
        kdist=kdist.kdist(klong,klat,vslices)
        kdistmod=5*np.log10(kdist)-5
        result=kdistmod.tolist()
        result.extend(1.0 for i in range (nslices))
    else:
        result=[np.random.randint(4,19) for i in range(nslices)]
        result[6]=np.random.uniform(8,11)
        result[10]=np.random.uniform(14,17)
        result.extend(1.0 for i in range(nslices))

    chain_ensemble=[]
                                 
    for j in range(0,runs):
        #run for independent chains for Gelman-Rubin, with starting positions perturbed 
        #slightly perturb the starting positions for each walker, in a ball centered around kinematic distances and literature gas-to-dust coefficient
        #perturb according to a Gaussian of mean 0 and variance 1
        starting_positions = [[result + np.random.randn(ndim) for i in range(nwalkers)] for j in range(ntemps)]

        #set up the sampler object
        sampler = emcee.PTSampler(ntemps, nwalkers, ndim, model.log_likelihood, model.log_prior, loglargs=(ldata), logpargs=[pdata])

        #burn in, and save final positions for all parameters, which we'll then set off our walkers at for the "real" thing
        post_burn_pos, prob, state = sampler.run_mcmc(starting_positions, 300)
        sampler.reset()

        #run the sampler
        sampler.run_mcmc(post_burn_pos, nsteps)

        #Extract the coldest [beta=1] temperature chain from the sampler object; discard first half of samples as burnin (required for GR diagnostic)
        samples_cold = sampler.chain[0,:,int(nsteps/2):,:]
        traces_cold = samples_cold.reshape(-1, ndim).T

        chain_ensemble.append(traces_cold)
        
    #convert into array of shape nchains x nsteps x ndim
    chain_ensemble=np.array(chain_ensemble)

    gr=gelman_rubin(chain_ensemble)
    return gr, chain_ensemble                                       
                                               


                                           

                                           
                                  
