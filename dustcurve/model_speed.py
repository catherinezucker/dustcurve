import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from dustcurve import pixclass
from dustcurve import globalvars as gv

#stellar posterior array properties
dmin, dmu, dr =(4.0,0.0625,0.01)
dbin_num=240

def get_line_integral(num_batch,dbins,redbins,p_b,dname):
    """
    returns line integral over stellar posterior for individual star

    Parameters:
    num_batch: the index of the batch of stars being used 
        unique_co: an array containing all unique sets of CO intensities 
    unique_post: list of length nbatch, containing a 3D array of 700x120 stellar posterior arrays for all stars in that batch (shape=nstarsx700x120)
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
    """

    #extract reddening profile along each array; creates 2D array of size nstarsx120 where nstars is the number of stars in the CO pixel
    flat_unique_post=gv.unique_post[dname][num_batch]

    #return an array of likelihoods for the batch of stars with a unique set of CO intensities; this is "L_good" in the context of mixture model
    unique_likelihoods=np.sum(flat_unique_post[:,np.repeat(np.insert(redbins[num_batch,:],0,0),np.ediff1d(dbins,to_end=dbin_num-dbins[-1],to_begin=dbins[0])),np.arange(0,dbin_num)],axis=1)

    #mixture model: L_total=P_b*L_bad + (1-P_b)*L_good w/P_b=10^-3, L_bad=1/N_E=1/700
    total_likelihoods=flat_unique_post.shape[0]*((p_b)/700.0) + (1.0-p_b)*unique_likelihoods

    #return total log-likelihood for all stars in this batch
    #mixture model: L_total=P_b*L_bad + (1-P_b)*L_good w/P_b=10^-3, L_bad=1/N_E=1/700
    return np.sum(np.log(total_likelihoods))

def log_prior(theta,bounds,nslices):
    """
    returns: prior for set of distances
    Parameters:
    theta: model parameters
    bounds: an array of the format [lowerbound, upperbound,sigma] specifying the bounds on distance and the sigma on the log-normal prior on coefficients
    """

    #flat prior on distance; check that distances (first half of theta array) are within acceptable bounds
    if ((theta[3:nslices+3] < bounds[0]) | (theta[3:nslices+3] > bounds[1])).any()==True:
        return -np.inf

    #check dust screen:
    if (theta[0] < 4) or (theta[0] > np.min(theta[3:nslices+3])) or (theta[1] < 0.00) or (theta[1]>6.99):
        return -np.inf
    
    #check prior proability that source is 'bad' aka P_b:
    if (theta[2] < 0.0) or (theta[2] > bounds[3]):
        return -np.inf

    #log-normal prior on coefficients:
    else:
        return np.sum(-1*(np.log(theta[nslices+3:])**2)/(2*bounds[2]**2))

def log_likelihood(theta,ratio,dname,nslices):
    """
    returns log of likelihood for all the stars in a single pixel

    Parameters:
        theta: model parameters specified as a tuple)
        unique_co: 2D (shape=nbatch x nslices) array of CO intensities for each batch
    indices: list of length nbatch, with each element being an array of star indices belonging to that batch
        unique_post: list of length nbatch, containing a 3D array of 700x120 stellar posterior arrays for all stars in that batch (shape=nstarsx700x120)
        ratio: the gas-to-dust coefficient
    """

    #sort distances in ascending order and sort the coefficient array+co array by this order
    order=np.argsort(theta[3:nslices+3])
    num_batches=np.arange(0,gv.unique_co[dname].shape[0])
    
    ordered_dist_arr=np.hstack((theta[0],theta[3:nslices+3][order]))

    co_arr=gv.unique_co[dname][:,order]
    coeff_arr=theta[nslices+3:][order]*ratio
    ordered_reddening_arr=np.hstack((np.array(np.ones(co_arr.shape[0])*theta[1]).reshape(-1,1),np.multiply(co_arr,coeff_arr[:,np.newaxis].T)))

    #convert actual distance to a distance bin in the stellar posterior array
    dbins=np.array((ordered_dist_arr-dmin)/dmu).astype(int)
    
    #convert co intensities to cumulative reddenings using gas-to-dust coefficients and then clip array according to reddening bounds on stellar posterior arrays
    #also convert cumulative reddenings to cumulative reddening bin indices
    redbins=np.array((np.cumsum(ordered_reddening_arr,axis=1).clip(min=0,max=6.99))/dr).astype(int)
    
    #construct an array holding the likelihood probability for each individual star and return the sum of the log likelihoods for all the stars
    return (np.sum(np.array(v_get_line_integral(num_batches,dbins,redbins,theta[2],dname))))

def log_posterior(theta,co_array,post_array,bounds,ratio,nslices):
    """
    returns log of posterior probability distribution for all the stars in a single pixel

    Parameters:
        theta: model parameters (specified as a tuple)
        co_array: 2D (shape=nstars x nslices) array of CO intensities for all stars
        post_array: set of 700x120 stellar posterior array all stars (shape=nstarsx700x120)
        n_stars: integer storing the number of stars in the file
        ratio: the gas-to-dust coefficient
    """
    return log_prior(theta,bounds,nslices) + log_likelihood(theta,unique_co,unique_post,ratio,nslices)

#vectorize only the stellar index component of get_line_integral; exclude all other components
v_get_line_integral=np.vectorize(get_line_integral,otypes=[np.float])
v_get_line_integral.excluded.add(1)
v_get_line_integral.excluded.add(2)
v_get_line_integral.excluded.add(3)
v_get_line_integral.excluded.add(4)



