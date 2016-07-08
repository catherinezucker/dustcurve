import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from dustcurve import pixclass
from dustcurve import globalvars as gv

#stellar posterior array properties
dmin, dmu=(4.0, 0.125)
rmin, dr=(0.0, 0.01)

def get_line_integral(num_batch,dist_array,co_coeff_array,h1_coeff_array,order,dname):
    """
    returns line integral over stellar posterior for individual star

    Parameters:
	num_batch: the index of the batch of stars being used 
        unique_co: an array containing all unique sets of CO intensities 
	indices: list of length nbatch, with each element being an array of star indices belonging to that batch (UNUSED)
	unique_post: list of length nbatch, containing a 3D array of 700x120 stellar posterior arrays for all stars in that batch (shape=nstarsx700x120)
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
    """
    #convert actual distance to a distance bin in the stellar posterior array
    dbins=((dist_array-dmin)/dmu).astype(int)
    
    #convert co intensities to cumulative reddenings using gas-to-dust coefficients and then clip array according to reddening bounds on stellar posterior arrays
    #convert cumulative reddenings to cumulative reddening bin indices
    #total cumulative reddening is a linear combination of the contributions from the CO and HI gas. 
    reddening=np.cumsum(gv.unique_co[dname][num_batch][order]*co_coeff_array)+np.cumsum(gv.unique_h1[dname][num_batch][order]*h1_coeff_array)
    redbins=((reddening.clip(min=0,max=6.99)-rmin)/dr).astype(int)

    #extract reddening profile along each array; creates 2D array of size nstarsx120 where nstars is the number of stars in the CO pixel
    flat_unique_post=gv.unique_post[dname][num_batch]

    #return an array of likelihoods for the batch of stars with a unique set of CO intensities; this is "L_good" in the context of mixture model
    unique_likelihoods=np.sum(flat_unique_post[:,np.repeat(np.insert(redbins,0,0),np.ediff1d(dbins,to_end=120-dbins[-1],to_begin=dbins[0])),np.arange(0,120)],axis=1)

    #mixture model: L_total=P_b*L_bad + (1-P_b)*L_good w/P_b=10^-3, L_bad=1/N_E=1/700
    total_likelihoods=np.add(1e-3/700, np.multiply(1-1e-3,unique_likelihoods))

    #return total log-likelihood for all stars in this batch
    #mixture model: L_total=P_b*L_bad + (1-P_b)*L_good w/P_b=10^-3, L_bad=1/N_E=1/700
    return np.sum(np.log(total_likelihoods))

def log_prior(theta,bounds):
    """
    returns: prior for set of distances
    Parameters:
    theta: model parameters
    bounds: an array of the format [lowerbound, upperbound, sigma] specifying the prior range on distance (flat prior)
            and the standard deviation of the log-normal prior on coefficients
    """
    #flat prior on distance
    #check that distances (first half of theta array) are within acceptable bounds
    if ((theta[0:nslices]<bounds[0]) | (theta[0:nslices]>bounds[1])).any()==True:
        return -np.inf
    #log-normal prior on coefficients
    else:
        return np.sum(-1*(np.log(theta[nslices:])**2)/(2*bounds[2]**2))

def log_likelihood(theta,co_ratio,h1_ratio,dname,nslices):
    """
    returns log of likelihood for all the stars in a single pixel

    Parameters:
        theta: model parameters (specified as a tuple)
        co_ratio: the gas-to-dust coefficient for CO
        h1_ratio: the gas to dust coefficient for HI
        nslices: the number of velocity slices used
    """

    dist_array=theta[0:nslices]
    co_coeff_array=np.multiply(theta[nslices:2*nslices],co_ratio)
    h1_coeff_array=np.multiply(theta[2*nslices:],h1_ratio)

    #sort distances in ascending order and sort the coefficient array+co array by this order
    order=np.argsort(dist_array)
    dist_array=dist_array[order]
    co_coeff_array=co_coeff_array[order]
    h1_coeff_array=h1_coeff_array[order]
    num_batches=np.arange(0,gv.unique_co[dname].shape[0])

    #construct an array holding the likelihood probability for each individual star and return the sum of the log likelihoods for all the stars
    return (np.sum(np.array(v_get_line_integral(num_batches,dist_array,co_coeff_array,h1_coeff_array,order,dname))))

def log_posterior(theta,co_array,post_array,bounds,ratio):
    """
    returns log of posterior probability distribution for all the stars in a single pixel

    Parameters:
        theta: model parameters (specified as a tuple)
        co_array: 2D (shape=nstarsx12) array of CO intensities for all stars
        post_array: set of 700x120 stellar posterior array all stars (shape=nstarsx700x120)
        n_stars: integer storing the number of stars in the file
        ratio: the gas-to-dust coefficient
    """
    return log_prior(theta,bounds) + log_likelihood(theta,unique_co,unique_post,ratio)

#vectorize only the stellar index component of get_line_integral; exclude all other components
v_get_line_integral=np.vectorize(get_line_integral,otypes=[np.float])
v_get_line_integral.excluded.add(1)
v_get_line_integral.excluded.add(2)
v_get_line_integral.excluded.add(3)
v_get_line_integral.excluded.add(4)
v_get_line_integral.excluded.add(5)


