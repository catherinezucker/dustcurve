import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from dustcurve import pixclass

#stellar posterior array properties
dmin, dmu=(4.0, 0.125)
rmin, dr=(0.0, 0.01)

def get_line_integral(num_batch,unique_co,indices,unique_post,dist_array,coeff_array):
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
    dbins=np.divide(np.subtract(dist_array,dmin),dmu).astype(int)

    #convert co intensities to cumulative reddenings using gas-to-dust coefficients and then clip array according to reddening bounds on stellar posterior arrays
    #convert cumulative reddenings to cumulative reddening bin indices
    redbins=np.divide(np.subtract(np.cumsum(np.multiply(unique_co[num_batch],coeff_array)).clip(min=0,max=6.99),rmin), dr).astype(int)

    #extract reddening profile along each array; creates 2D array of size nstarsx120 where nstars is the number of stars in the CO pixel
    unique_post=unique_post[num_batch]

    #return an array of likelihoods for the batch of stars with a unique set of CO intensities; this is "L_good" in the context of mixture model
    unique_likelihoods=np.sum(unique_post[:,np.repeat(np.insert(redbins,0,0),np.ediff1d(dbins,to_end=120-dbins[-1],to_begin=dbins[0])),np.arange(0,120)],axis=1)

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
    bounds: an array of the format [lowerbound, upperbound] specifying the prior range on distance

    """
    #check that distances (first half of theta array) are within acceptable bounds
    if ((theta[0:int(len(theta)/2)]<bounds[0]) | (theta[0:int(len(theta)/2)]>bounds[1])).any()==True:
        return -np.inf

    #check that coefficients (second half of theta array) are within acceptable bounds
    if ((theta[int(len(theta)/2):]<bounds[2]) | (theta[int(len(theta)/2):]>bounds[3])).any()==True:
        return -np.inf

    return 0

def log_likelihood(theta,unique_co,indices,unique_post,ratio):
    """
    returns log of likelihood for all the stars in a single pixel

    Parameters:
        theta: model parameters (specified as a tuple)
        unique_co: 2D (shape=nbatch x nslices) array of CO intensities for each batch
	indices: list of length nbatch, with each element being an array of star indices belonging to that batch
        unique_post: list of length nbatch, containing a 3D array of 700x120 stellar posterior arrays for all stars in that batch (shape=nstarsx700x120)
        ratio: the gas-to-dust coefficient
    """
    dist_array=theta[0:int(len(theta)/2)]
    coeff_array=np.multiply(theta[int(len(theta)/2):],ratio)

    #sort distances in ascending order and sort the coefficient array+co array by this order
    order=np.argsort(dist_array)
    dist_array=dist_array[order]
    coeff_array=coeff_array[order]
    unique_co=unique_co[:,order]
    num_batches=np.arange(0,unique_co.shape[0])

    #construct an array holding the likelihood probability for each individual star and return the sum of the log likelihoods for all the stars
    return (np.sum(np.array(v_get_line_integral(num_batches,unique_co,indices,unique_post,dist_array,coeff_array))))

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
    return log_prior(theta,bounds) + log_likelihood(theta,unique_co,indices,post_array,ratio)

#vectorize only the stellar index component of get_line_integral; exclude all other components
v_get_line_integral=np.vectorize(get_line_integral,otypes=[np.float])
v_get_line_integral.excluded.add(1)
v_get_line_integral.excluded.add(2)
v_get_line_integral.excluded.add(3)
v_get_line_integral.excluded.add(4)
v_get_line_integral.excluded.add(5)


