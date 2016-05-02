import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from dustcurve import pixclass

def get_line_integral(stellar_index,co_array,post_array,dist_array,coeff_array,order):
    """
    returns line integral over stellar posterior for individual star 
    
    Parameters:
        co_star: 1D array of CO intensities for an individual star (shape=1xnslices)
        post_star: 2D stellar posterior array for an individual star (shape=700x120)
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
    """
    dbins, redbins=convert_to_bins(co_array[stellar_index,:][order],dist_array, coeff_array)
    probpath=flatten_prob_path(post_array[stellar_index,:,:],dbins,redbins)
    
    #return log likelihood for individual star
    return np.log(np.sum(probpath))

def convert_to_bins(co_star, dist_array, coeff_array):
    """
    returns: dbins= an array of bin indices in post_star corresponding to the distances to each velocity slice 
             rbins= an array of bin indices in post_star corresponding to the reddening to each velocity slice 
             
    Parameters:
        co_star: 1D array of CO intensities for an individual star (shape=1xnslices)
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
    """
    #convert actual distance to a distance bin in the stellar posterior array
    dmin, dmu=(4.0, 0.125)
    dbins=np.divide(np.subtract(dist_array,dmin),dmu)
    dbins=dbins.astype(int)

    #convert co intensities to reddenings using gas-to-dust coefficients
    red_array=np.multiply(co_star, coeff_array)
    red_array=np.cumsum(red_array)

    #clip reddening values if too high or too low (to fit within bounds of stellar posterior array, with range 0-7 mags
    red_array=red_array.clip(min=0) #if any of the reddening values are negative (due to negative CO intensities, set to zero)
    red_array=red_array.clip(max=6.999) #if any of the reddening values are > 7 magnitudes, set to 7, because this is the max value our reddening axis in stellar posterior can hold

    #convert actual cumulative reddening to a reddening bin in the stellar posterior array 
    rmin, dr=(0.0, 0.01)
    redbins=np.divide(np.subtract(red_array,rmin), dr)
    redbins=redbins.astype(int)

    return dbins, redbins

def flatten_prob_path(post_star, dbins, redbins):
    """
    returns: 
    probpath: an array of probabilities flattened along the reddening axis, defined by the reddening profile
                 
    Parameters:
        post_star: 2D stellar posterior array for an individual star (shape=700x120)
        dbins: an array of bin indices in post_array corresponding to the distances to each velocity slice 
        rbins: an array of bin indicies in post_array corresponding to the reddening to each velocity slice 
    """
    nslices=12
    #flatten the reddening profile along the reddening axis 
    #store the probability bins corresponding to each reddening "ledge" 

    probpath=np.array([post_star[0, 0:dbins[0]]]) # add first reddening ledge; assume no extinction before first distance bin
    for i in range(0, nslices-1):
        probpath=np.append(probpath, post_star[redbins[i],dbins[i]:dbins[i+1]])
    probpath=np.append(probpath, post_star[redbins[-1],dbins[-1]:119]) #add reddening ledge from last distance bin to end of posterior array
    return probpath.flatten()

def log_prior(theta,bounds):
    """
    returns: prior for set of distances and coefficients

    Parameters:
    theta: model parameters
    lowerd: lower bound on distance
    upperd: upper bound on distance
    lowerc: lower bound on gas-to-dust conversion coefficients
    upperd: upper bound on gas-to-dust conversion coefficients
    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    dcheck=np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12])

   # order=np.argsort(dcheck)
   # if np.array_equal(dcheck,dcheck[order])==False:
   #     return -np.inf

    testd = (dcheck<bounds[0]) | (dcheck>bounds[1])
    if testd.any()==True:
        return -np.inf

    return 0

def log_likelihood(theta, co_array, post_array, nstars, ratio):
    """
    returns log of likelihood for all the stars in a single pixel
    
    Parameters:
        theta: model parameters (specified as a tuple)
        co_array: 2D (shape=nstarsx12) array of CO intensities for all stars 
        post_array: set of 700x120 stellar posterior arrays all stars (shape=nstarsx700x120)
        nstars: integer storing the number of stars in the file
    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta

    dist_array=np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12])
    
    coeff_array=np.ones((12))*ratio

    #sort distances in ascending order and sort the coefficient array by this order
    order=np.argsort(dist_array)
    dist_array=dist_array[order]
    coeff_array=coeff_array[order]

    v_get_line_integral=np.vectorize(get_line_integral,otypes=[np.float])

    v_get_line_integral.excluded.add(1)
    v_get_line_integral.excluded.add(2)
    v_get_line_integral.excluded.add(3)
    v_get_line_integral.excluded.add(4)
    v_get_line_integral.excluded.add(5)

    stellar_indices=np.arange(0,nstars,dtype='i')
    
    #construct an array holding the likelihood probability for each individual star 
    prob_ensemble=np.array(v_get_line_integral(stellar_indices,co_array,post_array,dist_array,coeff_array,order))
    
    if np.any(np.isfinite(prob_ensemble))==True:
        finite=np.where(np.isfinite(prob_ensemble)==True)
        prob_ensemble=prob_ensemble[finite]
        return(np.sum(prob_ensemble))
    else:
        return -np.inf

def log_posterior(theta,co_array,post_array,bounds,n_stars,ratio):
    """
    returns log of posterior probability distribution for all the stars in a single pixel 
    
    Parameters:
        theta: model parameters (specified as a tuple)
        co_array: 2D (shape=nstarsx12) array of CO intensities for all stars 
        post_array: set of 700x120 stellar posterior array all stars (shape=nstarsx700x120)
        nstars: integer storing the number of stars in the file

    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    #print('returned',log_prior(theta,bounds) + log_likelihood(theta, co_array, post_array, n_stars, ratio))
    return log_prior(theta,bounds) + log_likelihood(theta, co_array, post_array, n_stars, ratio)

