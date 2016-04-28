import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from dustcurve import pixclass

def get_line_integral(co_array, post_array, dist_array, coeff_array):
    """
    returns line integral over stellar posterior for individual star 
    
    Parameters:
        co_array: array of CO intensities (ndim=nslices) for an individual star 
        post_array: 700x120 stellar posterior array for an individual star
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
    """
    dbins, redbins=convert_to_bins(co_array, dist_array, coeff_array)
    probpath=flatten_prob_path(post_array,dbins,redbins)
    return np.sum(probpath) 
        
def convert_to_bins(co_array, dist_array, coeff_array):
    """
    returns: dbins= an array of bin indices in post_array corresponding to the distances to each velocity slice 
             rbins= an array of bin indices in post_array corresponding to the reddening to each velocity slice 
             
    Parameters:
        co_array: array of CO intensities (ndim=nslices) for an individual star 
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
    """
    #convert actual distance to a distance bin in the stellar posterior array
    dmin, dmu=(4.0, 0.125)
    dbins=np.divide(np.subtract(dist_array,dmin),dmu)
    dbins=dbins.astype(int)
    
    #convert co intensities to reddenings using gas-to-dust coefficients
    red_array=np.multiply(co_array, coeff_array)
    red_array=np.cumsum(red_array)
    
    #clip reddening values if too high or too low (to fit within bounds of stellar posterior array, with range 0-7 mags
    red_array=red_array.clip(min=0) #if any of the reddening values are negative (due to negative CO intensities, set to zero)
    red_array=red_array.clip(max=6.999) #if any of the reddening values are > 7 magnitudes, set to 7, because this is the max value our reddening axis in stellar posterior can hold
    
    #convert actual cumulative reddening to a reddening bin in the stellar posterior array 
    rmin, dr=(0.0, 0.01)
    redbins=np.divide(np.subtract(red_array,rmin), dr)
    redbins=redbins.astype(int)
    return dbins, redbins
    
def flatten_prob_path(post_array, dbins, redbins):
    """
    returns: 
    probpath: an array of probabilities flattened along the reddening axis, defined by the reddening profile
                 
    Parameters:
        post_array: 700x120 stellar posterior array for an individual star
         dbins: an array of bin indices in post_array corresponding to the distances to each velocity slice 
        rbins: an array of bin indicies in post_array corresponding to the reddening to each velocity slice 
    """
    nslices=12
    #flatten the reddening profile along the reddening axis 
    #store the probability bins corresponding to each reddening "ledge" 
    
    probpath=np.array([post_array[0, 0:dbins[0]]]) # first reddening ledge; assume no extinction before first distance bin
    for i in range(0, nslices-1):
        probpath=np.append(probpath, post_array[redbins[i],dbins[i]:dbins[i+1]])
    probpath=np.append(probpath, post_array[redbins[-1],dbins[-1]:119]) #reddening ledge from last distance bin to end of posterior array
    return probpath.flatten() 

def log_prior(theta):
    """
    returns log of prior probability distribution
    
    Parameters:
        theta: model parameters (specified as a tuple)
    """
    # unpack the model parameters (12 distance parameters for 12 velocity slices)
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12 = theta
    nslices=12    

    dcheck=np.array([d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12])
    ccheck=np.array([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
    
    #check to make sure each d and c is within the range specified by prior; if not, return -np.inf

    if np.any(dcheck) < 4.0 or np.any(dcheck) > 19.0:
        return -np.inf
        
    if np.any(ccheck) < .01 or np.any(ccheck) > 2.0:
        return -np.inf
    
    return 0.0

def log_likelihood(theta, co_array, post_array, nstars):
    """
    returns log of likelihood for all the stars in a single pixel
    
    Parameters:
        theta: model parameters (specified as a tuple)
        pixel: a string representing the pixel within the hdf5 file we're pulling data from
    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12 = theta
    coeff_array=np.array([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
    dist_array=np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12])
    
    #sort in ascending order 
    ascending=np.argsort(dist_array)
    dist_array=dist_array[ascending]
    coeff_array=coeff_array[ascending]
    
    probpix=np.empty((nstars))
    
    for i in range(0,nstars):
        co_array[i,:]=co_array[i,:][ascending] #sort CO slices in ascending order, according to distance estimates
        probpix[i]=np.log(get_line_integral(co_array[i,:], post_array[i,:,:], dist_array, coeff_array))
    probpix=np.sum(probpix)
    return(probpix)    

def log_posterior(theta,co_array,post_array,n_stars):
    """
    returns log of posterior probability distribution for all the stars in a single pixel 
    
    Parameters:
        theta: model parameters (specified as a tuple)
        pixel: a string representing the pixel within the hdf5 file we're pulling data from
    """
    
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12 = theta
    return log_prior(theta) + log_likelihood(theta, co_array, post_array, n_stars)
