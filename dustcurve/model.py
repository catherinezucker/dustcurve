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
    red_array=red_array.clip(max=7) #if any of the reddening values are > 7 magnitudes, set to 7, because this is the max value our reddening axis in stellar posterior can hold
    
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
    probpath=np.array([post_array[0:dbins[0],0]]) # first reddening ledge; assume no extinction before first distance bin
    for i in range(0, nslices-1):
        probpath=np.append(probpath, post_array[dbins[i]:dbins[i+1],redbins[i]])
    probpath=np.append(probpath, post_array[dbins[-1]:119,redbins[-1]]) #reddening ledge from last distance bin to end of posterior array
    return probpath.flatten() 

def log_prior(theta):
    """
    returns log of prior probability distribution
    
    Parameters:
        theta: model parameters (specified as a tuple)
    """
    # unpack the model parameters (12 distance parameters for 12 velocity slices)
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta 
    nslices=12    
    #check to make sure each d is within the range specified by prior; if not, return -np.inf
    dcheck=np.array([d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12])
    for i in range(0,nslices):
        if dcheck[i] < 4 or dcheck[i] > 19: 
              return -np.inf 
    print(d7,d11)
    #else return 0 
    return 0.0

def log_likelihood(theta, co_array, post_array, nstars):
    """
    returns log of likelihood for all the stars in a single pixel
    
    Parameters:
        theta: model parameters (specified as a tuple)
        pixel: a string representing the pixel within the hdf5 file we're pulling data from
    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    coeff_array=np.ones((12))
    dist_array=np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12])
    dist_array=np.sort(dist_array, axis=0)
    probpix=np.zeros((nstars))
    for i in range(0,nstars):
        probpix[i]=np.log(get_line_integral(co_array[i,:], post_array[i,:,:], dist_array, coeff_array))
    probpix=np.sum(probpix)
    print(probpix)
    return(probpix)    

def log_posterior(theta,fname,pixel):
    """
    returns log of posterior probability distribution for all the stars in a single pixel 
    
    Parameters:
        theta: model parameters (specified as a tuple)
        pixel: a string representing the pixel within the hdf5 file we're pulling data from
    """
    pixObj=pixclass.PixStars(fname,pixel)
    co_array=np.asarray(pixObj.get_co()[:,:])
    post_array=np.asarray(pixObj.get_p()[:,:,:])
    nstars=pixObj.get_n_stars()
    
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    return log_prior(theta) + log_likelihood(theta, co_array, post_array, nstars)
