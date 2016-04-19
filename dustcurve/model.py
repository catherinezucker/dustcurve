import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from dustcurve import pixclass

def get_line_integral(co_array, post_array, dist_array, coeff_array):
    """
    returns log of the line integral over stellar posterior for individual star 
    
    Parameters:
        co_array: array of CO intensities (ndim=nslices) for an individual star 
        post_array: 700x120 stellar posterior array for an individual star
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
    """
    dbins, redbins=convert_to_bins(co_array, post_array, dist_array, coeff_array)
    probpath, dpath=flatten_prob_path(post_array,dbins,redbins)
    
    psum=0.0
    for i in range(0,120):
        if dpath[0] >= i:
            psum+=np.sum(probpath[0:np.where(dpath==i)[0][0]])
    return psum
        
def convert_to_bins(co_array, post_array, dist_array, coeff_array):
    """
    returns: dbins= an array of bin indices in post_array corresponding to the distances to each velocity slice 
             rbins= an array of bin indices in post_array corresponding to the reddening to each velocity slice 
             
    Parameters:
        co_array: array of CO intensities (ndim=nslices) for an individual star 
        post_array: 700x120 stellar posterior array for an individual star
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
    
    #convert actual cumulative reddening to a reddening bin in the stellar posterior array 
    rmin, dr=(0.0, 0.01)
    redbins=np.divide(np.subtract(red_array,rmin), dr)
    redbins=redbins.astype(int)
    return dbins, redbins
    
def flatten_prob_path(post_array, dbins, redbins):
    """
    returns: 
    probpath: an array of probabilities flattened along the reddening axis 
    dpath: the distance bins corresponding to each reddening "ledge"
             
    Parameters:
        post_array: 700x120 stellar posterior array for an individual star
         dbins: an array of bin indices in post_array corresponding to the distances to each velocity slice 
        rbins: an array of bin indicies in post_array corresponding to the reddening to each velocity slice 
    """
    nslices=12
    probpath=np.array([0])
    dpath=np.array([0])
    #flatten the reddening profile along the reddening axis 
    #store the probability and distance bins corresponding to each reddening "ledge" 
    for i in range(0, nslices-1):
        probpath=np.hstack((probpath, post_array[dbins[i]:dbins[i+1],redbins[i]]))
        #dpath=np.hstack((dpath,np.arange(dbins[i],dbins[i+1]+0.01, 1)))
        dpath.np.ones(len(probpath))
    return probpath.flatten(), dpath.flatten()


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
    for i in range(0,nslices):
        if theta[i] < 4 or theta[i] > 19: 
              return -np.inf 
    #else return 0.0 
    return 0.0

    
def log_likelihood(theta, co_array, post_array, nstars):
    """
    returns log of likelihood for all the stars in a single pixel
    
    Parameters:
        theta: model parameters (specified as a tuple)
        PixStars: Holds all the stellar posteriors and CO emission for all the stars in a single pixel 
    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    coeff_array=np.ones((12))
    dist_array=np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12])
    probpix=np.zeros((nstars))
    for i in range(0,nstars):
        probpix[i]=np.log(get_line_integral(co_array[i,:], post_array[i,:,:], dist_array, coeff_array))
    probpix=np.sum(probpix)
	return(probpix)    

def log_posterior(theta, pixel):
    """
    returns log of posterior probability distribution for all the stars in a single pixel 
    
    Parameters:
        theta: model parameters (specified as a tuple)
        PixStars: Holds all the stellar posteriors and CO emission for all the stars in a single pixel 
    """
    pixObj=pixclass.PixStars('testdata.h5', pixel)
    co_array=np.asarray(pixObj.get_co()[:,:])
    post_array=np.asarray(pixObj.get_p()[:,:,:])
    nstars=pixObj.get_n_stars()
    
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    return log_prior(theta) + log_likelihood(theta, co_array, post_array, nstars)
