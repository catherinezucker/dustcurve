import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


#define global array for calculating the uncertainty on the reddening, which varies depending on the number of slices summed
nslices=12
weight0=np.zeros((12))
weight1=np.zeros((12))
weight2=np.zeros((12))

for i in range(0,nslices):
    #calculate standard deviation of the Gaussian describing the uncertainty in the reddening, which increases as number of slices increases 
    #std adds in quadrature with the std of the convolved Gaussian given by c^2=c1^2+c2^2
    #the uncertainty in a single CO intensity measure (converted to reddening) is 0.005 magnitudes
    std=np.sqrt((i+1)*(0.005)**2)

    #calculate the fraction of the probability that will fall within the center bin, +/- 1 bins
    bin0_frac=norm.cdf(0.005, 0, std)-norm.cdf(-0.005,0,std) #fraction of probability falling within center bin
    bin1_frac=norm.cdf(0.015,0,std)-norm.cdf(-0.015,0,std) #fraction of probability falling within +/-1 bin

    #convert this to the appropriate weighting each bin should get
    weight0[i]=bin0_frac
    weight1[i]=(bin1_frac-bin0_frac)/2.0
    weight2[i]=(1-bin1_frac)/2.0
    

def get_line_integral(stellar_index,co_array,post_array,dist_array,coeff_array,order):
    """
    returns line integral over stellar posterior for individual star 
    
    Parameters:
        stellar_index: The index of the star to retrieve the likelihood for 
        co_array: 2D array of CO intensities for all stars (shape=nstarsxnslices)
        post_array: 3D stellar posterior array containing posteriors for all stars (shape=nstarsx700x120)
        dist_array: array of distances to the slices from MCMC
        coeff_array: array of dust-to-gas coefficients for the slices from MCMC
        order: an array which sorts dist_array and co_array in ascending order
    """
    dbins, redbins=convert_to_bins(co_array[stellar_index,:][order],dist_array, coeff_array)
    prob=flatten_prob_path(post_array[stellar_index,:,:],dbins,redbins)
    
    #return log likelihood for individual star
    return np.log(prob)

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
    psum: probability summed over the flattened path extracted from the reddening profile 
              
    Parameters:
        post_star: 2D stellar posterior array for an individual star (shape=700x120)
        dbins: an array of bin indices in post_array corresponding to the distances to each velocity slice 
        rbins: an array of bin indicies in post_array corresponding to the reddening to each velocity slice 
    """
    nslices=12
    #flatten the reddening profile along the reddening axis 
    #store the probability bins corresponding to each reddening "ledge" 

    psum=0
    psum+=np.sum(post_star[0, 0:dbins[0]])  # sum first reddening ledge; assume no extinction before first distance bin and thus no uncertainty

    for i in range(0,nslices-1):
        #determine the probability on the reddening ledge, +/- 1 bin from the reddening ledge and +/-2 bins from reddening ledge (with appropriate weighting to each)
        center_bin=np.sum(np.multiply(post_star[redbins[i],dbins[i]:dbins[i+1]],weight0[i]))
        pm_one=np.sum(np.multiply(post_star[redbins[i]+1,dbins[i]:dbins[i+1]],weight1[i]))+np.sum(np.multiply(post_star[redbins[i]-1,dbins[i]:dbins[i+1]],weight1[i]))
        pm_two=np.sum(np.multiply(post_star[redbins[i]+2,dbins[i]:dbins[i+1]],weight2[i]))+np.sum(np.multiply(post_star[redbins[i]-2,dbins[i]:dbins[i+1]],weight2[i]))

        psum+=center_bin+pm_one+pm_two

    #for final reddening ledge
    center_bin=np.sum(np.multiply(post_star[redbins[-1],dbins[-1]:119],weight0[-1]))
    pm_one=np.sum(np.multiply(post_star[redbins[-1]+1,dbins[-1]:119],weight1[-1]))+np.sum(np.multiply(post_star[redbins[-1]-1,dbins[i]:119],weight1[-1]))
    pm_two=np.sum(np.multiply(post_star[redbins[-1]+2,dbins[-1]:119],weight2[-1]))+np.sum(np.multiply(post_star[redbins[-1]-2,dbins[i]:119],weight2[-1]))
    
    psum+=center_bin+pm_one+pm_two

    return psum

        
def log_prior(theta,bounds):
    """
    returns: prior for set of distances

    Parameters:
    theta: model parameters
    bounds: an array of the format [lowerbound, upperbound] specifying the prior range on distance
    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    dcheck=np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12])

    #check if any elements are outside specified bounds
    testd = (dcheck<bounds[0]) | (dcheck>bounds[1])
    if testd.any()==True:
        return -np.inf

    return 0

def log_likelihood(theta, co_array, post_array, nstars, ratio):
    """
    returns log of likelihood for all the stars in a single pixel
    
    Parameters:
        theta: model parameters (specified as a tuple)
        co_array: 2D (shape=nstarsxnslices) array of CO intensities for all stars 
        post_array: set of 700x120 stellar posterior arrays all stars (shape=nstarsx700x120)
        nstars: integer storing the number of stars in the file
        ratio: the gas-to-dust coefficient
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
    
    #remove values that aren't a finite number from likelihood array; only sum finite numbers
    if np.any(np.isfinite(prob_ensemble))==True:
        finite=np.where(np.isfinite(prob_ensemble)==True)
        prob_ensemble=prob_ensemble[finite]
        return(np.sum(prob_ensemble))
    #if no finite numbers, return probability of -inf
    else:
        return -np.inf

def log_posterior(theta,co_array,post_array,bounds,n_stars,ratio):
    """
    returns log of posterior probability distribution for all the stars in a single pixel 
    
    Parameters:
        theta: model parameters (specified as a tuple)
        co_array: 2D (shape=nstarsx12) array of CO intensities for all stars 
        post_array: set of 700x120 stellar posterior array all stars (shape=nstarsx700x120)
        n_stars: integer storing the number of stars in the file
        ratio: the gas-to-dust coefficient

    """
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = theta
    return log_prior(theta,bounds) + log_likelihood(theta, co_array, post_array, n_stars, ratio)

