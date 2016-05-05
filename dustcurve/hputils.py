#code pulled and partially adapted from github.com/ggreen/bayestar/scripts/hputils.py, available under GNU public license

import numpy as np
import healpy as hp


def lb2pix(nside, l, b, nest=True):
    '''
    Convert (l, b) to pixel index.
    '''
    
    theta = np.pi/180. * (90. - b)
    phi = np.pi/180. * l
    
    return hp.pixelfunc.ang2pix(nside, theta, phi, nest=nest)


def pix2lb(nside, ipix, nest=True, use_negative_l=False):
    '''
    Convert pixel index to (l, b).
    '''
    
    theta, phi = hp.pixelfunc.pix2ang(nside, ipix, nest=nest)
    
    l = 180./np.pi * phi
    b = 90. - 180./np.pi * theta
    
    idx = (l > 180.)
    l[idx] = l[idx] - 360.
    
    return l, b

def pix2lb_scalar(nside, ipix, nest=True, use_negative_l=False):
	'''
	Convert pixel index to (l, b).
	
	Takes scalar input (no arrays).
	'''
	
	theta, phi = hp.pixelfunc.pix2ang(nside, ipix, nest=nest)
	
	l = 180./np.pi * phi
	b = 90. - 180./np.pi * theta
	
	if l > 180.:
		l -= 360.
	
	return l, b


def lb_in_bounds(l, b, bounds):
    '''
    Determine whether the given (l, b) coordinates
    are within the provided bounds.
    
    The bounds are in the format:
    
        l_0, l_1, b_0, b_1
    
    l and b can be either floats or numpy arrays. In the
    first case, a boolean will be returned, and in the
    second, a numpy boolean array will be returned.
    '''
    
    l_0 = np.mod(bounds[0], 360.)
    l_1 = np.mod(bounds[1] - l_0, 360.)
    
    l_p = np.mod(l - l_0, 360.)
    
    return (l_p >= 0.) & (l_p <= l_1) & (b >= bounds[2]) & (b <= bounds[3])

def list_indices(nside, bounds):
    '''
    Return all healpy pixels of side nside within a given area
    
    The bounds are in the format:
    
        l_0, l_1, b_0, b_1
    
        '''
    npix=hp.pixelfunc.nside2npix(nside)
    indices=np.arange(0,npix)
    
    inbounds=lb_in_bounds(pix2lb(nside,indices)[0],pix2lb(nside,indices)[1],bounds)
    
    hplist=np.where(inbounds==True)

    return indices[hplist]
