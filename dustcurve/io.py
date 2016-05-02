from dustcurve import pixclass
import numpy as np

def fetch_args(fname, bounds, ratio):
    """
    returns: ldata= the likelihood arguments for running emcee's PTSampler ("loglargs")
             pdata= the prior arguments for running emcee's PTSampler ("logpargs")
             
    Parameters:
        fname: filename of the pixel you want to read in
        the bounds you'd like for your uniform prior
        ratio: the gas-to-dust conversion coefficient you'd like to adopt
    """
    pixObj=pixclass.PixStars(fname)
    co_array=np.asarray(pixObj.get_co()[:,:])
    post_array=np.asarray(pixObj.get_p()[:,:,:])
    n_stars=pixObj.get_n_stars()
    bounds=np.asarray([bounds[0],bounds[1]], dtype='f')
    ratio=ratio
    pdata=(bounds)
    ldata=(co_array,post_array,n_stars,ratio)
    return ldata,pdata
