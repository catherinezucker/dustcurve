from dustcurve import pixclass
import numpy as np
import multiprocessing
import ctypes
import sys
import mmap
from multiprocessing.sharedctypes import Array

nslices=6

def fetch_args(fnames):
    """
    returns: ldata= the likelihood arguments for running emcee's PTSampler ("loglargs")
             pdata= the prior arguments for running emcee's PTSampler ("logpargs")
             
    Parameters:
        fname: filenames of the pixels you want to read in (a string or array of strings)
        the bounds you'd like for your uniform prior on distance and the sigma value you'd like for the log-normal prior on coefficients; format=[lower_dist, upper_dist, sigma]        
        ratio: the gas-to-dust conversion coefficient you'd like to adopt (a float)

    """
    #one healpix pixel
    if isinstance(fnames, str)==True:
        pixObj=pixclass.PixStars('/n/fink1/czucker/Data/' + fnames)
        co_array=np.asarray(pixObj.get_co()[:,:])
        post_array=np.asarray(pixObj.get_p()[:,:,:])
        nstars=pixObj.get_n_stars()

        #find unique co intensity arrays
        sorted_co = co_array[np.lexsort(co_array.T)]
        unique_co=sorted_co[np.concatenate(([True], np.any(sorted_co[1:] != sorted_co[:-1],axis=1)))]

        #find the indices in stacked stellar posterior array corresponding to these intensities
        indices=[]
        for i in range(unique_co.shape[0]):
            unique_stars=np.asarray(np.where(np.equal(unique_co[i], co_array).all(axis=1)))
            indices.append(unique_stars)

        #group posterior arrays with same CO intensity together and store them in a list
        unique_post = []
        for i in range(unique_co.shape[0]):
            unique_post.append(post_array[indices[i],:,:][0])

        print("Total number of stars used in analysis:", nstars)

        #save arrays to files, so they can then be read back in as global variables in the globalvars module
        np.savez("unique_post", unique_post)
        np.save("unique_co", np.array(unique_co))

    #several healpix pixels
    else:
        co_all=np.empty((0,nslices))
        post_all=np.empty((0,700,120))
        nstars_all=0

        for file in fnames:
            pixObj=pixclass.PixStars('/n/fink1/czucker/Data/'+file)
            co_array=np.asarray(pixObj.get_co()[:,:])
            post_array=np.asarray(pixObj.get_p()[:,:,:])
            nstars=pixObj.get_n_stars()

            co_all=np.vstack((co_all,co_array))
            post_all=np.vstack((post_all,post_array))
            nstars_all+=nstars
        
        #find unique co intensity arrays
        sorted_co = co_all[np.lexsort(co_all.T)]
        unique_co=sorted_co[np.concatenate(([True], np.any(sorted_co[1:] != sorted_co[:-1],axis=1)))]
        unique_co=np.array(unique_co)

        #find the indices in stacked stellar posterior array corresponding to these intensities
        indices=[]
        for i in range(unique_co.shape[0]):
            unique_stars=np.asarray(np.where(np.equal(unique_co[i], co_all).all(axis=1)))
            indices.append(unique_stars)

        #group posterior arrays with same CO intensity together and store them in a list
        unique_post = []
        for i in range(unique_co.shape[0]):
            unique_post.append(post_all[indices[i],:,:][0])

        print("Total number of stars used in analysis:", nstars)

        #save arrays to files, so they can then be read back in as global variables in the globalvars module
        np.savez("unique_post", unique_post)
        np.save("unique_co", np.array(unique_co))
        
def get_nstars(fnames):
    """
    returns: total number of stars in all files
    Parameters:
        fnames: filenames of the pixels you want to read in (a string or array of strings)

    """
    if len(fnames)==1:
        pixObj=pixclass.PixStars('/n/fink1/czucker/Data/' + fnames)
        nstars=pixObj.get_n_stars()
        print("Total number of stars used in analysis:", n_stars)
        return nstars

    else:

        nstars_all=0

        for file in fnames:
            pixObj=pixclass.PixStars('/n/fink1/czucker/Data/'+file)
            nstars=pixObj.get_n_stars()
            nstars_all+=nstars

        print("Total number of stars used in analysis:", nstars_all)
        return nstars_all



    
