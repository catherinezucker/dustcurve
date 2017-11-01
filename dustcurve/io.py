from dustcurve import pixclass
from dustcurve import globalvars
import numpy as np
import multiprocessing
import ctypes
import sys
import mmap
from multiprocessing.sharedctypes import Array
from astropy.io import fits

def fetch_args(fnames, dname, nslices,region):
    """
    returns: ldata= the likelihood arguments for running emcee's PTSampler ("loglargs")
             pdata= the prior arguments for running emcee's PTSampler ("logpargs")
             
    Parameters:
        fname: filenames of the pixels you want to read in (a string or array of strings)
        the bounds you'd like for your uniform prior on distance and the sigma value you'd like for the log-normal prior on coefficients; format=[lower_dist, upper_dist, sigma]        
        ratio: the gas-to-dust conversion coefficient you'd like to adopt (a float)

    """

    #noise_threshold=fits.getheader('/n/fink2/czucker/CO_Input/Downsampled_'+region+'.fits')['CDELT3']*0.35/0.06350219

    #one healpix pixel
    if isinstance(fnames, str)==True:
        pixObj=pixclass.PixStars(globalvars.path + fnames)
        co_array=np.asarray(pixObj.get_co()[:,:])
        post_array=np.asarray(pixObj.get_p()[:,:,:])
        coords_array=np.asarray(pixObj.get_coords()[:,:])
        nstars=pixObj.get_n_stars()

        #mask out pixels with negative CO
        #pos_co=np.any(np.isfinite(co_array)==False,axis=1)
        #pos_co=np.array(np.sum(co_array,axis=1) < -999)
        #co_array=co_array[pos_co]
        #post_array=post_array[pos_co]
        #coords_array=coords_array[pos_co]

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

        unique_coords=[]
        for i in range(unique_co.shape[0]):
            unique_coords.append(coords_array[indices[i],:][0])

        print("Total number of stars used in analysis:", nstars)
        #print("Total stars={}, masked stars={}, final star count={}".format(nstars,nstars-co_array.shape[0],co_array.shape[0]))

        #save arrays to files, so they can then be read back in as global variables in the globalvars module
        np.savez(globalvars.path+dname+"_unique_post", unique_post)
        np.save(globalvars.path+dname+"_unique_co", np.array(unique_co))
        np.savez(globalvars.path+dname+"_unique_coords",unique_coords)

    #several healpix pixels
    else:
        co_all=np.empty((0,nslices))
        post_all=np.empty((0,700,240))
        coords_all=np.empty((0,2))
        nstars_all=0

        for file in fnames:
            pixObj=pixclass.PixStars(globalvars.path+file)
            co_array=np.asarray(pixObj.get_co()[:,:])
            post_array=np.asarray(pixObj.get_p()[:,:,:])
            coords_array=np.asarray(pixObj.get_coords()[:,:])
            nstars=pixObj.get_n_stars()

            #pos_co=~np.any(co_array<noise_threshold,axis=1)
            #pos_co=~np.array(np.sum(co_array,axis=1) < -999)
            #co_array=co_array[pos_co]
            #post_array=post_array[pos_co]
            #coords_array=coords_array[pos_co]

            co_all=np.vstack((co_all,co_array))
            coords_all=np.vstack((coords_all,coords_array))
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

       #group posterior arrays with same CO intensity together and store them in a list
        unique_coords = []
        for i in range(unique_co.shape[0]):
            unique_coords.append(coords_all[indices[i],:][0])

        print("Total number of stars used in analysis:", nstars_all)
        #print("Total stars={}, masked stars={}, final star count={}".format(nstars_all,nstars_all-co_all.shape[0],co_all.shape[0]))

        #save arrays to files, so they can then be read back in as global variables in the globalvars module
        np.savez(globalvars.path+dname+"_unique_post", unique_post)
        np.save(globalvars.path+dname+"_unique_co", np.array(unique_co))
        np.savez(globalvars.path+dname+"_unique_coords",unique_coords)
        
def get_nstars(fnames):
    """
    returns: total number of stars in all files
    Parameters:
        fnames: filenames of the pixels you want to read in (a string or array of strings)

    """
    if len(fnames)==1:
        pixObj=pixclass.PixStars(globalvars.path + fnames)
        nstars=pixObj.get_n_stars()
        print("Total number of stars used in analysis:", n_stars)
        return nstars

    else:

        nstars_all=0

        for file in fnames:
            pixObj=pixclass.PixStars(globalvars.path+file)
            nstars=pixObj.get_n_stars()
            nstars_all+=nstars

        print("Total number of stars used in analysis:", nstars_all)
        return nstars_all



    
