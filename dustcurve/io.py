from dustcurve import pixclass
import numpy as np

nslices=6

def fetch_args(fnames, bounds, ratio):
    """
    returns: ldata= the likelihood arguments for running emcee's PTSampler ("loglargs")
             pdata= the prior arguments for running emcee's PTSampler ("logpargs")
             
    Parameters:
        fname: filenames of the pixels you want to read in (a string or array of strings)
        the bounds you'd like for your uniform prior (given in the format [lowerbound, upperbound]); must be confined to within 4-19 in distance modulus
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

        indices=[]
        for i in range(unique_co.shape[0]):
            unique_stars=np.asarray(np.where(np.equal(unique_co[i], co_array).all(axis=1)))
            indices.append(unique_stars)


        #make unique posterior arrays
        unique_post = []
        for i in range(unique_co.shape[0]):
            unique_post.append(post_array[indices[i],:,:][0])

        #package and return data
        ldata=(unique_co,indices,unique_post,ratio)
        pdata=(np.asarray(bounds,dtype='f'))
        print("Total number of stars used in analysis:", nstars)
        return ldata,pdata    

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

        indices=[]
        for i in range(unique_co.shape[0]):
            unique_stars=np.asarray(np.where(np.equal(unique_co[i], co_all).all(axis=1)))
            indices.append(unique_stars)

        unique_post = []
        for i in range(unique_co.shape[0]):
            unique_post.append(post_all[indices[i],:,:][0])

        #package and return data
        ldata=(unique_co,indices,unique_post,ratio)
        pdata=(np.asarray(bounds,dtype='f'))
        print("Total number of stars used in analysis:", nstars_all)
        return ldata,pdata

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



    
