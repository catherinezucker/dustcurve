#this script repackages stellar posterior data from its current format on Harvard Odyssey, 
#and combines it with CO intensity data from the Dame et al. 2001 spectral cubes
#Originally, posteriors are packaged in healpy nside 1024 pixels, but we want everything in nside 128 pixels (larger area chunks)

import os
from spectral_cube import SpectralCube
import numpy as np
import healpy as hp
import h5py
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

nslices=6

def get_pix_integers(pixstr):
    result=pixstr.split(' ')[1].rsplit('-')
    return np.array(result)
    
def repackage(index,nside=128):

    #open up Perseus CO and HI data files
    w=WCS("/n/fink1/czucker/Perseus_Coordinates_Final.fits")

    co_data=fits.getdata("/n/fink1/czucker/DAME_slab_Jul21.fits")
    h1_data=fits.getdata("/n/fink1/czucker/GALFA_slab_Jul21.fits")
 
    pdf_array=np.empty((0,700,120))
    co_array=np.empty((0,nslices))
    h1_array=np.empty((0,nslices))
    coord_array=np.empty((0,2))

    #save all the h5 files in the directory to an array
    filelist=np.array([])
    for file in os.listdir("/n/fink1/czucker/bayestar/output/"):
        if file.endswith(".h5"):
            filelist=np.append(filelist,file)

    #iterate over each file        
    for file in filelist:
        
        #save name of all the pixels in the h5 file to an array
        f=h5py.File("/n/fink1/czucker/bayestar/output/"+file)
        allpix=[key for key in f.keys()]

        file_nsides=np.empty((0))
        file_indices=np.empty((0))
        
        #extract the nside values and indices from the file allpix strings 
        for i in range(0,np.array(allpix).shape[0]):
            result=allpix[i].split(' ')[1].rsplit('-')
            file_nsides=np.append(file_nsides,result[0])
            file_indices=np.append(file_indices,result[1])

        unique_nsides=np.unique(file_nsides)
        unique_ratios=np.power(np.divide((unique_nsides.astype(int)),nside),2)

        for i in range(0,len(unique_nsides)):

            #indices of the unique nside we need to find
            need_indices=np.arange(index*unique_ratios[i],index*unique_ratios[i]+unique_ratios[i]).astype(int)

            #indices of the unique nside that we have from the current file
            fetch_indices=np.where(file_nsides==unique_nsides[i])
            have_indices=file_indices[fetch_indices]

            dsets=np.intersect1d(need_indices,have_indices).astype(int)
            #extract info from individual pixel datasets
            for j in range(0, len(dsets)):
                pix_str='pixel '+unique_nsides[i]+'-' + str(dsets[j]) + '/'
                good_stars=f[pix_str+'stellar chains'].attrs['converged'][:]==1
                pdf_array=np.vstack((pdf_array,f[pix_str+'stellar pdfs'][:,:,:][good_stars]))
                
                #open file containing coordinate information on stars and then save all the coordinates
                fin=h5py.File("/n/home09/ggreen/BMK/input/" + file)
                l,b=(fin['photometry/' + pix_str]['l'][:][good_stars], fin['photometry/' + pix_str]['b'][:][good_stars])
                coord_array=np.vstack((coord_array, np.vstack((l,b)).T))
                
                #grab the CO and HI intensity values of the stars given their precise l,b coordinates
                #HI file has been convolved/regridded to CO resolution/grid so the physical coords are same for both
                px,py=w.wcs_world2pix(l,b,1)

                co_array=np.vstack((co_array,co_data[:,py.astype(int),px.astype(int)].T)) #axis order is vel, lat, lon
                h1_array=np.vstack((h1_array,h1_data[:,py.astype(int),px.astype(int)].T))

                fin.close()
                
        f.close()
                

    nstars=pdf_array[:,0,0].shape[0]

    print(nstars,pdf_array.shape,co_array.shape,h1_array.shape,coord_array.shape)

    #write out all the packaged data into a new h5 file                    
    fwrite = h5py.File("/n/fink1/czucker/Data/"+str(index)+ ".h5", "w")

    pdfdata = fwrite.create_dataset("/stellar_pdfs", (nstars,700,120), dtype='f')
    pdfdata[:,:,:]=pdf_array
    
    co_data = fwrite.create_dataset("/co_data", (nstars,nslices), dtype='f')
    co_data[:,:]=co_array

    h1_data = fwrite.create_dataset("/h1_data", (nstars,nslices), dtype='f')
    h1_data[:,:]=h1_array

    coord_data = fwrite.create_dataset("/coord_data", (nstars,2), dtype='f')
    coord_data[:,:]=coord_array

    
    fwrite.close()
                
            
            
    
        
        
        

