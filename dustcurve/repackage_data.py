#this script repackages stellar posterior data from its current format on Harvard Odyssey, 
#and combines it with CO intensity data from the Dame et al. 2001 spectral cubes
#Originally, posteriors are packaged in healpy nside 1024 pixels, but we want everything in nside 128 pixels
#This file preferably needs to be run in /n/home09/ggreen/BMK/output-savesurfs 

import os
from spectral_cube import SpectralCube
import numpy as np
import healpy as hp
import h5py

greg_nside=1024.0 #nside of the originally packaged data

nslices=12

#read in spectral cube
cube=SpectralCube.read('Cepheus.fits')
cube_data=np.asarray(cube)
dl=0.25 #change in longitude per spectral cube pixel
db=0.25 #change in latitude per spectral cube pixel 

def get_pix_integers(pixstr):
    return int(pixstr.rsplit('-')[1])
    
def get_co_array(l,b):
    il=np.divide(np.subtract(111.125,l),dl)
    ib=np.divide(np.subtract(b,11.875),db)
    co_array=cube_data[:,ib,il]    
    return(co_array)
    
def repackage(index,nside=128):
    nside_ratio=(greg_nside/nside)**2
    
    #convert nside 128 index to the 1024 indices it corresponds to in the BMK files
    #these are the nside 1024 pixels whose contents we need to extract from the BMK h5 files  
    need_indices=np.arange(index*nside_ratio,index*nside_ratio+nside_ratio)
    
    pdf_array=np.empty((0,700,120))
    co_array=np.empty((0,12))
    coord_array=np.empty((0,2))

    #save all the h5 files in the directory to an array
    filelist=np.array([])
    for file in os.listdir("/n/home09/ggreen/BMK/output-savesurfs/"):
        if file.endswith(".h5"):
            filelist=np.append(filelist,file)
    
    #iterate over each file        
    for file in filelist:
        
        #save name of all the pixels in the h5 file to an array
        f=h5py.File("/n/home09/ggreen/BMK/output-savesurfs/"+file)
        allpix=[key for key in f.keys()]
        
        #extract the nside 1024 index from the allpix strings 
        v_get_pix_integers=np.vectorize(get_pix_integers)
        allpix_int=v_get_pix_integers(allpix)
        
        #check if any nside 1024 pixels are located within the nside 128 pixel of interest
        if len(np.intersect1d(allpix_int,need_indices)) > 0:
            dsets=np.intersect1d(allpix_int,need_indices).astype(int)
            
            #extract stellar pdfs from individual pixel datasets
            for i in range(0, len(dsets)):
            
                pix_str='pixel 1024-' + str(dsets[i]) + '/'
                good_stars=f[pix_str+'stellar chains'].attrs['converged'][:]==1
                pdf_array=np.vstack((pdf_array,f[pix_str+'stellar pdfs'][:,:,:][good_stars]))
                
                #open file containing coordinate information on stars and then save all the coordinates
                #fin=h5py.File('/n/home09/ggreen/BMK/input/' + filelist)
                fin=h5py.File("/n/home09/ggreen/BMK/input/" + file)
                l,b=(fin['photometry/' + pix_str]['l'][:][good_stars], fin['photometry/' + pix_str]['b'][:][good_stars])
                coord_array=np.vstack((coord_array, np.vstack((l,b)).T))
                
                #grab the CO intensity values of the stars given their precise l,b coordinates
                for i in range(0,l.shape[0]):
                    co_array=np.vstack((co_array, get_co_array(l[i],b[i])))

                fin.close()
                
        f.close()
                
    nstars=pdf_array[:,0,0].shape[0]
    #write out all the packaged data into a new h5 file    
    print(pdf_array.shape)
                
    fwrite = h5py.File(str(index)+ ".h5", "w")
    pdfdata = fwrite.create_dataset("/stellar_pdfs", (nstars,700,120), dtype='f')
    pdfdata[:,:,:]=pdf_array
    
    print(co_array.shape)
    co_data = fwrite.create_dataset("/co_data", (nstars,nslices), dtype='f')
    co_data[:,:]=co_array

    print(coord_array.shape)
    coord_data = fwrite.create_dataset("/coord_data", (nstars,2), dtype='f')
    coord_data[:,:]=coord_array

    fwrite.close()
                
            
            
        
        
        
        

