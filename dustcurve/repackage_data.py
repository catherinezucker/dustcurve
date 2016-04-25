#this script repackages stellar posterior data from its current format on Harvard Odyssey, 
#and combines it with CO intensity data from the Dame et al. 2001 spectral cubes
#Originally, posteriors are packaged in healpy nside 1024 pixels, but we want everything in nside 128 pixels
#This file preferably needs to be run in /n/home09/ggreen/BMK/output-savesurfs 

import os
from spectral_cube import SpectralCube
import numpy as np
import healpy as hp
import h5py

greg_nside=1024 #nside of the originally packaged data

nslices=12

#read in spectral cube
cube=SpectralCube.read('Cepheus.fits')
cube_data=np.asarray(cube)
dl=0.25 #change in longitude per spectral cube pixel
db=0.25 #change in latitude per spectral cube pixel 

def get_pix_integers(pixstr):
    return int(pixstr.rsplit('-')[1])
    
def convert_nside():
    mult_factor = (nside_cat/nside)**2
    
def get_co_array(l,b):
    il=(111.125-l)/dl
    ib=(b-11.875)/db
    co_array=cube_data[:,ib,il]
    
def repackage(index,nside=128):
    
    nside_ratio=greg_nside/nside
    
    #convert nside 128 index to the 1024 indices it corresponds to in the BMK files
    #these are the nside 1024 pixels whose contents we need to extract from the BMK h5 files  
    needs_indices=np.arange(index*nside_ratio,index*nside_ratio+nside_ratio)

    pdf_array=np.array([])
    co_array=np.array([])
    coord_array=np.array([])

    #save all the h5 files in the directory to an array
    filelist=np.array([])
    for file in os.listdir("/n/home09/ggreen/BMK/output-savesurfs"):
        if file.endswith(".h5"):
            filelist=np.append(filelist,file)
    
    #iterate over each file        
    for f in filelist:
        
        #save name of all the pixels in the h5 file to an array
        f=h5py.File(f)
        allpix=[key for key in f.keys()]
        
        #extract the nside 1024 index from the allpix strings 
        v_get_pix_integers=np.vectorize(get_pix_integers)
        allpix_int=v_get_pix_integers(allpix)
        
        #check if any nside 1024 pixels are located within the nside 128 pixel of interest
        if len(np.intersect(allpix_int,need_indices)) > 0:
            dsets=np.intersect(allpix_int,need_indices)
            
            #extract stellar pdfs from individual pixel datasets
            for i in range(0, len(dsets)):
            
                pix_str='pixel 1024-' + str(dsets[i]) + '/'
                good_stars=f[pix_str+'stellar chains'].attrs['converged'][:]==1
                pdf_array=np.append(pdf_array,f[pix_str+'stellar pdfs'][:,:,:][goodstars])
                
                #open file containing coordinate information on stars and then save all the coordinates
                fin=h5py.File('/n/home09/ggreen/BMK/input/' + filelist)
                dset_coord=(fin['photometry/' + pix_str]['l'][:][good_stars], fin['photometry/' + pix_str]['b'][:][good_stars])
                coord_array=np.append(coord_array, dset_coord)
                
                #grab the CO intensity values of the stars given their precise l,b coordinates
                v_get_co_array=np.vectorize(get_co_array)
                co_array=np.append(co_array, v_get_co_array(dset_coord[:,0],dset_coord[:,1]))
     
                fin.close()
                f.close()
                
    nstars=co_array.shape[0]
    #write out all the packaged data into a new h5 file            
    fwrite = h5py.File(str(index)+ ".h5", "w")
    pdfdata = fwrite.create_dataset("/stellar_pdfs", (nstars,700,120), dtype='f')
    pdfdata[:,:,:]=pdf_array

    co_data = fwrite.create_dataset("/co_data", (nstars,nslices), dtype='f')
    co_data[:,:]=co_array

    coord_data = fwrite.create_dataset("/coord_data", (nstars,2), dtype='f')
    coord_data[:,:]=coord_array

    fwrite.close()
                
            
            
        
        
        
        

