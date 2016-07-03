#this script repackages stellar posterior data from its current format on Harvard Odyssey, 
#and combines it with CO intensity data from the Dame et al. 2001 spectral cubes
#Originally, posteriors are packaged in healpy nside 1024 pixels, but we want everything in nside 128 pixels (larger area chunks)

import os
from spectral_cube import SpectralCube
import numpy as np
import healpy as hp
import h5py

nslices=6

#read in spectral cube
cube=SpectralCube.read('Cepheus_6slice.fits')
cube_data=np.asarray(cube)
dl=0.125 #change in longitude per spectral cube pixel
db=0.125 #change in latitude per spectral cube pixel 

def get_pix_integers(pixstr):
    result=pixstr.split(' ')[1].rsplit('-')
    return np.array(result)
    
def get_co_array(l,b):
    il=int(np.divide(np.subtract(112.0625,l),dl)) #112.0625 deg=lower left hand corner longitude value
    ib=int(np.divide(np.subtract(b,10.9375),db)) #10.9375 deg=lower left hand corner latitude value
    co_array=cube_data[:,ib,il]    
    return(co_array)
    
def repackage(index,nside=128):
 
    pdf_array=np.empty((0,700,120))
    co_array=np.empty((0,nslices))
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
                
                #grab the CO intensity values of the stars given their precise l,b coordinates
                for k in range(0,l.shape[0]):
                    co_array=np.vstack((co_array, get_co_array(l[k],b[k])))

                fin.close()
                
        f.close()
                
    nstars=pdf_array[:,0,0].shape[0]
    #write out all the packaged data into a new h5 file    
    print(pdf_array.shape)
                
    fwrite = h5py.File("/n/fink1/czucker/Data/"+str(index)+ ".h5", "w")
    pdfdata = fwrite.create_dataset("/stellar_pdfs", (nstars,700,120), dtype='f')
    pdfdata[:,:,:]=pdf_array
    
    print(co_array.shape)
    co_data = fwrite.create_dataset("/co_data", (nstars,nslices), dtype='f')
    co_data[:,:]=co_array

    print(coord_array.shape)
    coord_data = fwrite.create_dataset("/coord_data", (nstars,2), dtype='f')
    coord_data[:,:]=coord_array
    
    fwrite.close()
                
            
            
    
        
        
        

