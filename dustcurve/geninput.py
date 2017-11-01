from astropy.wcs import WCS
import numpy as np
import pickle
import gzip
import sys
import os
sys.path.append('/n/fink2/czucker/dustmaps')
from dustmaps.sfd import SFDQuery
import dustmaps
from astropy.coordinates import SkyCoord
import h5py
from astropy.io import fits
from dustcurve import hputils
from astropy import wcs
from astropy import units as u
import healpy as hp
import pyregion
import shapely
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
import healpy as hp
from dustcurve import hputils
import os
from dustcurve import io
from shapely import affinity

#os.environ['maskfile']='/n/fink2/czucker/CO_Input/PerA_Extn2MASS_F_Gal_REGRIDDED_mask_FINAL.fits'
os.environ['maskfile']='/n/fink2/czucker/CO_Input/PerA_Extn2MASS_F_Gal_REGRIDDED_mask_PUBTEST.fits'
os.environ['infilepath']='/n/home09/ggreen/BMK/input/'
os.environ['outfilepath']='/n/fink2/czucker/bayestar/output/Perseus/'
os.environ['writepath']="/n/fink2/czucker/Perseus_Input/"

dbin_num=240


#read in master bayestar pickle file containing info on every pixel in BMK map
fname='/n/fink2/czucker/bayestar/output/status_BMK_greg.pickle'
with open(fname, 'rb') as f:
    data = pickle.load(f, encoding='latin1') 

master_bayestar_arr=np.vstack((data['nside'],data['pix_idx'])).T


def repackage(target_nside,new_indices,nslices,cubefile,mdwarfs=False): 

    """
    Extract stellar posterior and CO emission information for stars inside given polygonal boundary. Package all stars inside boundary
    into a new h5 files for input into emcee. Writes out the new repackaged files. 

    Parameters:
    -----------
    target_nside: int
    the new nside for repackaged files. Should be significantly smaller than polygonal area of interest. Default: 1024

    new_indices: list
    A list of the healpix indices at nside target_nside within which you want to pull stars from the BMK files

    nslices: int
    the number of CO slices used in analysis

    region: str
    Name of the Perseus region you are repackaging

    mdwarfs: boolean
    Whether or not you only want to repackage mdwarfs (default=False)

    """
    
    ntot_stars=0
    os.environ['cubefile']="/n/fink2/czucker/CO_Input/"+cubefile

    #define all nsides available in bayestar
    nsides_in_bayestar=[64,128,256,512,1024,2048]

    #retrieve and repackage all stars within each pixel at target nside, write out to file
    for index in new_indices:

        pdf_array=np.empty((0,700,120))
        co_array=np.empty((0,nslices))
        coord_array=np.empty((0,2))
        
        #get vertices defining the boundary of the healpix pixel at new nside
        vertices=hp.boundaries(target_nside,index,nest=True).T
    
        #fetch all possible nside and pixel index combos which could overlap with the new pixel; do this for all nsides used in bayestar
        for nside in nsides_in_bayestar:

            #find all pixel indices of given nside that lie within bounds of new pixel
            pix_at_nside=hp.query_polygon(nside,vertices,inclusive=True,nest=True)
            nside_pix_arr=np.vstack((np.ones((pix_at_nside.shape))*nside,pix_at_nside)).T.astype(int)
                                 
            nrows, ncols = master_bayestar_arr.shape
            dtype={'names':['f{}'.format(i) for i in range(ncols)],'formats':ncols * [master_bayestar_arr.dtype]}

            #find out which of these are actually found in bayestar and store their positions in master array
            overlap = np.in1d(master_bayestar_arr.copy().view(dtype), nside_pix_arr.copy().view(dtype))
            extract_indices=np.where(overlap==True)[0]
                
            unique_infiles,unique_index=np.unique(data['infile_idx'][extract_indices],return_index=True)
            
            for infile_idx in unique_infiles:
                fin=h5py.File(os.environ['infilepath'] + "BeautifulMonkeyKing."+str(infile_idx)+'.h5')
                fout=h5py.File(os.environ['outfilepath'] + "BeautifulMonkeyKing."+str(infile_idx)+'.h5')
            
                dset_pixidx_nside=master_bayestar_arr[extract_indices][np.where(data['infile_idx'][extract_indices]==infile_idx)]
                for dset in dset_pixidx_nside:
                    
                    pix_str='pixel '+str(dset[0])+'-' + str(dset[1]) + '/'
                    convergence=fout[pix_str+'stellar chains'].attrs['converged'][:]
                    lstars,bstars=(fin['photometry/' + pix_str]['l'][:], fin['photometry/' + pix_str]['b'][:])
                        
                    #mask stars that are non-converged or outsides bounds of new target pixel
                    star_mask=(convergence==1) & (hputils.lb2pix(target_nside,lstars,bstars)==index)                    
                    
                    #implement M-dwarf selection criteria from Eddie's 2014 paper
                    if mdwarfs==True:
                        mags=(fin['photometry/' + pix_str]['mag'][:,:])
                        coords = SkyCoord(lstars*u.deg, bstars*u.deg, frame='galactic')
                        sfd = SFDQuery()
                        ebv = sfd(coords)
                        A_g=3.172*ebv
                        A_r=2.271*ebv
                        A_i=1.682*ebv
                        g=mags[:,0]
                        r=mags[:,1]
                        i=mags[:,2]
                        
                        dwarf_mask = ((g-(A_g/(A_g-A_r))*(g-r-1.2)) < 20) & ((r-i-((A_r-A_i)/(A_g-A_r))*(g-r)) > 0)
                        star_mask=star_mask & dwarf_mask
                        
                    co_subset,star_mask_final=get_co_array(lstars,bstars,star_mask,nslices)
                    co_array=np.vstack((co_array, co_subset))
                    pdf_array=np.vstack((pdf_array,fout[pix_str+'stellar pdfs'][:,:,:][star_mask_final]))
                    coord_array=np.vstack((coord_array, np.vstack((lstars[star_mask_final],bstars[star_mask_final])).T))
 
                fin.close()
                fout.close()
            
        #write out all the packaged data into a new h5 file    
        nstars=pdf_array.shape[0]       
        ntot_stars+=nstars
        
        if nstars > 0:

            #interpolate distance axis of pdf array by factor of two:
            interp_post=np.zeros((pdf_array.shape[0],700,dbin_num))
            interp_post[:,:,::2]=pdf_array[:,:,:]
            interp_post[:,:,1:-1:2]=0.5*(pdf_array[:,:,0:-1:1]+pdf_array[:,:,1::1])
            interp_post[:,:,-1]=pdf_array[:,:,-1]
            
            fwrite = h5py.File(os.environ['writepath']+str(index)+ ".h5", "w")
            pdfdata = fwrite.create_dataset("/stellar_pdfs", (nstars,700,dbin_num), dtype='f')
            pdfdata[:,:,:]=interp_post
    
            co_data = fwrite.create_dataset("/co_data", (nstars,nslices), dtype='f')
            co_data[:,:]=co_array

            coord_data = fwrite.create_dataset("/coord_data", (nstars,2), dtype='f')
            coord_data[:,:]=coord_array
        
            print(index,pdfdata.shape,co_data.shape,coord_data.shape)
    
            fwrite.close()
        
    print("NSTARS_TOTAL",ntot_stars)


def get_co_array(wx,wy,star_mask,nslices):
    """
    extract CO information for set of spectral cube slices, given a set of Galactic coordinates
    
    Parameters:
    ----------
    wx: float array
    The Galactic longitudes for the stars

    wy: float array
    The Galactic latitudes for the stars

    star_mask: boolean array
    Indicates whether star is non-converged or outside of pixel bounds, so no need to fetch CO info
   
    nslices: int
    The number of slices to extract CO info for. Must be less than or equal to depth of cubefile (the CO spectral cube)

    """
    cubefile=fits.getdata(os.environ['cubefile'])

    #Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(os.environ['maskfile'])

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)

    cgal = SkyCoord(l=wx*u.degree, b=wy*u.degree, frame='galactic')

    px,py=w.wcs_world2pix(cgal.l.value,cgal.b.value,1)
    
    #check to make sure all the pixels are inside boundary of image
    inbounds=(px > 0) & (px < hdulist[0].data.shape[1]) & (py > 0) & (py < hdulist[0].data.shape[0])
    
    px[~inbounds]=0
    py[~inbounds]=0
        
    px=px.astype(int)
    py=py.astype(int)
    
    mask_high_av=((hdulist[0].data[py,px]!=12) & (np.isnan(hdulist[0].data[py,px])==True))
    
    star_mask=star_mask & inbounds & ~mask_high_av

    if np.sum(star_mask) > 0:
        co_array=cubefile[:,py,px][:,np.where(star_mask==True)[0]]

        return co_array.T,star_mask
    
    else:
        return np.zeros((0,nslices)),star_mask


def region_setup(cubefile,regionfile,savename,target_nside,nslices):
    '''
    Repackage and gather input files for a given region. Setup the region so it can be called from globalvars

    Parameters:
    -----------
    cubefile: the 3D cube fits file with CO data for region of choice

    regionfile: a string indicating the name of the region file within which to repackage data

    target_nside: the healpix nside you want for each of the files

    nslices: the number of velocity slices you want stored in the repackaged input files. Must be same as depth of cubefile

    savename: the name you want the CO/posterior arrays to be named in globalvars

    '''
    target_nside=1024

    rboxes = pyregion.open('/n/fink2/czucker/CO_Input/'+regionfile)
    centerl,centerb,length,height,angle=rboxes[0].coord_list

    boxpoly=shapely.geometry.box(centerl-(length/2.0), centerb-(height/2.0), centerl+(length/2.0), centerb+(height/2.0), ccw=True)
    boxpoly = affinity.rotate(boxpoly,360.0-angle)
    
    bb=boxpoly.bounds
    bboxl=[bb[0],bb[0],bb[2],bb[2]]
    bboxb=[bb[1],bb[3],bb[3],bb[1]]
    
    patch=PolygonPatch(boxpoly)

    vec=hp.pixelfunc.ang2vec(bboxl,bboxb,lonlat=True)
    new_indices=hp.query_polygon(target_nside,vec,inclusive=False,nest=True)

    lpix,bpix=hputils.pix2lb(target_nside*np.ones((new_indices.shape)).astype(int),new_indices)

    inside_poly_mask=[]
    for k in range(0,len(new_indices)):
        inside_poly_mask.append(patch.contains_point([lpix[k],bpix[k]],radius=0))
    
    inside_poly_mask=np.asarray(inside_poly_mask)
    new_indices=new_indices[inside_poly_mask]
        
    repackage(target_nside,new_indices,nslices,cubefile) 
        
    fnames=[str(j)+'.h5' for j in new_indices]
    fnames_final=[]
    for f in fnames:
        if os.path.isfile(os.environ['writepath']+f)==True:
            fnames_final.append(f)
            
    #setup the CO and stellar posterior arrays for input into emcee
    io.fetch_args(fnames_final,savename,nslices,cubefile)

    print("Done")
