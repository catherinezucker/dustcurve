import h5py

#adapted from class in https://github.com/bayestar/hdf5io.py
#this class reads in data from a repackaged hdf5 file and stores the stellar posteriors and co data from stars in the given pixel
class PixStars:
	def __init__(self, fname=None, pix=None):
		f = None
		close = False
		
		if fname != None:
			if type(pix) != str:
				raise TypeError("If 'fname' is provided, 'dset' must be "
				                "a string containing the name of the dataset.")
			if type(fname) == h5py._hl.files.File:
				f = fname
			else:
				f = h5py.File(fname, 'r')
				close = True
			dpdfs = f[pix+'/stellar_pdfs']
			dco = f[pix+'/co_data']
		
		if pix == None:
			raise ValueError('A dataset name or object must be provided.')
		
		self.load(dpdfs,dco)
		
		if close:
			f.close()
	
	def load(self, dpdfs, dco):
		self.nImages = dpdfs.shape[0]
		self.p = dpdfs[:,:,:]
		self.co=dco[:,:]
	
	#return the stellar posteriors for all the stars in the pixel
	def get_p(self, imgIdx=None):
		if imgIdx == None:
			return self.p
		else:
			return self.p[imgIdx]
			
	#return the co_data for all the stars in the pixel
	def get_co(self, imgIdx=None):
		if imgIdx == None:
			return self.co
		else:
			return self.co[imgIdx]
	
	#return the number of stars in the pixel
	def get_n_stars(self):
		return self.nImages
