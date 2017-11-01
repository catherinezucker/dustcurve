import numpy as np

path = '/n/fink2/czucker/Perseus_Input/'

unique_co = {}
unique_post = {}
unique_coords = {}

def setup(dname):
	unique_co[dname] = np.load(path+dname+'_unique_co.npy')
	unique_coords[dname] = np.load(path+dname+'_unique_coords.npz')['arr_0']
	unique_post[dname] = np.load(path+dname+'_unique_post.npz')['arr_0']
	print("Globalvars setup complete")


