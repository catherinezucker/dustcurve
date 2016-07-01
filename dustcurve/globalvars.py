import numpy as np

path = '/n/fink1/sportillo/dustcurve/Data/'

unique_co = {}
unique_post = {}

def setup(dname):
	unique_co[dname] = np.load(path+dname+'unique_co.npy')
	unique_post[dname] = np.load(path+dname+'unique_post.npz')['arr_0']
