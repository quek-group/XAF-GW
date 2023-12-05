import h5py
import numpy as np
# Put flag for optional chi0
chi0a = h5py.File("chi0mat_a.h5","r+")
chi0b = h5py.File("chi0mat_b.h5","r+")
chi0tot = h5py.File("chi0mat_tot.h5","w")

chia = h5py.File("chimat_a.h5","r+")
chib = h5py.File("chimat_b.h5","r+")
chitot = h5py.File("chimat_tot.h5","w")

# Check that some of the datasets are the same!
#assert 
#nq = chia[''][()]

# Copy headers, prepare datasets, add 
chi0tot.copy(chi0a['mf_header'],'mf_header')
chi0tot.copy(chi0a['eps_header'],'eps_header')
chi0tot.create_group('mats')
name = 'mats/matrix'
shape = np.array(chi0a[name].shape)
zeros = np.zeros(shape)
chi0tot.create_dataset(name, shape, data=zeros, chi0a[name].dtype)

chitot.copy(chia['mf_header'],'mf_header')
chitot.copy(chia['eps_header'],'eps_header')
chitot.create_group('mats')
name = 'mats/matrix'
shape = np.array(chia[name].shape)
zeros = np.zeros(shape)
chitot.create_dataset(name, shape, data=zeros, chia[name].dtype)
name = 'mats/matrix'
chi0tot[name] = chi0a[name] + chi0b[name]
chitot[name] = chia[name] + chib[name]
