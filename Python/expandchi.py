import h5py
import sys
import os
from optparse import OptionParser
import numpy as np

def findigs(ref,qg_sc):
    nmtx = ref['eps_header/gspace/nmtx_max'][()]
    for igs in range(nmtx):
        gp_e2r  = ref['eps_header/gspace/gind_eps2rho'][0,igs]
        gp_ref  = ref['mf_header/gspace/components'][gp_e2r,:]
        if np.linalg.norm(gp_ref-qg_sc) < 1e-6:
            return igs

usage = '%prog: [options]'
parser = OptionParser(usage=usage)
parser.add_option('-q', '--scqmesh', nargs=2, dest=scqmesh, type='int',
                  help="supercell q mesh, give two arguments.")
parser.add_option('-s', '--scsize', nargs=4, dest=scsize, type='int',
                  help="supercell size, give four arguments in the format
                        n1 n2 m1 m2.
                        [n1 n2
                         m1 m2]")
parser.add_option('--ismetal',default=False, action='store_true',
                  help="Turn on if your system is metallic. BGW handles the
                        number of k-points differently for metals and non-metals.
                        Default False if flag not included on input.")
parser.add_option('--isspin',default=False, action='store_true',
                  help="Turn on if your system is spin-polarized. 
                        Default False if flag not included on input.")
# Implement subsampling later

(options, args) = parser.parse_args()
if len(args)<6:
    parser.error("wrong number of arguments")

# Setup variables
p,q = options.scqmesh[0], options.scqmesh[1]
n1,n2,m1,m2 = options.scsize[0],options.scsize[1],options.scsize[2],options.scsize[3]

mm = n1*m2 - n2*m1
if (mm < 0):
    print("Be careful that the determinant is negative\n")
    mm = -mm
nq, nqs = mm*p*q, p*q
nq0 = 1
print("total number of q points in unit  cell: " + str(nq))
print("total number of q points in super cell: " + str(nqs))

chi0sc  = h5py.File("chi0mat_sc.h5","w")
chisc   = h5py.File( "chimat_sc.h5","w")
chi0uc  = h5py.File("chi0mat_uc.h5","r+")
chiuc   = h5py.File( "chimat_uc.h5","r+")
chi0ref = h5py.File("chi0mat_rf.h5","r+")
chiref  = h5py.File( "chimat_rf.h5","r+")

# Read in chi0uc
print("Copy headers from ref chimat\n")

# Strategy: copy over headers from ref chimat, then initialize blank matrices for the chimat
chi0sc.copy(chi0ref['mf_header'], 'mf_header')
chisc.copy(  chiref['mf_header'], 'mf_header')
chi0sc.copy(chi0ref['eps_header'], 'eps_header')
chisc.copy(  chiref['eps_header'], 'eps_header')

# Initialise headers for chi0 and chi matrices
chi0sc.create_group('mats')
name = 'mats/matrix-diagonal'
shape = np.array(chi0ref[name].shape)
chi0sc.create_dataset(name, shape, chi0ref[name].dtype)
name = 'mats/matrix'
shape = np.array(chi0ref[name].shape)
chi0sc.create_dataset(name, shape, chi0ref[name].dtype)

chisc.create_group('mats')
name = 'mats/matrix-diagonal'
shape = np.array(chiref[name].shape)
chisc.create_dataset(name, shape, chiref[name].dtype)
name = 'mats/matrix'
shape = np.array(chiref[name].shape)
chisc.create_dataset(name, shape, chiref[name].dtype)

print("Copy done\n")

# Do chi0 first. Ideally I should do this for both metal and non-metals...
# Also to deal with the case of subsampling.
print("Dealing with chi0\n")
nmtx_uc = chi0uc['eps_header/gspace/nmtx'][0]
q_uc = chi0uc['/eps_header/qpoints/qpts'][0,:]
qg_sc, qgp_sc= [0.0,0.0,0.0], [0.0,0.0,0.0]

# Deal with the metal case later.
name = 'mats/matrix'
for ig in range(nmtx_uc):
    # q+G in unit cell...
    # Dangerous. epsmat is sorted wrt epsmat g-space.
    # Get g-vec from mf_header first, then map to eps.
    g_e2r = chi0uc['eps_header/gspace/gind_eps2rho'][0,ig]
    g_uc  = chi0uc['mf_header/gspace/components'][g_e2r,:]
    qg_sc[0] = n1 * g_uc[0] + n2 * g_uc[1] + round(n1 * q_uc[0] + m2 * q_uc[1] )
    qg_sc[1] = m1 * g_uc[0] + m2 * g_uc[1] + round(m1 * q_uc[0] + m2 * q_uc[1] )
    qg_sc[2] = g_uc[2]
    igs = findigs(chi0ref,qg_sc) 
    for igp in range(nmtx_uc):
        gp_e2r = chi0sc['eps_header/gspace/gind_eps2rho'][0,igp]
        gp_uc  = chi0sc['mf_header/gspace/components'][gp_e2r,:]
        qgp_sc[0] = n1 * gp_uc[0] + n2 * gp_uc[1] + round(n1 * q_uc[0] + m2 * q_uc[1] )
        qgp_sc[1] = m1 * gp_uc[0] + m2 * gp_uc[1] + round(m1 * q_uc[0] + m2 * q_uc[1] )
        qgp_sc[2] = gp_uc[2]
        igps = findigs(chi0ref,qgp_sc) 
        chi0sc[name][0,:,:,igps,igp,:] = chi0uc[name][0,:,:,igp,ig,:]

print("chi0 done\n")

print("Dealing with chi\n")
for iqs in range(nqs):
    print("Dealing with q=%d/%d"%(iqs+1, nqs))
    q_sc = chiref['/eps_header/qpoints/qpts'][iqs,:]
    qg_sc, qgp_sc= [0.0,0.0,0.0], [0.0,0.0,0.0]
    for iq in range(nq):
        q_uc = chiuc['/eps_header/qpoints/qpts'][iq,:]
        if (np.mod(abs(n1 * q_uc[0] + n2 * q_uc[1] - q_sc[0] ) + 0.0001,1.0) > 0.01) or (np.mod(abs(m1 * q_uc[0] + m2 * q_uc[1] - q_sc[1] ) + 0.0001,1.0) > 0.01):
            continue
        nmtx_uc = chiuc['eps_header/gspace/nmtx'][iq]
        qg_sc, qgp_sc= [0.0,0.0,0.0], [0.0,0.0,0.0]
        
        for ig in range(nmtx_uc):
            g_e2r = chiuc['eps_header/gspace/gind_eps2rho'][iq,ig]
            g_uc  = chiuc['mf_header/gspace/components'][g_e2r,:]
            qg_sc[0] = n1 * g_uc[0] + n2 * g_uc[1] + round(n1 * q_uc[0] + m2 * q_uc[1] - q_sc[0] )
            qg_sc[1] = m1 * g_uc[0] + m2 * g_uc[1] + round(m1 * q_uc[0] + m2 * q_uc[1] - q_sc[1] )
            qg_sc[2] = g_uc[2]
            igs = findigs(chiref,qg_sc) 
            for igp in range(nmtx_uc):
                gp_e2r = chiuc['eps_header/gspace/gind_eps2rho'][iq,igp]
                gp_uc  = chiuc['mf_header/gspace/components'][gp_e2r,:]
                qgp_sc[0] = n1 * gp_uc[0] + n2 * gp_uc[1] + round(n1 * q_uc[0] + m2 * q_uc[1] - q_sc[0] )
                qgp_sc[1] = m1 * gp_uc[0] + m2 * gp_uc[1] + round(m1 * q_uc[0] + m2 * q_uc[1] - q_sc[1] )
                qgp_sc[2] = gp_uc[2]
                igps = findigs(chiref,qgp_sc) 
                name = 'mats/matrix'
                chisc[name][iqs,:,:,igps,igp,:] = chiuc[name][iq,:,:,igp,ig,:]
print("chi done\n")
print("Job Done!\n")
