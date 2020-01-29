import os
import sys
import numpy as np

import h5py


DataFolder = f'/Users/stephenwilkins/research/simulations/FLARES/data/'

data = h5py.File(f'{DataFolder}/flares_noparticles.hdf5', 'r')

# --- walk through HDF5 file

# Sim/Snapshot/

print(data.keys()) # snapshots
print(data['00'].keys()) # redshifts
print(data['00/005_z010p000'].keys()) # properties
print(data['00/005_z010p000/Galaxy'].keys()) # properties


print(data['00/005_z010p000/Galaxy']['SFR'].keys())
print(data['00/005_z010p000/Galaxy']['SFR']['SFR_10'].shape)

print(data['00/005_z010p000/Galaxy']['SubhaloMass'].shape)
print(data['00/005_z010p000/Galaxy']['Mstar_30'].shape)


print(data['00/005_z010p000/Galaxy/BPASS/ChabrierIMF/Luminosity'].keys()) # properties
print(data['00/005_z010p000/Galaxy/BPASS/ChabrierIMF/Luminosity/DustModelI'].keys()) # properties

# def get_name_shape(name, item):
#     shape = ''
#     if hasattr(item, 'value'):
#         shape = item.shape
#     print(name, shape, list(item.attrs.items()))
#
# data.visititems(get_name_shape)
