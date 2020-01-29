import os
import sys
import numpy as np

import h5py


DataFolder = f'../data'


old = h5py.File(f'{DataFolder}/flares.hdf5', 'r')

new = h5py.File(f'{DataFolder}/flares_noparticles.hdf5', 'w') # --- create new


print(old['00'].keys())

for sim in old.keys():
    for snap in old[sim].keys():

        print(old[f'{sim}/{snap}'].keys())
        new.create_group(f'{sim}/{snap}')
        old.copy(f'{sim}/{snap}/Galaxy', new[f'{sim}/{snap}'])
