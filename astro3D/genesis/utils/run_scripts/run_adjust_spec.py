"""
This file contains an example for running the `adjust_spec` function within
`adjust_spec.py`.  Refer to the documentation of that function for full
explanation of each variable.

The default parameters are chosen to match the ASTRO3D Genesis trees as
produced by VELOCIraptor + Treefrog.
"""
#!/usr/bin:env python
from __future__ import print_function
from astro3D.genesis.utils import adjust_spec as adjust 

import numpy as np


if __name__ == '__main__':

    fname_in = "/fred/oz009/pelahi/GENESIS/L105_N2048/newtree/VELOCIraptor.tree.t4.0-189.unifiedforest.hdf5"
    fname_out = "/fred/oz070/jseiler/astro3d/nov2018/L105_N2048_hosthaloIDfixed.hdf5"
    haloID_field = "ID"
    FirstHaloInFOFgroup_field = "hostHaloID"
    index_mult_factor = 1e12

    adjust.adjust_spec(fname_in, fname_out, haloID_field,
                       FirstHaloInFOFgroup_field, index_mult_factor)
