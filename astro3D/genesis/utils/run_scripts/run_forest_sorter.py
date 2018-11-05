"""
This file contains an example for running the `forest_sorter` function within
`forest_sorter.py`.  Refer to the documentation of that function for full
explanation of each variable.

The default parameters are chosen to match the ASTRO3D Genesis trees as
produced by VELOCIraptor + Treefrog.
"""
#!/usr/bin:env python
from __future__ import print_function
from astro3D.genesis.utils import forest_sorter as fs

import numpy as np


if __name__ == '__main__':

    fname_in = "/fred/oz009/N1024/unifiedcatalogs/VELOCIraptor.tree.t4.unifiedhalotree.links.snap.hdf.data"
    fname_out = "/fred/oz004/jseiler/astro3d/nov2018/N1024_sorted.hdf5"
    haloID_field = "ID"
    sort_fields = ["ForestID", "hostHaloID", "Mass_200mean"]
    sort_direction = np.array([1, 1, -1])
    ID_fields = ["Head", "Tail", "RootHead", "RootTail", "ID", "hostHaloID"]
    index_mult_factor = 1e12

    fs.forest_sorter(fname_in, fname_out, haloID_field,
                     sort_fields, sort_direction, ID_fields,
                     index_mult_factor)
