"""
"""
#!/usr/bin:env python
from __future__ import print_function
from astro3D.genesis.utils import reorganize_trees as reorg 

import sys

if __name__ == '__main__':

    fname_in = "/fred/oz070/jseiler/astro3d/jan2019/L500_N2160_take2/lhalo/converted"
    fname_out = "/fred/oz070/jseiler/astro3d/jan2019/L500_N2160_take2/lhalo/sub_volume/converted"
    num_files_in = 128
    N_side = 4
    boxsize = 500
    binary_flag = 1
    debug = 1

    reorg.reorganize_trees(fname_in, fname_out, num_files_in, N_side=N_side,
                           boxsize=boxsize, binary_flag=binary_flag,
                           debug=debug)
