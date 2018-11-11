"""
This file contains an example for running the `convert_indices` and
`write_lhalo_binary` functions within`converter.py`.  Refer to the
documentation of those function for full explanation of each variable.

The default parameters are chosen to match the ASTRO3D Genesis trees as
produced by VELOCIraptor + Treefrog.
"""
#!/usr/bin:env python
from __future__ import print_function
from astro3D.genesis.utils import treefrog_to_lhalo as frog_to_l

import sys

if __name__ == '__main__':

    fname_in = "/fred/oz070/jseiler/astro3d/nov2018/N1024_lhalo_indices.hdf5"
    #fname_out = "/fred/oz070/jseiler/astro3d/nov2018/N1024_converted"
    fname_out = "/fred/oz070/jseiler/astro3d/nov2018/tmp"
    haloID_field = "ID"
    forestID_field = "ForestID"
    Nforests = None 
    write_binary_flag = 1 
    fname_alist = "/fred/oz070/jseiler/astro3d/nov2018/N1024_alist"
    dry_run = 1
    debug = 0

    frog_to_l.treefrog_to_lhalo(fname_in, fname_out,
                                haloID_field, forestID_field, Nforests,
                                write_binary_flag, fname_alist=fname_alist,
                                dry_run=dry_run, debug=debug)
