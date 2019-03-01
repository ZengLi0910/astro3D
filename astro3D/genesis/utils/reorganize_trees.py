"""
Authors: Jacob Seiler, Manodeep Sinha
"""
#!/usr/bin/env python
from __future__ import print_function
from astro3D.genesis.utils import treefrog_to_lhalo as frog

import numpy as np
import tqdm
from tqdm import trange
from tqdm import tqdm

import time

__all__ = ("reorganize_trees", )


def reorganize_trees(fname_in, fname_out, num_files_in, boxsize, N_side=3,
                     binary_flag=1, debug=0):

    # First get the MPI parameters.
    try:
        from mpi4py import MPI
    except ImportError:
        rank = 0
        size = 1
    else:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

    partition_length = boxsize / N_side
 
    for partition_idx in range(rank, pow(N_side, 3), size):

        # The bounds we allow depends on what partition we're working with
        x_bound, y_bound, z_bound= determine_task_bounds(N_side, partition_length,
                                                         partition_idx)

        if debug:
            print("partition_idx {0}\tx {1}\ty {2}\tz {3}".format(partition_idx, x_bound,
                                                                  y_bound, z_bound))

        NTrees, NHalos_task, NHalosPerTree_task = handle_trees(fname_in, num_files_in,
                                                               x_bound, y_bound, z_bound,
                                                               debug=debug)

        task_fname_out = "{0}.{1}".format(fname_out, partition_idx)
        frog.write_header(task_fname_out, NTrees, NHalos_task,
                          NHalosPerTree_task, binary_flag) 

        _, _, _ = handle_trees(fname_in, num_files_in, x_bound, y_bound, z_bound,
                               fname_out=task_fname_out, debug=debug)



def handle_trees(fname_in, num_files_in, x_bound, y_bound, z_bound, fname_out=None,
                 debug=0):

    if fname_out:
        f_out = open(fname_out, "ab")
    else:
        NTrees_total = 0
        NHalos_total = 0
        NHalosPerTree_total = []

    LHalo_Desc, _ = frog.get_LHalo_datastruct()

    for file_in_idx in range(num_files_in):

        print("Reading file {0} on rank {1}".format(file_in_idx, rank))

        fname = "{0}.{1}".format(fname_in, file_in_idx)

        with open(fname, "rb") as f_in:

            # First get header info from the binary file.
            NTrees = np.fromfile(f_in, np.dtype(np.int32), 1)[0]
            NHalos = np.fromfile(f_in, np.dtype(np.int32), 1)[0]
            NHalosPerTree = np.fromfile(f_in,
                                        np.dtype((np.int32, NTrees)), 1)[0]

            for tree_idx in range(NTrees):

                binary_tree = np.fromfile(f_in, LHalo_Desc,
                                          NHalosPerTree[tree_idx])

                x = binary_tree["Posx"][0]
                y = binary_tree["Posy"][0]
                z = binary_tree["Posz"][0]

                if x >= x_bound[0] and x < x_bound[1] and \
                   y >= y_bound[0] and y < y_bound[1] and \
                   z >= z_bound[0] and z < z_bound[1]:

                    if fname_out:
                        binary_tree.tofile(f_out)
                    else:
                        NTrees_total += 1
                        NHalos_total += NHalosPerTree[tree_idx]
                        NHalosPerTree_total.append(NHalosPerTree[tree_idx])

    if debug and not fname_out:
        print("NTrees_total {0}\tNHalos_total {1}".format(NTrees_total,
                                                          NHalos_total))

    if fname_out:
        print("Write to file {0}".format(fname_out))
        close(f_out)

    if fname_out:
        return None, None, None
    else:
        return NTrees_total, NHalos_total, NHalosPerTree_total


def determine_task_bounds(N_side, partition_length, partition_idx):
    
    # Bounds will be a half open interval, [x0, x1) etc.
    x = np.zeros(2, dtype=np.float64)
    y = np.zeros(2, dtype=np.float64)
    z = np.zeros(2, dtype=np.float64)

    # We build up partitions from z->y->x. 

    # z increments every step, resetting after `N_side` iterations.
    z[0] = divmod(partition_idx, N_side)[1]*partition_length
    z[1] = z[0] + partition_length 

    # y increments every `N_side` steps, resetting after `N_side^2` iterations.
    y[0] = divmod(divmod(partition_idx, N_side)[0], N_side)[1]*partition_length
    y[1] = y[0] + partition_length

    # x increments every `N_side^2` iterations, never resetting.
    x[0] = divmod(partition_idx, N_side*N_side)[0]*partition_length
    x[1] = x[0] + partition_length

    return x, y, z


if __name__ == "__main__":

    fname = "/fred/oz070/jseiler/astro3d/jan2019/L500_N2160_take2/lhalo/sub_volume/converted"
    check_trees(fname, 64)
