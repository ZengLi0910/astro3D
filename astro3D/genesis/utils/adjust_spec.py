"""
Authors: Jacob Seiler, Manodeep Sinha
"""
#!/usr/bin:env python
from __future__ import print_function
from astro3D.genesis.utils import common as cmn

import numpy as np
import h5py
from tqdm import tqdm
import time

__all__ = ("adjust_spec", "adjust_hostHaloID", )


def adjust_hostHaloID(f_out, haloID_field, FirstHaloInFOFgroup_field,
                      Snap_Keys, Snap_Nums, index_mult_factor):
    """
    Adjusts the `hostHaloID` field in the output HDF5 file.

    In the original trees, if a halo is the main background FoF halo, its value
    of `hostHaloID` is set to `-1`.  Under the LHaloTree specs, the property
    correpsonding to this field (`FirstHaloInFOFgroup`) can never be -1.
    Instead, halos in these instances should point to themselves. 

    Parameters
    ----------

    f_out : Open HDF5 file. 
        The HDF5 trees we're adjusting.

    haloID_field : String, optional
        Field name within the HDF5 file that corresponds to the unique halo ID.

    FirstHaloInFOFgroup_field : String, optional
        Field name within the HDF5 file that corresponds to
        `FirstHaloInFOFgroup` in the LHaloTree structure.

    Snap_Keys : List of strings.
        Names of the snapshot keys within the passed keys.

    Snap_Nums : Dictionary of integers keyed by `Snap_Keys`.
        Snapshot number of each snapshot key.

    index_mult_factor: Integer, optional
        Multiplication factor to generate a temporally unique halo ID.

    Returns
    ----------

    None.
    """

    for snap_key in tqdm(Snap_Keys):
        original = f_out[snap_key][FirstHaloInFOFgroup_field][:]

        NHalos = len(original)

        if NHalos == 0:
            continue

        # Find those halos with -1 in the field. 
        inds_to_adjust = np.where(original == -1)[0]

        # Since we're using the original trees, need to generate temporally
        # unique IDs. The halos will point to themselves so use the indices
        # themselves. 
        IDs = cmn.index_to_temporalID(inds_to_adjust, Snap_Nums[snap_key],
                                      index_mult_factor) 

        # Then adjust!
        f_out[snap_key][FirstHaloInFOFgroup_field][list(inds_to_adjust)] = IDs 

def adjust_spec(fname_in, fname_out, haloID_field="ID",
                FirstHaloInFOFgroup_field="hostHaloID",
                index_mult_factor=int(1e12)):
    """
    Adjusts some fields of the VELOCIraptor trees to match the LHaloTree Specs.

    Currently calls the following functions:
        1. :py:mod:`astro3D.genesis.utils.`adjust_hostHaloID`

    Parameters
    ----------

    fname_in, fname_out : String
        Path to the input HDF5 trees and path to where the updated trees will be
        saved.

    haloID_field : String, optional
        Field name within the HDF5 file that corresponds to the unique halo ID.

    FirstHaloInFOFgroup_field : String, optional
        Field name within the HDF5 file that corresponds to
        `FirstHaloInFOFgroup` in the LHaloTree structure.

    index_mult_factor: Integer, optional
        Multiplication factor to generate a temporally unique halo ID.

    Returns
    ----------

    None.

    Notes
    ----------

    The default parameters are chosen to match the ASTRO3D Genesis trees as
    produced by VELOCIraptor + Treefrog.
    """

    print("")
    print("=================================")
    print("Running adjust_spec")
    print("Input Trees: {0}".format(fname_in))
    print("Output Trees: {0}".format(fname_out))
    print("HaloID Field: {0}".format(haloID_field))
    print("Index Mult Factor: {0}".format(index_mult_factor))
    print("=================================")
    print("")

    start_time = time.time()

    with h5py.File(fname_in, "r") as f_in, \
         h5py.File(fname_out, "w") as f_out:

        Snap_Keys, Snap_Nums = cmn.get_snapkeys_and_nums(f_in.keys())

        print("Copying the file over")
        for key in tqdm(f_in.keys()):
            # First copy all the groups to the output file.
            cmn.copy_group(f_in, f_out, key)
        print("Done!")

        # Then adjust all the things we need to adjust.
        print("Adjusting the hostHaloIDs to point to themselves rather than -1.")
        adjust_hostHaloID(f_out, haloID_field, FirstHaloInFOFgroup_field,
                          Snap_Keys, Snap_Nums, index_mult_factor)
        print("Done!")
 
    end_time = time.time()
    print("Took {0:3f} seconds".format(end_time - start_time))
