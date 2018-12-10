import time
import numpy as np
import h5py

from tqdm import tqdm

from astro3D.genesis.utils import common as cmn


__all__ = ("convert_indices", )


def convert_indices(fname_in, fname_out,
                    haloID_field="ID", forestID_field="ForestID",
                    ID_fields=["Head", "Tail", "RootHead", "RootTail", "ID",
                               "hostHaloID"], index_mult_factor=int(1e12),
                    debug=0):
    """
    Converts temporally unique tree IDs to ones that are forest-local as
    required by the LHalo Trees format.

    The data-structure of the Treefrog trees is assumed to be HDF5 File ->
    Snapshots -> Halo Properties at each snapshot.

    A new HDF5 file is saved out with the updated IDs.

    .. note::
        We require the input trees to be sorted via the forest ID
        (``forestID_field``) and suggest to also sub-sort on ``hostHaloID`` and
        mass. Sorting can be done using :py:mod:`astro3D.genesis.utils.forest_sorter`.

    Parameters
    ----------

    fname_in, fname_out : String
        Path to the input HDF5 VELOCIraptor + treefrog trees and the path
        where the LHalo correct trees will be saved.

    haloID_field : String, optional
        Field name within the HDF5 file that corresponds to the unique halo ID.

    forestID_field : String, optional
        Field name within the HDF5 file that corresponds to forest ID.

    ID_fields : List of strings, optional
        The HDF5 field names that correspond to properties that use halo IDs.
        As the halo IDs are updated to match the required LHalo Tree format,
        these must also be updated.

    index_mult_factor : Integer, optional
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
    print("Converting Treefrog IDs to LHalo Tree indices.")
    print("Input Trees: {0}".format(fname_in))
    print("Output LHalo ID Trees: {0}".format(fname_out))
    print("Foest ID Field Names: {0}".format(forestID_field))
    print("ID Field Names: {0}".format(ID_fields))
    print("ForestID Field Name: {0}".format(forestID_field))
    print("Index Mult Factor: {0}".format(index_mult_factor))
    print("=================================")
    print("")

    with h5py.File(fname_in, "r") as f_in, \
         h5py.File(fname_out, "w") as f_out:

        Snap_Keys, Snap_Nums = cmn.get_snapkeys_and_nums(f_in.keys())

        NHalos_forest_per_snap, NHalos_forest_per_snap_offset = \
            cmn.get_halos_per_forest_per_snap(f_in, Snap_Keys, haloID_field, forestID_field)

        print("Now creating a dictionary that maps the old, global indices to "
              "ones that are forest-local.")

        start_time = time.time()

        # NHalos_processed is a dictionary, indexed by the forest number, and
        # represents the number of Halos in this forest that have already been
        # processed.
        NHalos_processed = {}
        for forest in NHalos_forest_per_snap:
            NHalos_processed[forest] = 0

        ID_maps = {}

        if debug:
            print("There are a total of {0} "
                  "forests.".format(len(NHalos_forest_per_snap.keys())))

        for snap_key in tqdm(Snap_Keys[::-1]):
            try:
                NHalos = len(f_in[snap_key][haloID_field])
                if NHalos == 0:
                    continue
            except KeyError:
                continue

            oldIDs_global = []
            newIDs_global = []

            forests_thissnap = np.unique(f_in[snap_key][forestID_field][:])

            oldIDs = f_in[snap_key][haloID_field][:]

            if debug:
                print("For {0} the forests are {1}".format(snap_key,
                                                           forests_thissnap))
                print("The oldIDs for the halos in this snapshot are "
                      "{0}".format(oldIDs))

            for forest in forests_thissnap:

                NHalos_snapshot = NHalos_forest_per_snap[forest][snap_key]
                offset = NHalos_forest_per_snap_offset[forest][snap_key]

                NHalos_been_processed = NHalos_processed[forest]

                idx_lower = offset
                idx_upper = NHalos_snapshot + offset

                oldIDs_thisforest = oldIDs[idx_lower:idx_upper]
                newIDs_thisforest = np.arange(NHalos_been_processed,
                                              NHalos_been_processed + NHalos_snapshot)

                for val1, val2 in zip(oldIDs_thisforest, newIDs_thisforest):
                    oldIDs_global.append(int(val1))
                    newIDs_global.append(int(val2))

                NHalos_processed[forest] += NHalos_snapshot

            oldIDs_to_newIDs = dict(zip(list(oldIDs_global),
                                        list(newIDs_global)))
            ID_maps[Snap_Nums[snap_key]] = oldIDs_to_newIDs

        # For some ID fields (e.g., NextProgenitor), the value is -1.
        # When we convert from the temporalID to a snapshot number, we
        # subtract 1 and divide by the multiplication factor (default 1e12)
        # then cast to an integer.  Hence -2 divided by a huge number will
        # be less than 1 and when it's cast to an integer will result in 0.
        # So the 'Snapshot Number' for values of -1 will be 0.  We want to
        # preserve these -1 flags so we map -1 to -1.

        # However we also need to preserve the dictionary for `Snap_000`...

        try:
            oldID_maps_zero_keys = list(ID_maps[0].keys())
            oldID_maps_zero_values = list(ID_maps[0].values())
        except KeyError:
            # Catch for case where no halos at Snapshot 0 (hence ID_maps[0]
            # won't exist).
            ID_maps[0] = {-1: -1}
        else:
            ID_maps[0] = dict(zip(oldID_maps_zero_keys + [-1],
                                  oldID_maps_zero_values + [-1]))

        end_time = time.time()
        print("Creation of dictionary took {0:3f} "
              "seconds.".format(end_time - start_time))

        print("Copying the old tree file to a new one.")
        for key in tqdm(f_in.keys()):
            cmn.copy_group(f_in, f_out, key)

        print("Now going through all the snapshots and updating the IDs.")
        start_time = time.time()

        for snap_key in tqdm(Snap_Keys):
            try:
                NHalos = len(f_in[snap_key][haloID_field])
                if NHalos == 0:
                    continue
            except KeyError:
                continue

            forests_thissnap = np.unique(f_in[snap_key][forestID_field][:])

            for field in ID_fields:  # If this field has an ID...

                oldID = f_in[snap_key][field][:]
                snapnum = cmn.temporalID_to_snapnum(oldID,
                                                    index_mult_factor)

                # We now want to map the oldIDs to the new, forest-local
                # IDs.  However because you can't hash a dictionary with a
                # numpy array, this needs to be done manually in a `for`
                # loop.

                newID = [ID_maps[snap][ID] for snap, ID in zip(snapnum,
                                                               oldID)]
                f_out[snap_key][field][:] = newID

            # Field Loop.
        # Snapshot Loop.

        end_time = time.time()
        print("Updating took {0:.3f} seconds.".format(end_time - start_time))
