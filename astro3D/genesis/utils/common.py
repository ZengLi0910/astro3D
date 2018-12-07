#!/usr/bin/env python

from __future__ import print_function

__all__ = ()


def snap_key_to_snapnum(snap_key):
    """
    Given the name of a snapshot key, finds the associated snapshot number.

    This is necessary because the 0th snapshot key may not be snapshot 000 and
    there could be missing snapshots. This function searches backwards for a
    group of digits that identify the snapshot number.  If there are numbers
    outside of this cluster they will be disregarded and a warning raised.

    For example, if the key is "Snap1_030", the function will return 30 and
    issue a warning that there were digits ignored.

    Parameters
    ----------

    snap_key: String.
        The name of the snapshot key.

    Returns
    ----------

    snapnum: Integer.
        The snapshot number that corresponds to the snapshot key.

    Examples
    ----------

    >>> snap_key_to_snapnum('Snap_018')
    18

    >>> snap_key_to_snapnum('018_Snap')
    18

    >>> snap_key_to_snapnum('Sn3p_018')
    --WARNING--
    For Snapshot key 'Sn3p_018' there were numbers that were not \
clustered together at the end of the key.
    We assume the snapshot number corresponding to this key is 18; \
please check that this is correct.
    18
    """

    snapnum = ""
    reached_numbers = None

    for letter in reversed(snap_key):  # Go backwards through the key.
        if letter.isdigit():
            if reached_numbers == False and len(snapnum):
                print("--WARNING--")
                print("For Snapshot key '{0}' there were numbers that were not"
                      " clustered together at the end of the key.\nWe assume "
                      "the snapshot number corresponding to this key is {1}; "
                      "please check that this is correct."
                      .format(snap_key, int(snapnum[::-1])))
                break 
            # When a number is found, we concatenate it with the others and
            # flag that we have encountered a cluster of numbers.
            snapnum = "{0}{1}".format(snapnum, letter)
            reached_numbers = True 

        else:
            # When we reach something that's not a number, turn flag off.
            reached_numbers = False

    snapnum = snapnum[::-1]  # We searched backwards so flip the string around.

    return int(snapnum)  # Cast as integer before returning.


def index_to_temporalID(index, snapnum, index_mult_factor):
    """
    Takes snapshot-local halo index and converts into temporally unique ID.

    Note: IDs start counting at 1.  So the index 0 gets mapped to an ID of 1.

    Parameters
    ----------

    index: array-like of integers, or integer.
        Array or single value that describes the snapshot-local haloID.

    snapnum: Integer.
        Snapshot that the halo/s are/is located at.

    index_mult_factor: Integer.
        Factor to convert a the snapshot-unique halo index to a temporally
        unique halo ID.

    Returns
    ----------
    
    index: array-like of integers, or integer.
        Array or single value that contains the temporally unique haloID.

    Examples
    ----------

    >>> index_to_temporalID(23, 18, 1e12)
    18000000000024
    """

    temporalID = snapnum*int(index_mult_factor) + index + 1

    return temporalID


def temporalID_to_snapnum(temporalID, index_mult_factor):
    """
    Given a temporalID, return the corresponding snapshot number.

    Parameters
    ----------

    ID: array-like of integers, or integer.
        Array or single value that describes the temporalID/s.

    index_mult_factor: integer.
        Factor to convert to from temporally-unique halo ID to snap-shot unique
        halo index.

    Returns
    ----------

    snapnum: array-like of integers, or integer.
        Array or single value that contains the snapshot number corresponding
        to the temporal ID.

    Examples
    ----------

    >>> temporalID_to_snapnum(-1, 1e12)
    0

    >>> temporalID_to_snapnum(18000000000001, 1e12)
    18

    >>> test_list = [18000000000001, 20000000000050, 134000000000005]
    >>> temporalID_to_snapnum(test_list, 1e12)
    array([ 18,  20, 134])


    >>> import numpy as np
    >>> test_array = np.array([20000000000050, 134000000000005])
    >>> temporalID_to_snapnum(test_array, 1e12)
    array([ 20, 134])
    """

    import numpy as np

    if isinstance(temporalID, list) or isinstance(temporalID, np.ndarray):
        snapnum = ((np.subtract(temporalID,1)) / index_mult_factor).astype(int)
    else:
        snapnum = int((temporalID - 1) / index_mult_factor)

    return snapnum


def get_snapkeys_and_nums(file_keys):
    """
    Gets names of snapshot keys and snapshot numbers.

    We assume that the snapshot data keys are named to include the word
    "snap" (case insensitive). We also assume that the snapshot number
    for each snapshot key will be in a single cluster towards the end
    of the key. If this is not the case we issue a warning showing what
    we believe to be the corresponding snapshot number.

    Parameters
    ----------

    file_keys : Keys.
        Keys from a given file or dataset.

    Returns
    ----------

    Snap_Keys : List of strings.
        Names of the snapshot keys within the passed keys.

    Snap_Nums : Dictionary of integers keyed by `Snap_Keys`.
        Snapshot number of each snapshot key.
    """

    Snap_Keys = [key for key in file_keys if ("SNAP" in key.upper())]
    Snap_Nums = dict()
    for key in Snap_Keys:
        Snap_Nums[key] = snap_key_to_snapnum(key)

    return Snap_Keys, Snap_Nums


def copy_group(file_in, file_out, key):
    """
    Copies HDF5 group into a new HDF5 file (with same data-structure).

    Parameters
    ----------

    file_in, file_out: Open HDF5 files.
        HDF5 files for the data being copied (file_in) and the file the
        data is being copied to (file_out).

    key: String. 
        Name of the HDF5 group being copied.

    Returns
    ----------
    None
        None
    """

    group_path = file_in[key].parent.name  # Name of the group path.
    group_id = file_out.require_group(group_path)  # Create the group.
    name = "{0}".format(key)  # Name the group.
    file_in.copy(name, group_id, name=key)  # Copy over the data.


def get_halos_per_forest(f_in, Snap_Keys, haloID_field="ID",
                         forestID_field="ForestID", is_mpi=0,
                         debug=0):
    """
    Determines the number of halos in each forest.

    .. note::
        The default parameters are chosen to match the ASTRO3D Genesis trees as
        produced by VELOCIraptor + Treefrog.    

    Parameters
    ----------

    f_in: Open HDF5 file. 
        HDF5 file that contains the sorted trees.

    Snap_Keys: List of strings.
        List of keys that correspond to the fields containing the snapshot
        data.

    haloID_field: String. Default: 'ID'.
        Field name within the HDF5 file that corresponds to the unique halo ID.

    forestID_field: String. Default: 'ForestID'.
        Field name within the HDF5 file that corresponds to forest ID.

    is_mpi : Integer 
        Flag to denote whether this function has been called with multiple
        processors. If not 0, turns off the progress bar from ``tqdm``.

    debug : Integer
        Flag to denote whether extra debugging information should be printed to
        ``stdout``.
 
    Returns
    ----------

    NHalos_forest: Dictionary
        Dictionary that contains the number of halos for each Forest. Key is
        the ForestID.
    """

    import numpy as np
    import time
    import os 
    from tqdm import tqdm

    start_time = time.time()

    print("")
    print("Generating the dictionary for the number of halos in each forest")

    NHalos_forest = {}

    # If this function has been called with multiple processors, then don't use
    # the tqdm progress bar.
    if is_mpi:
        snap_key_loop = enumerate(Snap_Keys)
    else:
        snap_key_loop = enumerate(tqdm(Snap_Keys))

    for count, snap_key in snap_key_loop: 
        if len(f_in[snap_key][haloID_field]) == 0:  # Skip empty snapshots.
            continue

        halos_counted = 0
        halo_forestids = f_in[snap_key][forestID_field][:]

        # First get the number of halos in each forest then grab the indices
        # (i.e., the forestID as we start from 0) of the forests that have
        # halos.
        if debug:
            print("{0}\thalo_forestIDs {1}".format(snap_key, halo_forestids))

        forestIDs, halos_in_forest = np.unique(halo_forestids,
                                               return_counts=True)

        for forest_num, forest_id in enumerate(forestIDs):
            this_snap_NHalos = halos_in_forest[forest_num]

            # The first time a forest appears it won't have a corresponding key
            try:
                NHalos_forest[forest_id] += this_snap_NHalos
            except KeyError:
                NHalos_forest[forest_id] = 0 


    return NHalos_forest


def get_halos_per_forest_per_snap(f_in, Snap_Keys, haloID_field="ID",
                                  forestID_field="ForestID", is_mpi=0,
                                  debug=0, forests_to_process=None):
    """
    Determines the number of halos in each forest at each snapshot.

    The resulting Dictionary is nested with the outer-key given by the ForestID
    and the inner-key given by the snapshot field name.

    We also generate the offset for each Forest at each snapshot.  This is
    necessary because whilst Forest 5 may saved first at snapshot 20, it isn't
    necessarily saved first at snapshot 21.

    .. note::
        The default parameters are chosen to match the ASTRO3D Genesis trees as
        produced by VELOCIraptor + Treefrog.    

    Parameters
    ----------

    f_in: Open HDF5 file. 
        HDF5 file that contains the sorted trees.

    Snap_Keys: List of strings.
        List of keys that correspond to the fields containing the snapshot
        data.

    haloID_field: String. Default: 'ID'.
        Field name within the HDF5 file that corresponds to the unique halo ID.

    forestID_field: String. Default: 'ForestID'.
        Field name within the HDF5 file that corresponds to forest ID.

    is_mpi : Integer 
        Flag to denote whether this function has been called with multiple
        processors. If not 0, turns off the progress bar from ``tqdm``.

    debug : Integer
        Flag to denote whether extra debugging information should be printed to
        ``stdout``.

    forests_to_process : List of integers
        If not ``None``, then only the specified forests will be used.
        Otherwise all forests are used.
 
    Returns
    ----------

    NHalos_forest_per_snap: Nested Dictionary
        Nested dictionary that contains the number of halos for each Forest at
        each snapshot.  Outer-key is the ForestID and inner-key is the snapshot
        key.

    NHalos_forest_per_snap_offset: Nested Dictionary
        Nested dictionary that contains the offset for each Forest at each
        snapshot. Outer-key is the ForestID and inner-key is the snapshot key.
        This is required because whilst the tree is sorted by ForestID, the
        relative position of the tree can change from snapshot to snapshot.
    """

    import numpy as np
    import time
    import os 
    from tqdm import tqdm

    start_time = time.time()

    print("")
    print("Generating the dictionary for the number of halos in each forest "
          "at each snapshot.")

    if forests_to_process is not None:
        print("Doing this for a total of {0} "
              "forests.".format(len(forests_to_process)))

    NHalos_forest_per_snap = {}
    NHalos_forest_per_snap_offset = {}

    # If this function has been called with multiple processors, then don't use
    # the tqdm progress bar.
    if is_mpi:
        snap_key_loop = enumerate(Snap_Keys)
    else:
        snap_key_loop = enumerate(tqdm(Snap_Keys))

    for count, snap_key in snap_key_loop: 
        if not f_in[snap_key][haloID_field]:  # Skip empty snapshots.
            continue

        halos_counted = 0
        halo_forestids = f_in[snap_key][forestID_field][:]

        # First get the number of halos in each forest then grab the indices
        # (i.e., the forestID as we start from 0) of the forests that have
        # halos.
        if debug:
            print("{0}\thalo_forestIDs {1}".format(snap_key, halo_forestids))

        # If we've been given a specific set of forests to use, count how many
        # halos there are at this snapshot for those forests. Otherwise, count
        # for all forests.
        if forests_to_process is not None:
            forestIDs = forests_to_process

            # TODO: This is quite sloppy. Probably a more elegant way...
            NHalos_forest_snap = []
            for forest_id in forestIDs:
                halos_forest_snap = np.where(halo_forestids == forest_id)[0]
                NHalos_forest_snap.append(len(halos_forest_snap))

        else:
            forestIDs, NHalos_forest_snap = np.unique(halo_forestids,
                                                      return_counts=True)


        if debug:
            print("{0}\tforestIDs {1}\tNHalos_forest {2}".format(snap_key,
                                                                 forestIDs,
                                                                 NHalos_forest_snap))

        # Now go through each forest and add its halos to the dictionary.
        # Necessary to do it in a for loop due to dictionary restrictions.
        for forest_num, forest_id in enumerate(forestIDs):
            NHalos = NHalos_forest_snap[forest_num]

            # The first time a forest appears it won't have a corresponding key
            # in the nested dictionary so create it if it's the case.
            try:
                NHalos_forest_per_snap[forest_id][snap_key] = NHalos 
                NHalos_forest_per_snap_offset[forest_id][snap_key] = halos_counted
            except KeyError:
                NHalos_forest_per_snap[forest_id] = {snap_key: NHalos}
                NHalos_forest_per_snap_offset[forest_id] = {snap_key: halos_counted}

            halos_counted += NHalos

        if debug:
            print("{0}\tNHalos_forest {1}".format(snap_key,
                                                  NHalos_forest_per_snap))

    end_time = time.time()
    print("Creation of number of halos per forest took {0:3f} seconds."
          .format(end_time - start_time))
    print("")

    return NHalos_forest_per_snap, NHalos_forest_per_snap_offset


def search_dict_of_lists(value, dictionary):
    """
    Search through a dictionary of lists for a given value.

    Parameters
    ----------

    value: Any data-type.
        The value that we are searching for in the lists.

    dictionary: Dictionary of lists.
        A dictionary of lists we're searching through. 
 
    Returns
    ----------

    True
        If the value is in the dictionary.

    False
        Otherwise.

    Examples
    ----------

    >>> my_dict = {'People' : ['John', 'Mary', 'Joseph'],
    ...            'Age'    : [21, 8, 87],
    ...            'Height' : [186.4, 203.1, 87.8]}
    >>> search_dict_of_lists("John", my_dict)
    True

    >>> search_dict_of_lists("Carol", my_dict)
    False

    >>> search_dict_of_lists(87, my_dict)
    True

    >>> search_dict_of_lists(5, my_dict)
    False

    >>> search_dict_of_lists(186.4, my_dict)
    True

    >>> search_dict_of_lists(186.9, my_dict)
    False
    """
    
    for key in dictionary.keys():
        if value in dictionary[key]:
            return True

    return False


if __name__ == "__main__":
    import doctest
    import numpy as np
    doctest.testmod()
