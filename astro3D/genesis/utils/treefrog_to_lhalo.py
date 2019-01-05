"""
Authors: Jacob Seiler, Manodeep Sinha
"""
#!/usr/bin/env python
from __future__ import print_function
from astro3D.genesis.utils import common as cmn

import numpy as np
import h5py
import tqdm
from tqdm import tqdm
from tqdm import trange

import time

try:
    from mpi4py import MPI
except ImportError:
    rank = 0
    size = 1
else:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

__all__ = ("treefrog_to_lhalo", )


def get_LHalo_datastruct():
    """
    Generates the LHalo numpy structured array.

    Parameters
    ----------

    None.

    Returns
    ----------

    LHalo_Desc : Numpy structured array
        Structured array for the LHaloTree data format.

    mutltipledim_names : Dictionary of lists
        Specifies the field names for multi-dimensional arrays and the 1D LHalo
        components.

    Notes
    ----------

    Ideally an LHalo tree would specify properties such as position as a single
    Nx3 array. However due to the way ``h5py`` indexes its elements, slicing
    the input data in such a way is not possible. Instead we have used
    individual 1D arrays which when read by a binary reader will correctly
    correspond to an Nx3 array.

    If you specify a HDF5 file to be written, the resulting file will contains
    positions/velocity/spin as single Nx3 arrays.
    """

    LHalo_Desc_full = [
        ('Descendant',          np.int32),
        ('FirstProgenitor',     np.int32),
        ('NextProgenitor',      np.int32),
        ('FirstHaloInFOFgroup', np.int32),
        ('NextHaloInFOFgroup',  np.int32),
        ('Len',                 np.int32),
        ('M_Mean200',           np.float32),
        ('Mvir',                np.float32),
        ('M_TopHat',            np.float32),
        ('Posx',                np.float32),
        ('Posy',                np.float32),
        ('Posz',                np.float32),
        ('Velx',                np.float32),
        ('Vely',                np.float32),
        ('Velz',                np.float32),
        ('VelDisp',             np.float32),
        ('Vmax',                np.float32),
        ('Spinx',               np.float32),
        ('Spiny',               np.float32),
        ('Spinz',               np.float32),
        ('MostBoundID',         np.int64),
        ('SnapNum',             np.int32),
        ('Filenr',              np.int32),
        ('SubHaloIndex',        np.int32),
        ('SubHalfMass',         np.float32)
                         ]

    names = [LHalo_Desc_full[i][0] for i in range(len(LHalo_Desc_full))]
    formats = [LHalo_Desc_full[i][1] for i in range(len(LHalo_Desc_full))]
    LHalo_Desc = np.dtype({'names': names, 'formats': formats}, align=True)

    multipledim_names = {"Pos": ["Posx", "Posy", "Posz"],
                         "Vel": ["Velx", "Vely", "Velz"],
                         "Spin": ["Spinx", "Spiny", "Spinz"]}

    return LHalo_Desc, multipledim_names


def fix_nextprog(forest_halos, forestID, debug=0):
    """
    Walks the descendants of a single forest to generate the `NextProgenitor`
    field.

    Parameters
    ----------

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The halos within a single forest.

    forestID : Integer
        The forest ID for the forest we're updating.

    debug : Integer
        Flag to denote whether extra debugging information should be printed to
        ``stdout``.

    Returns
    ----------

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The forest halos with updated ``NextProgenitor`` field.
    """

    all_descendants = forest_halos["Descendant"][:]
    for ii, d in enumerate(all_descendants):
        if d == ii:
            continue
        curr = forest_halos["FirstProgenitor"][d]

        if debug:
            print("Descendant {0}\tCurr {1}".format(d, curr))

        if curr == ii:
            continue

        try:
            tmp = forest_halos["NextProgenitor"][curr]
        except IndexError:
            print("When attempting to fix the NextProgenitor values for "
                  "forest {0} on Rank {1}, I ran into an "
                  "error.".format(forestID, rank))
            print("However, we only have {0} Halos in this "
                  "forest.".format(len(forest_halos)))
        else:
            pass


        while forest_halos["NextProgenitor"][curr] != -1:

            if debug:
                print("Loop: Descendant {0}\tCurr {1}\tNext {2}".format(d, curr,
                                                                        forest_halos["NextProgenitor"][curr]))

            prev = curr
            curr = forest_halos["NextProgenitor"][curr]

            try:
                tmp = forest_halos["NextProgenitor"][curr]
            except IndexError:
                print("When attempting to fix the NextProgenitor values for "
                      "forest {0} on Rank {1}, I ran into an "
                      "error.".format(forestID, rank))
                print("The NextProgenitor value of Halo {0} was found to be "
                      "{1}".format(prev, curr))
                print("However, we only have {0} Halos in this "
                      "forest.".format(len(forest_halos)))
            else:
                pass

        assert(forest_halos["NextProgenitor"][curr] == -1)
        forest_halos["NextProgenitor"][curr] = ii

    return forest_halos


def fix_flybys(forest_halos, NHalos_root):
    """
    Fixes flybys for a single forest.

    Under the LHalo tree data structure, multiple FoFs at the root redshift are
    allowed IF AND ONLY IF all `FirstHaloInFOFgroup` values point to the same
    FoF. This is not enforced in the Treefrog data structure hence we must fix
    it here.

    We designate the most massive FoF halo at the root redshift to be the
    'True' FoF halo and update the Treefrog-equivalent field of
    ``FirstHaloInFOFgroup`` to point to this most massive halo.

    .. note::
        We pass all halos within the forest to this function but only those at
        the root snapshot are altered.

    Parameters
    ----------

    forest_halos: Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The halos within a single forest.

    NHalos_root: Integer
        The number of halos at the root snapshot.

    Returns
    ----------

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The forest halos with updated ``FirstHaloInFOFgroup`` and
        ``NextHaloInFOFgroup`` fields.

    global_true_fof_idx : Integer
        The forest-local index corresponding to the 'true' FoF halo for this
        forest.

    global_flyby_ind : List of integers.
        The forest-local indices corresponding to those FoF halos that were
        flagged as flybys.
    """

    # Since we're at the root snapshot, the indexing will start from 0.
    fof_halo_inds = np.unique(forest_halos["FirstHaloInFOFgroup"][0:NHalos_root])
    fof_halos = forest_halos[fof_halo_inds]

    # If there is only one FoF Halo, no changes need to be made.
    if len(fof_halo_inds) == 1:
        return forest_halos, None, None

    # Find the 'true' FoF Halo.
    true_fof_idx = np.argmax(fof_halos["Mvir"])
    global_true_fof_idx = fof_halo_inds[true_fof_idx]

    # Find those halos whose `FirstHaloInFOFgroup` is incorrect
    flyby_ind = np.where(fof_halos["FirstHaloInFOFgroup"]
                         != global_true_fof_idx)[0]
    global_flyby_ind = fof_halo_inds[flyby_ind]

    # Update the flybys and flip the `MostBoundID` to flag the flyby.
    forest_halos["FirstHaloInFOFgroup"][0:NHalos_root] = global_true_fof_idx
    forest_halos["MostBoundID"][global_flyby_ind] *= -1

    return forest_halos, global_true_fof_idx, global_flyby_ind


def fix_nextsubhalo(forest_halos, fof_groups, offset, NHalos, forestID,
                    snap_num):
    """
    Fixes the ``NextHaloInFOFgroup`` field for a single forest at a single
    snapshot.

    .. note::
        We pass all halos within the forest to this function but only alter
        those within a single snapshot.

    Parameters
    ----------

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The halos within a single forest.

    fof_groups : List of integers
        The ID of the host halo for each FoF group.

    offset : Integer
        The (global) offset for the halos within the snapshot we're updating.

    NHalos : Integer
        The number of halos being updated.

    forestID : Integer
        The forest ID for the forest we're updating.

    snap_num : Integer
        The snapshot we're altering.

    Returns
    ----------

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The forest halos with updated ``NextHaloInFOFgroup`` field.

    Errors 
    ----------

    RuntimeError
        We assume that all halos in a FoF group are stored contiguously in the
        input trees.  If this is not the case, a RuntimeError is raised.
    """

    # Every FoF group will point to a single halo, so loop over the FoF groups.
    for fof in fof_groups:

        # Find those halos within the snapshot we're altering in this FoF group.
        # We search only over the snapshot (offset -> offset+NHalos) to
        # increase efficiency for very large trees.
        halos_in_fof = np.where(forest_halos["FirstHaloInFOFgroup"][offset:offset+NHalos] == fof)[0]
        halos_in_fof_global_inds = halos_in_fof + offset

        # The first halo will point to index 1 and so on.
        nexthalo = np.arange(min(halos_in_fof_global_inds)+1,
                             min(halos_in_fof_global_inds)+len(halos_in_fof)+1)

        # We make the assumption that halos within a FoF group are stored
        # contiguously.  Ensure this is the case and that the indices of the
        # halos within the FoF group are simply an arange.
        if not np.allclose(halos_in_fof_global_inds, nexthalo-1):
            print("When attempting to fix the `NextHaloInFOFgroup` field for "
                  "ForestID {0} at snapshot {1} we encountered substructure that"
                  " was not stored contiguously.".format(forestID, snap_num))
            print("The offset value is {0} and NHalos is {1}".format(offset,
                                                                     NHalos))
            print("That is, the `hostHaloID` values were not stored "
                  "contiguously.") 
            print("For hostHaloID {0}, the halos within this FoF group were "
                  "{1}".format(fof, halos_in_fof_global_inds))
            print("The value of `nexthalo` is {0}. This should be equal to the"
                  " above list + 1".format(nexthalo))

            print(forest_halos[:][offset:offset+NHalos])
            raise RuntimeError

        # Check passed so we can update the halos.
        forest_halos["NextHaloInFOFgroup"][halos_in_fof_global_inds] = nexthalo

        # The final halo terminates with -1.
        forest_halos["NextHaloInFOFgroup"][halos_in_fof_global_inds[-1]] = -1

        #print("FoF group {0}\tNextHaloInFOFgroup {1}".format(fof, forest_halos["NextHaloInFOFgroup"][halos_in_fof_global_inds]))

    return forest_halos


def treefrog_to_lhalo(fname_in, fname_out, haloID_field="ID",
                      forestID_field="ForestID", Nforests=None,
                      write_binary_flag=1, fname_alist=None, dry_run=0,
                      total_files=size, files_processed=0, debug=0):
    """
    Takes the Treefrog trees that have had their IDs corrected to be in LHalo
    format and saves them in LHalo binary format.

    The data-structure of the Treefrog trees is assumed to be HDF5 File ->
    Snapshots -> Halo Properties at each snapshot.

    .. note::
        We require the input trees to be sorted via the forest ID
        (``forestID_field``) and suggest to also sub-sort on hostHaloID and mass.
        Sorting can be done using :py:mod:`astro3D.genesis.utils.forest_sorter`.

        We also require the input trees to have IDs that are LHalo compatible.
        See :py:mod:`astro3D.genesis.utils.convert_indices`.

    Parameters
    ----------

    fname_in, fname_out : String
        Path to the input HDF5 VELOCIraptor + treefrog trees and the path
        where the LHalo binary file will be saved.

    haloID_field : String, optional
        Field name within the HDF5 file that corresponds to the unique halo ID.

    forestID_field : String, optional
        Field name within the HDF5 file that corresponds to forest ID.

    Nforests : Integer, optional
        The number of forests to be processed. If ``None`` is passed then all
        forests are processed.

    write_binary_flag : Integer, optional
        Flag to decide whether to write to a binary or HDF5 file.
        0: HDF5 file only.
        1: Binary file only.
        2: Both binary and HDF5 file.

    fname_alist : String, optional
        If not ``None``, creates a file with path ``fname_alist`` that contains
        the scale factors at each snapshot (from lowest to highest). 

    total_files : Integer, optional
        The total number of files we splitting the trees across. Can be
        greater than then number of processors.

    files_processed : Integer, optional
         The number of files that have been processed by previous executions of
        the code.  Used to set the forest assignment offset and file number.

    debug : Integer
        Flag to denote whether extra debugging information should be printed to
        ``stdout``.

    Returns
    ----------

    None.

    Notes
    ----------

    The default parameters are chosen to match the ASTRO3D Genesis trees as
    produced by VELOCIraptor + Treefrog.
    """

    if debug:
        np.set_printoptions(threshold=np.nan)

    if rank == 0:
        print("")
        print("=================================")
        print("Going through the LHalo indices corrected Treefrog trees and "
              "saving in LHalo binary format.")
        print("Input Trees: {0}".format(fname_in))
        print("Output LHalo Trees: {0}".format(fname_out))
        print("ForestID Field Name: {0}".format(forestID_field))
        print("Number of Processors: {0}".format(size))
        print("Number of forests to process: {0}".format(Nforests))
        print("Write Binary Flag: {0}".format(write_binary_flag))
        print("=================================")
        print("")

    LHalo_Desc, multipledim_names = get_LHalo_datastruct()

    with h5py.File(fname_in, "r") as f_in:
        Snap_Keys, Snap_Nums = cmn.get_snapkeys_and_nums(f_in.keys())

        # Create a txt file that contains all the scalefactors.
        if fname_alist and rank == 0:
            write_alist(f_in, fname_alist, Snap_Keys)
            print("Saved alist to {0}".format(fname_alist))

        # If we're running in parallel, don't use the tqdm progress bar.
        if size > 1:
            is_mpi = True
        else:
            is_mpi = False

        NHalos_forest = cmn.get_halos_per_forest(f_in, Snap_Keys, haloID_field,
                                                 forestID_field, is_mpi, debug)

        # Find the max value of the object where the compared
        # values are returned via the "key". In this case,
        # compares the integer Snapshot number values, and then
        # returns the "Snapshot_group" key in the Snap_Keys dictionary.
        # Note to future Jacob: Manodeep wrote this comment so that's why it
        # doesn't make sense. Code Taken from
        # https://stackoverflow.com/questions/268272/getting-key-with-maximum-value-in-dictionary
        last_snap_key = max(Snap_Nums, key=Snap_Nums.get)

        total_forests_to_process = np.unique(f_in[last_snap_key][forestID_field][:])
        print("Total forests {0}".format(len(total_forests_to_process)))

        if Nforests:
            total_forests_to_process = total_forests_to_process[0:int(Nforests)]

        # If we're running in parallel, determine what forest IDs each
        # processor is handling.
        if size > 1:
            forests_to_process = determine_forests(NHalos_forest,
                                                   total_forests_to_process,
                                                   total_files=total_files,
                                                   files_processed=files_processed)
        else:
            forests_to_process = total_forests_to_process

        NHalos_rank = 0
        for forestID in forests_to_process:
            NHalos_rank += NHalos_forest[forestID] 

        print("I am rank {0} and I am processing {1} Forests with a total of "
              "{2} halos".format(rank, len(forests_to_process), NHalos_rank))

        # Now that each processor knows what Forests they're processing,
        # they need to know the number of halos per snapshot per forest.
        NHalos_forest_snap, NHalos_forest_snap_offset = \
            cmn.get_halos_per_forest_per_snap(f_in, Snap_Keys, haloID_field,
                                              forestID_field, is_mpi, debug,
                                              forests_to_process=forests_to_process)

        filenr = rank + files_processed

        # We first want to determine the number of forests, total number of
        # halos in these forests, and number of halos per forest for each
        # forest we are processing.
        totNHalos = 0
        global_halos_per_forest = []

        for forestID in forests_to_process:
            # NHalos_forest_snap is a nested dictionary accessed by each forestID.
            halos_per_forest = sum(NHalos_forest_snap[forestID].values())
            global_halos_per_forest.append(halos_per_forest)
            totNHalos += halos_per_forest

        # Write out the header with all this info.
        print("Rank {0} writing {1} forests containing a total of {2} "
              "halos.".format(rank, len(forests_to_process), totNHalos))
        if write_binary_flag == 1 or write_binary_flag == 2:
            my_fname_out = "{0}.{1}".format(fname_out, filenr)
        else:
            my_fname_out = "{0}.{1}.hdf5".format(fname_out, filenr)

        write_header(my_fname_out, len(forests_to_process), totNHalos,
                     global_halos_per_forest, write_binary_flag)

        start_time = time.time()
        hubble_h = get_hubble_h(f_in)
        # Now for each forest we want to populate the LHalos forest struct, fix
        # any IDs (e.g., flybys) and then write them out.
        if write_binary_flag == 1 or write_binary_flag == 2:
            f_out = open(my_fname_out, "ab")
        else:
            f_out = h5py.File(my_fname_out, "a")

        sum_NHalos_root = 0
        for count, forestID in enumerate(forests_to_process):
            #if count % int(len(forests_to_process)/10) == 0:
            #    print("Rank {0} processed {1} Forests ({2:.2f} seconds "
            #          "elapsed).".format(rank, count,
            #                             time.time()-start_time))
            #    print("Rank {0} processed {1} Root "
            #          "Halos".format(rank, sum_NHalos_root))

            NHalos = sum(NHalos_forest_snap[forestID].values())
            NHalos_root = NHalos_forest_snap[forestID][last_snap_key]
            sum_NHalos_root += NHalos_root

            forest_halos = np.zeros(NHalos, dtype=LHalo_Desc)
            forest_halos = populate_forest(f_in, forest_halos, Snap_Keys,
                                           Snap_Nums, forestID,
                                           NHalos_forest_snap,
                                           NHalos_forest_snap_offset,
                                           filenr, hubble_h)

            forest_NHalos_root = np.where(forest_halos["SnapNum"] ==
                                          Snap_Nums[last_snap_key])[0]

            if len(forest_NHalos_root) != NHalos_root:
                print("Using the NHalos_forest dictionary there should be {0} "
                      "halos at the root snapshot.".format(NHalos_root))
                print("However, after filling the forest halos from the HDF5 "
                      "file, there were only {0} halos at the root "
                      "snapshot.".format(forest_NHalos_root))
                print("Rank {0}, Forest {1}".format(rank, forestID))
                raise RuntimeError

            if dry_run:
                continue

            forest_halos, true_fof_idx, flyby_inds = fix_flybys(forest_halos,
                                                                NHalos_root)

            # Now if there were any flybys, we need to update the
            # `NextHaloInFOFgroup` chain to account for them.
            tmp = 0
            if true_fof_idx:
                # We do this by starting at the main FoF group and moving until we reach
                # the end (`NextHaloInFOFgroup = -1`).  Then we attach the first flyby halo
                # onto the end.  We then move down THAT flyby's chain and repeat the
                # process.
                next_in_chain = forest_halos["NextHaloInFOFgroup"][true_fof_idx]
                curr_halo = true_fof_idx

                for flyby_ind in flyby_inds:
                    while next_in_chain != -1:
                        curr_halo = next_in_chain
                        next_in_chain = forest_halos["NextHaloInFOFgroup"][next_in_chain]
                       
                        if debug:
                            print("Curr Halo: {0}\tNext In Chain: "
                                  "{1}\tflyby_inds: {2}\tFlyby Ind "
                                  "{3}\ttrue_fof_idx {4}" \
                                  .format(curr_halo, next_in_chain, flyby_inds,
                                          flyby_ind,
                                          true_fof_idx)) 

                            print(forest_halos["NextHaloInFOFgroup"])

                    forest_halos["NextHaloInFOFgroup"][curr_halo] = flyby_ind
                    next_in_chain = flyby_ind

                # After this there should only be one halo at the root
                # snapshot with `NextHaloInFOFgroup == -1`.
                assert(len(np.where(forest_halos["NextHaloInFOFgroup"][0:NHalos_root] \
                                    == -1)[0]) == 1)
            # Flybys and `NextHaloInFOFgroup` now fixed.

            forest_halos = fix_nextprog(forest_halos, forestID)

            # The VELOCIraptor + Treefrog trees point to themselves when
            # they terminate.  However LHalo Trees requires these to be -1,
            # so find the instances where  `FirstProgenitor` point to
            # themselves and adjust them to -1.
            w = np.arange(NHalos)
            #NextProg_tofix = [x for x in w if x == forest_halos["NextProgenitor"][x]]
            FirstProg_tofix = [x for x in w if x == forest_halos["FirstProgenitor"][x]]

            #forest_halos["NextProgenitor"][NextProg_tofix] = -1
            forest_halos["FirstProgenitor"][FirstProg_tofix] = -1

            # All done! Append to the file.

            if write_binary_flag == 1 or write_binary_flag == 2:
                forest_halos.tofile(f_out)
            else:
                group_name = "tree_{0:03d}".format(forestID)
                f_out.create_group(group_name)

                # Need to be careful here.  Some properties (e.g., 'Position') we
                # want to save as Nx3 arrays (instead of 3 Nx1 arrays).  So
                # first save only those arrays that aren't multi-dimensional.
                for subgroup_name in LHalo_Desc.names:
                    if not cmn.search_dict_of_lists(subgroup_name, multipledim_names):
                        f_out[group_name][subgroup_name] = forest_halos[subgroup_name]

                # Then go through all the multi-dimensional arrays and save
                # an Nx3 array that contains all the data.
                for name in multipledim_names:

                    # Initialize an Nx3 array.
                    Ndim = len(multipledim_names[name])
                    array = np.zeros((len(forest_halos), Ndim))

                    # Then populate that array.
                    for dim, dim_name in enumerate(multipledim_names[name]):
                        array[:, dim] = forest_halos[dim_name]

                    # Finally save it.
                    f_out[group_name][name] = array

        # End of Forests Loop.
        f_out.close()

    # Input HDF5 file closed.

    print("Rank {0} has finished writing out {1} Forests to "
          "{2}".format(rank, len(forests_to_process), my_fname_out))
    print("Rank {0} wrote a total of {1} Root Halos".format(rank, sum_NHalos_root))
    print("Total time elapsed: {0:.2f} Seconds.".format(time.time()-start_time))

    # If the user set `write_binary_flag == 2`, then convert the binary file to
    # a HDF5 one.

    if write_binary_flag == 2 and not dry_run:
        hdf5_fname_out = "{0}.{1}.hdf5".format(fname_out, rank)
        convert_binary_to_hdf5(my_fname_out, hdf5_fname_out)
        print("Binary file also converted to HDF5.")



def determine_forests(NHalos_forest, all_forests, total_files=size,
                      files_processed=0):
    """
    Load balances the number of halos across processors.

    We split the halos into ``total_files`` number of files. This can be
    greater than the number of processors currently running.  If it is, we
    still split the assignment to ``total_files`` but only process
    ``num_processors`` this execution. Only following executions, we use
    ``execution_number`` to properly offset the assignments.

    .. note::
        Since we do not split trees across processors, this function will
        result in processors having an unequal number of trees but similar
        total number of halos across processors.

    Parameters
    ----------

    NHalos_forest : Dictionary
        Dictionary that contains the number of halos for each Forest. Key is
        the ForestID.

    all_forests : List of integers
        List of Forest IDs to be processed across all files.

    total_files : Integer, optional
        The total number of files we splitting the trees across. Can be
        greater than then number of processors.

    files_processed : Integer, optional
         The number of files that have been processed by previous executions of
        the code.  Used to set the forest assignment offset.

    Returns
    ----------

    my_forest_assignment : List of integers
        Rank-unique list of forest IDs to be processed by this rank.
    """

    # We may be only using a subset of the total forests, so get these from the
    # dictionary.  If we're using all the forests, can do this quickly.
    # We want to get the forestIDs and the number of halos in the forest for
    # those forests we're interested in. 
    if len(all_forests) == len(NHalos_forest):
        NHalos_forest_interest = list(NHalos_forest.values())
        forestID_forest_interest = list(NHalosS_forest.keys())           
    else:
        NHalos_forest_interest = []
        forestID_forest_interest = all_forests
        for forestID in all_forests:
            NHalos_forest_interest.append(NHalos_forest[forestID])

    # Want to split the total number of halos across files.
    NHalos_total = sum(NHalos_forest_interest)
    NHalo_target = np.ceil(NHalos_total / total_files)

    # We want to partition the forests into buckets with number of halos
    # `NHalo_target`. To do this, we will generate a list with length
    # `total_files + 1`. This list will contain the index slice for the forests
    # processed for each file.
    assignment_idx = []
    assignment_idx.append(0)  # Want to start pulling from index 0.
    halos_counted = 0

    print("NHalos {0}\tTarget {1}".format(NHalos_total, NHalo_target))
    #print("NHalos_forest_interest {0}".format(NHalos_forest_interest))

    # Now go through all of the forests. When we have gathered `NHalo_target`
    # halos, remember the index.  Note: We do this manually instead of using a
    # cumulative sum because VERY large forests skew the assignment.
    for count, (forestID, this_NHalos) in enumerate(zip(forestID_forest_interest,
                                                        NHalos_forest_interest)):
        
        halos_counted += this_NHalos 
        if halos_counted >= NHalo_target:
            # Remember to +1 because slicing has an open right domain.
            # That is, if we want elements 0-8, the slice is [0:9].
            assignment_idx.append(count+1)
            halos_counted = 0

    # The last slice end will be the final element.
    assignment_idx.append(count+1)
 
    print("Assignment idx {0}".format(assignment_idx))
    # Assign the trees depending on the processor rank and how many files we've
    # processed so far.
    this_rank_idx = rank+files_processed
    forest_idx_low = assignment_idx[this_rank_idx] 
    forest_idx_high = assignment_idx[this_rank_idx+1]

    print("idx_low {0}\tidx_high {1}\tthis_rank_idx {2}".format(forest_idx_low, forest_idx_high, this_rank_idx))

    # We have the indices required to properly slice up the forests for each
    # processor. Do it! 
    forest_assignment = forestID_forest_interest[forest_idx_low:forest_idx_high]

    return forest_assignment


def write_header(fname_out, Nforests, totNHalos, halos_per_forest,
                 write_binary_flag):
    """
    Creates the LHalo Tree file and writes the header information.

    Parameters
    ----------

    fname_out : String
        Path to where the LHalo tree binary will be saved.

    Nforest : Integer
        Number of forests that will be saved in this file.

    totNHalos : Integer
        Total number of halos that will be saved in this file.

    halos_per_forest : List of integers
        The number of halos within each forest that will be saved in this file.

    write_binary_flag : Integer
        Flag to decide whether to write to a binary or HDF5 file.
        0: HDF5 file only.
        1: Binary file only.
        2: Both binary and HDF5 file.

    Returns
    ----------

    None.

    Notes
    ----------

    The header has the following data structure:

        Number of Forests within this file (Nforests): 4-byte integer.
        Total number of halos within this file: 4-byte integer.
        Number of halos within each forest for this file: Nforests*4-bytes integers
    """

    if write_binary_flag == 1 or write_binary_flag == 2:

        with open(fname_out, "wb") as f_out:
            f_out.write(np.array(Nforests, dtype=np.int32).tobytes())
            f_out.write(np.array(totNHalos, dtype=np.int32).tobytes())
            f_out.write(np.array(halos_per_forest, dtype=np.int32).tobytes())

    else:

        with h5py.File(fname_out, "w") as f_out:
            f_out.create_group("Header")
            f_out["Header"].attrs.create("Ntrees", Nforests, dtype=np.int32)
            f_out["Header"].attrs.create("totNHalos", totNHalos,
                                         dtype=np.int32)
            f_out["Header"].attrs.create("TreeNHalos", halos_per_forest,
                                         dtype=np.int32)

    return


def populate_forest(f_in, forest_halos, Snap_Keys, Snap_Nums, forestID,
                    NHalos_forest_snap, NHalos_forest_snap_offset, filenr,
                    hubble_h):
    """
    Takes an empty ``forest_halos`` structure and fills it with the halos
    corresponding to the passed ``forestID``.

    Parameters
    ----------

    f_in : Open HDF5 file
        The open HDF5 tree file we're reading from.

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        Data structure to hold all the halos in this forest.  Has been
        initalized to zeros with enough length to hold all halos within this
        forest.

    Snap_Keys : List of string
        List of keys that correspond to the fields containing the snapshot
        data.

    Snap_Nums : Dictionary of integers keyed by ``Snap_Keys``
        Snapshot number of each snapshot key.

    forestID : Integer
        The forest ID we're populating the halos for.

    NHalos_forest_snap : Nested Dictionary
        Nested dictionary that contains the number of halos for each Forest at
        each snapshot.  Outer-key is ``forestID`` and inner-key is the snapshot
        key.

    NHalos_forest_snap_offset : Nested Dictionary
        Nested dictionary that contains the offset for each Forest at each
        snapshot. Outer-key is ``forestID`` and inner-key is the snapshot key.

    filenr : Integer
        The output file number this tree will be written to.

    Returns
    ----------

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The filled data structure with the halos for the forest.
    """
    halos_offset = 0  # We need to slice the halos into the forest array in the
                      # proper place.

    # Start at the root redshift and work our way up the tree.
    for snap_key in Snap_Keys[::-1]:

        # Get the number, index offset and the corresponding indices for halos
        # at this snapshot.        
        try:
            NHalos = NHalos_forest_snap[forestID][snap_key]
        except KeyError:
            continue

        # Just another catch for empty snapshots.
        if NHalos == 0:
            continue

        halos_forest_offset = NHalos_forest_snap_offset[forestID][snap_key]
        halos_forest_inds = list(np.arange(halos_forest_offset,
                                           halos_forest_offset + NHalos))

        #print("Filling {0} Halos".format(len(halos_forest_inds)))
        # Then go to the HDF5 file and grab all the required properties.
        forest_halos, halos_offset = fill_LHalo_properties(f_in[snap_key],
                                                           forest_halos,
                                                           halos_forest_inds,
                                                           halos_offset,
                                                           Snap_Nums[snap_key],
                                                           filenr, hubble_h,
                                                           forestID)

    return forest_halos


def fill_LHalo_properties(f_in, forest_halos, halo_indices, current_offset,
                          snap_num, filenr, hubble_h, forestID):
    """
    Grabs all the properties from the input HDF5 file for the current snapshot.

    .. note::
        The passed parameter ``f_in`` should correspond to the snapshot group
        that contains the halo properties datasets.  E.g., if we are filling
        the halos for Snapshot 43 with key ``Snap_043``, this function should be
        called with ``f_in['Snap_043']``.

    Parameters
    ----------

    f_in: Open HDF5 file
        The snapshot level HDF5 group that we're reading from.

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        Data structure to hold all the halos for this forest.

    halo_indices : List of integers
        List of indices that correspond to the halos we're filling for this
        snapshot.

    current_offset : Integer
        The index into the ``forest_halo`` structure where this snapshot starts.

    snap_num : Integer
        The snapshot number we're filling for.

    filenr : Integer
        The output file number this tree will be written to.

    forestID : Integer
        The forest ID for the forest we're updating.

    Returns
    ----------

    forest_halos : Numpy structured array with data structure defined by
                  :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The filled data structure with the halos for the forest.
    """

    NHalos_thissnap = len(halo_indices)
    scale_factor = f_in.attrs['scalefactor']

    forest_halos["Descendant"][current_offset:current_offset+NHalos_thissnap] = f_in["Head"][halo_indices]
    forest_halos["FirstProgenitor"][current_offset:current_offset+NHalos_thissnap] = f_in["Tail"][halo_indices]

    # Initialize `NextProgenitor` with -1 and then fix it later with function
    # `fix_nextprog()`.
    forest_halos["NextProgenitor"][current_offset:current_offset+NHalos_thissnap] = -1

    # `FirstHaloInFOFgroup` is -1 for main FoF halos.  However in the LHaloTree
    # structure, this should point to itself. 
    forest_halos["FirstHaloInFOFgroup"][current_offset:current_offset+NHalos_thissnap] = f_in["hostHaloID"][halo_indices]
    #print(forest_halos["FirstHaloInFOFgroup"][current_offset:current_offset+NHalos_thissnap])

    # First find out what the FoF groups are.  Then go through the FoF groups
    # and update the `NextHaloInFOFgroup` pointer.
    forest_halos["NextHaloInFOFgroup"][current_offset:current_offset+NHalos_thissnap] = -1
    all_hosthalo_inds = f_in["hostHaloID"][halo_indices]
    fof_groups = np.unique(all_hosthalo_inds)
    forest_halos = fix_nextsubhalo(forest_halos, fof_groups, current_offset,
                                   NHalos_thissnap, forestID, snap_num)

    # All merger pointers are now done, read in the halo properties.
    forest_halos["Len"][current_offset:current_offset+NHalos_thissnap] = f_in["npart"][halo_indices]
    forest_halos["M_Mean200"][current_offset:current_offset+NHalos_thissnap] = f_in["Mass_200mean"][halo_indices] * hubble_h
    forest_halos["Mvir"][current_offset:current_offset+NHalos_thissnap] = f_in["Mass_200crit"][halo_indices] * hubble_h
    forest_halos["M_TopHat"][current_offset:current_offset+NHalos_thissnap] = f_in["Mass_200crit"][halo_indices] * hubble_h

    # Positions are in Co-moving Units #
    forest_halos["Posx"][current_offset:current_offset+NHalos_thissnap] = f_in["Xc"][halo_indices] * hubble_h / scale_factor
    forest_halos["Posy"][current_offset:current_offset+NHalos_thissnap] = f_in["Yc"][halo_indices] * hubble_h / scale_factor
    forest_halos["Posz"][current_offset:current_offset+NHalos_thissnap] = f_in["Zc"][halo_indices] * hubble_h / scale_factor

    # Velocities are in Physical Units #
    forest_halos["Velx"][current_offset:current_offset+NHalos_thissnap] = f_in["VXc"][halo_indices]
    forest_halos["Vely"][current_offset:current_offset+NHalos_thissnap] = f_in["VYc"][halo_indices]
    forest_halos["Velz"][current_offset:current_offset+NHalos_thissnap] = f_in["VZc"][halo_indices]

    forest_halos["VelDisp"][current_offset:current_offset+NHalos_thissnap] = f_in["sigV"][halo_indices]
    forest_halos["Vmax"][current_offset:current_offset+NHalos_thissnap] = f_in["Vmax"][halo_indices]

    # The 'spin' parameter in LHalo Tree is the Angular Momentum divided by the
    # total mass.
    M_tot = f_in["Mass_tot"][halo_indices]

    forest_halos["Spinx"][current_offset:current_offset+NHalos_thissnap] = \
    f_in["Lx"][halo_indices] * hubble_h * hubble_h / M_tot

    forest_halos["Spiny"][current_offset:current_offset+NHalos_thissnap] = \
    f_in["Ly"][halo_indices] * hubble_h * hubble_h / M_tot

    forest_halos["Spinz"][current_offset:current_offset+NHalos_thissnap] = \
    f_in["Lz"][halo_indices] * hubble_h * hubble_h / M_tot

    forest_halos["MostBoundID"][current_offset:current_offset+NHalos_thissnap] = f_in["oldIDs"][halo_indices]
    forest_halos["SnapNum"][current_offset:current_offset+NHalos_thissnap] = snap_num
    forest_halos["Filenr"][current_offset:current_offset+NHalos_thissnap] = filenr

    forest_halos["SubHaloIndex"][current_offset:current_offset+NHalos_thissnap] = -1
    forest_halos["SubHalfMass"][current_offset:current_offset+NHalos_thissnap] = -1
    current_offset += NHalos_thissnap

    return forest_halos, current_offset


def convert_binary_to_hdf5(fname_in, fname_out):
    """
    Converts a binary LHalo Tree file to HDF5 format.

    Parameters
    ----------

    fname_in, fname_out : String
        Path to the input LHalo binary tree and the path
        where the HDF5 tree will be saved.

    Returns
    ----------

    None.

    Notes
    ----------

    Multi-dimensional arrays (e.g., 'Position') are saved as Nx3 arrays.
    """

    LHalo_Desc, multipledim_names = get_LHalo_datastruct()

    with open(fname_in, "rb") as binary_file, \
         h5py.File(fname_out, "w") as hdf5_file:

        # First get header info from the binary file.
        NTrees = np.fromfile(binary_file, np.dtype(np.int32), 1)[0]
        NHalos = np.fromfile(binary_file, np.dtype(np.int32), 1)[0]
        NHalosPerTree = np.fromfile(binary_file,
                                    np.dtype((np.int32, NTrees)), 1)[0]

        print("For file {0} there are {1} trees with {2} total halos"
              .format(fname_in, NTrees, NHalos))

        # Write the header information to the HDF5 file.
        hdf5_file.create_group("Header")
        hdf5_file["Header"].attrs.create("Ntrees", NTrees, dtype=np.int32)
        hdf5_file["Header"].attrs.create("totNHalos", NHalos, dtype=np.int32)
        hdf5_file["Header"].attrs.create("TreeNHalos", NHalosPerTree,
                                         dtype=np.int32)

        # Now loop over each tree and write the information to the HDF5 file.
        for tree_idx in trange(NTrees):
            binary_tree = np.fromfile(binary_file, LHalo_Desc,
                                      NHalosPerTree[tree_idx])

            tree_name = "tree_{0:03d}".format(tree_idx)
            hdf5_file.create_group(tree_name)

            # Need to be careful here.  Some properties (e.g., position) we
            # want to save as Nx3 arrays (instead of 3 Nx1 arrays).  So
            # first save only those arrays that aren't multi-dimensional.
            for subgroup_name in LHalo_Desc.names:
                if not cmn.search_dict_of_lists(subgroup_name, multipledim_names):
                    hdf5_file[tree_name][subgroup_name] = binary_tree[subgroup_name]

            # Then go through all the multi-dimensional arrays and save
            # an Nx3 array that contains all the data.
            for name in multipledim_names:

                # Initialize an Nx3 array.
                Ndim = len(multipledim_names[name])
                array = np.zeros((len(binary_tree), Ndim))

                # Then populate that array.
                for dim, dim_name in enumerate(multipledim_names[name]):
                    array[:, dim] = binary_tree[dim_name]

                # Finally save it.
                hdf5_file[tree_name][name] = array


def get_hubble_h(f_in):
    """
    Gets the value of Hubble little h.

    Parameters
    ----------

    f_in : Open HDF5 File
        The open HDF5 file we're reading the data from.

    Returns
    ----------

    hubble_h : Float
        The value of Hubble little h for the given cosmology.
    """

    #hubble_h = f_in["Header"]["Cosmology"].attrs["h_val"]
    hubble_h = 0.6751

    return hubble_h


def write_alist(f_in, fname_alist, Snap_Keys):
    """
    Saves the scale factor values of the snapshots to a file.

    Parameters
    ----------

    f_in : Open HDF5 File
        The open HDF5 file we're reading the data from.

    fname_alist : String
        Full path to the file we're writing to..

    Snap_Keys : List of strings
        Names of the snapshot keys used to access the HDF5 file. 

    Returns
    ----------

    None.
    """

    alist = []

    for key in Snap_Keys:
        snap_a = f_in[key].attrs["scalefactor"]

        alist.append(snap_a)

    np.savetxt(fname_alist, alist)    
