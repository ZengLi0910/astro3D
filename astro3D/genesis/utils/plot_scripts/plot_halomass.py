"""
This file plots the Halo Mass Function of LHalo Binary/HDF5 trees. 
"""
#!/usr/bin:env python
from __future__ import print_function

from astro3D.genesis.utils import treefrog_to_lhalo as frog_to_l

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py

colors = ["k", "#dd1c77", "#3182bd", "#f03b20", "#31a354"] 
dashes = ['',
          [3, 3, 3, 3],
          [7, 1, 1, 1],
          [1, 1, 1, 1],
          [5, 1, 5, 1]]
fontsize = 22
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize=fontsize)
plt.rc("ytick", labelsize=fontsize)
plt.rc("axes", labelsize=fontsize)


def load_LHaloTree_halos(fname):
    """
    Loads the Halos of an LHaloTree binary file.

    Parameters
    ----------

    fname : String.
        The full path to the file.

    Returns
    ----------

    halos : Numpy structured array with data structure defined by
            :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The halos from the file.
    """

    LHalo_Desc, _ = frog_to_l.get_LHalo_datastruct()

    with open(fname, "rb") as f_in:

        Nforests, totNHalos, halos_per_forest = read_LHaloTree_header(f_in)

        halos = np.empty(totNHalos, dtype=LHalo_Desc)
        halos = np.fromfile(f_in, LHalo_Desc, totNHalos) 

    return halos


def load_genesis_halo_mass(fname, snapnum):
    """
    Loads the Halos of the Genesis HDF5 file at a specific snapshot. 

    Parameters
    ----------

    fname : String.
        The full path to the file.

    snapnum : Integer.
        The snapshot number we're reading.

    Returns
    ----------

    halo_mass : List of integers. 
        The mass of the halos at the specified snapshot. Units are log10(Msun).
    """

    snap_key = "Snap_{0:03d}".format(snapnum)

    with h5py.File(fname, "r") as f_in:

        hubble_h = f_in["Header"]["Cosmology"].attrs["h_val"]

        mass = np.log10(f_in[snap_key]["Mass_200mean"][:] * 1.0e10)# / hubble_h)

    return mass


def read_LHaloTree_header(f_in):
    """
    Reads the header information of an *open* LHaloTree binary file.

    Parameters
    ----------

    f_in : Open file. 
        The file we're reading from.

    Returns
    ----------

    Nforests : Integer.
        Number of forests in the file.

    totNHalos : Integer.
        Total number of halos in the file.

    halos_per_forest : List of integers.
        The number of halos within each forest in the file.

    Notes
    ----------

    The header has the following data structure:

        Number of Forests within this file (Nforests): 4-byte integer.
        Total number of halos within this file: 4-byte integer.
        Number of halos within each forest for this file: Nforests*4-bytes integers
    """

    Nforests = np.fromfile(f_in, np.dtype(np.int32), 1)[0]
    totNHalos = np.fromfile(f_in, np.dtype(np.int32), 1)[0]
    halos_per_forest = np.fromfile(f_in, np.dtype(np.int32), int(Nforests))[0]

    return Nforests, totNHalos, halos_per_forest
 

def update_HMF(current_HMF, halos, mstar_bins, hubble_h, mass=None, snapnum=None):
    """
    Updates the current Halo Mass Function.

    Parameters
    ----------

    current_HMF : List of integers.
        The count of halos in each mass bin (defined by `mstar_bins`).

    halos : Numpy structured array with data structure defined by
            :py:mod:`astro3D.genesis.utils.treefrog_to_lhalo.get_LHalo_datastruct`
        The halos used to update the function.

    mstar_bins : List of integers.
        The mass bins we bin the halos into.  Units are log10(Msun).

    hubble_h : Float.
        The little h Hubble parameter associated with these halos.

    snapnum : Integer, optional.
        The snapshot number we're updating the HMF for. Only halos at this
        snapshot will be counted.

    Returns
    ----------

    updated_HMF : List of integers.
        The count of halos in each mass bin after adding the current `halos`.

    Notes
    ----------

    The input halos are assumed to have mass in units of 1.0e10/`hubble_h`
    Msun, consistent with the `LHaloTree` specifications.  

    """ 

    if snapnum:
        w = np.where(halos["SnapNum"] == snapnum)[0]
    else:
        w = np.arange(len(halos))

    if mass is None:
        mass = np.log10(halos["Mvir"][w] * 1.0e10 / hubble_h)

    this_HMF = np.histogram(mass, bins=mstar_bins)

    updated_HMF = current_HMF + this_HMF[0]

    return updated_HMF

def plot_HMF(mstar_bins, mstar_bin_width, HMF, model_tags, output_dir,
             output_tag, output_format):
    """
    Plots the halo mass function. That is, the number count of halos 
    binned on mass. 

    Parameters
    ----------

    mstar_bins : List of floats
        Mass bins that the data is binned on. Units are log10(Msun).

    mstar_bin_width : Float
        The bin separation between the mass bins. Units are log10(Msun). 

    HMF : 2D nested lists of floats. Outer length is number of models, next is
          the number of ``mstar_bins``.
        The halo mass function at each snapshot for each model.  That is,
        the number of galaxies within each mass bin (given by `mstar_bins`).
        Normalized by the volume of the simulation, units are `Count*Mpc^-3`. 

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig1, ax = plt.subplots(nrows=1, ncols=1,
                            sharex='col', sharey='row', figsize=(16,6))

    for model_number in range(len(HMF)):

        label = model_tags[model_number]
        ax.plot(mstar_bins[:-1] + mstar_bin_width*0.5,
                HMF[model_number][0],
                color=colors[model_number], 
                dashes=dashes[model_number],
                label=label)
            
        ax.set_xlim([min(mstar_bins)-0.3, max(mstar_bins)+0.3])
        ax.set_xlabel(r'$\mathbf{log_{10} \: M \:[M_{\odot}]}$', 
                      fontsize=fontsize)

    # Since y-axis is shared, only need to do this once.
    ax.set_ylabel(r'$\mathbf{log_{10} \: \Phi\ [Mpc^{-3}\: dex^{-1}]}$', 
                  fontsize=fontsize) 
    ax.set_yscale('log', nonposy='clip')

    leg = ax.legend(loc='lower left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(12)

    plt.tight_layout()

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


if __name__ == '__main__':

    fname_base = ["/fred/oz070/jseiler/astro3d/nov2018/N1024_converted",
                  "/fred/oz009/N1024/unifiedcatalogs/VELOCIraptor.tree.t4.unifiedhalotree.links.snap.hdf.data",
                  "/fred/oz070/jseiler/astro3d/nov2018/N1024_hosthaloIDfixed.hdf5",
                  "/fred/oz070/jseiler/astro3d/nov2018/N1024_sorted.hdf5",
                  "/fred/oz070/jseiler/astro3d/nov2018/N1024_lhalo_indices.hdf5"]
    hubble_h = [0.671, 0.671, 0.671, 0.671, 0.671]
    boxsize = [500, 500, 500, 500, 500] # Mpc/h.
    snapnum = [100, 100, 100, 100, 100]
    file_type = ["binary", "hdf5", "hdf5", "hdf5", "hdf5"]

    model_tags = ["Mine N1024", "Base", "hostHaloIDfixed", "sorted",
                  "lhalo-indices"]
    output_dir = "."
    output_format = "png"

    num_models = 1
    num_files = 32

    # Parameters for the binning.
    mstar_bin_low = 9.0
    mstar_bin_high = 14.0
    mstar_bin_width = 0.2
    mstar_Nbins = int((mstar_bin_high - mstar_bin_low) / mstar_bin_width)
    mstar_bins = np.arange(mstar_bin_low,
                           mstar_bin_high + mstar_bin_width,
                           mstar_bin_width)

    # Halo Mass Function.
    HMF_allmodels = []

    for model_number in range(num_models):

        HMF_allmodels.append([])
        HMF_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                    dtype=np.float32))

        if file_type[model_number] == "binary": 

            NHalos = 0
            for filenr in range(num_files):

                fname = "{0}.{1}".format(fname_base[model_number], filenr)
                halos = load_LHaloTree_halos(fname)
                NHalos += len(halos)

                HMF_allmodels[model_number] = update_HMF(HMF_allmodels[model_number],
                                                         halos, mstar_bins,
                                                         hubble_h[model_number],
                                                         snapnum=snapnum[model_number])

            print("Mine: ALL NHalos {0}".format(NHalos))
            print("Mine: Snap {0} NHalos {1}".format(snapnum, np.sum(HMF_allmodels[model_number]))) 

        elif file_type[model_number] == "hdf5":

            fname = fname_base[model_number]

            NHalos = 0
            for snap in np.arange(0, 101):
                mass = load_genesis_halo_mass(fname, snap) 
                NHalos += len(mass)
                print("Snap {0} NHalos {1}".format(snap, len(mass))) 
            print("TOTAL Halos {0}".format(NHalos))
            print("")
            HMF_allmodels[model_number] = update_HMF(HMF_allmodels[model_number],
                                                     [], mstar_bins,
                                                     hubble_h[model_number],
                                                     mass=mass)


        else:
            print("The only valid options for `file_type` is 'binary' or "
                  "'hdf5'. For Model {0} you entered {1}".format(model_number,
                  file_type[model_number]))
        print("Model {0} HMF {1} NHalos {2}".format(model_number, HMF_allmodels[model_number],
                                                    np.sum(HMF_allmodels[model_number])))
        HMF_allmodels[model_number] = np.divide(HMF_allmodels[model_number],
                                                pow(boxsize[model_number], 3) /
                                                pow(hubble_h[model_number], 3) *
                                                mstar_bin_width)



    plot_HMF(mstar_bins, mstar_bin_width, HMF_allmodels, model_tags,
             output_dir, "HMF", output_format)
