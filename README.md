#genesis/utils

This directory contains a number of useful tools and utilities for handling the data from the
ASTRO3D Genesis simulations.

#Forest Sorter


This utility takes input HDF5 merger trees that have not been saved in any specific order and sorts
them on two fields; usually an ID field and a mass field.  The user can specify on which two fields 
the sorting should occur. In the default case, the trees are sorted first on the ForestID each halo 
belongs to and then the 200mean mass of the halo.  This results in halos that are sorted in
ascending order according to their ForestID, then within each Forest, sorted in ascending order 
according to their 200mean mass.

##Tests

First please run the basic tests on the default test data provided by invoking ``pytest``.  If this
default test does not pass, please email jseiler@swin.edu.au 

You are also able to run the test on your own provided data.  To do this run ``python3
tests/forest_sorter_test.py`` with command line arguments ``-f`` to specify your data file name and
``-n`` many halos you wish to test on (default is 10,000).  **If you are not using the Treefrog
merger tree data format you may need to specify more command line arguments.  Refer to the section
other merger tree formats**

``$ python3 tests/forest_sorter_test.py --fname_in=/Path/To/my_data.hdf5 --NHalos_test=1000``

If the default test passes but your specific test fails please ensure that your data file is not
corrupt.  Importantly, check that the snapshot keys are named appropriately.  We require the
snapshot fields to include the word **snap** (case insensitive) and assume that the snapshot number
corresponding to the snapshot key is included as a single cluster towards the end of the key;
**snap53_04** should correspond to snapshot number 04 for example. 

**If the snapshot fields are named correctly and your data can be otherwise read in via h5py, please
email jseiler@swin.edu.au**

###Testing on sorted products

Once you have run `forest_sorter` on your HDF5 trees you may wish to ensure
that the sorted data file is correct. By default, `tests/forest_sorter_test.py`
generates a small subset of halos from the provided HDF5 trees to test on.
This behaviour can be turned off by specifying `--gen_data=0`.  In this case,
you must specify `--fname_in` to be the **original** HDF5 trees and
`--fname_out` to be the **sorted** HDF5 data file. 

##Usage

``forest_sorter`` is run using command line arguments.  A full list is included below; to show this
help message in terminal run ``python3 forest_sorter.py -h``.  The two required arguments is the 
path for the input HDF5 data file, ``--fname_in=/Path/To/my_data.hdf5`` and the path for the output 
sorted HDF5 data file, ``--fname_out=/Path/To/output_sorted_data.hdf5``.

``$ python3 forest_sorter.py --fname_in=/Path/To/my_data.hdf5
--fname_out=/Path/To/output_sorted_data.hdf5``

#Command Line Options and Other Merger Tree Formats

The default settings of the command line options were chosen to run on the Treefrog
merger tree format.  If you are using a different merger tree format, or the field names have been
changed, you will need to specify more command line arguments at execution.  

If you are running with different a merger tree format, we strongly recommend that you run the
testing script with your data.  The field name command line arguments required for running the 
testing case and ``forest_sorter.py`` are identical. 

When specifying multiple arugments for an option (e.g., for ``sort_fields``),
use a comma to separate the arguments **with no white space**.  E.g.,

``$ python3 forest_sorter.py -s ForestID,Mass_200mean,NextSubhalo``

or

``$ python3 forest_sorter.py --sort_fields=ForestID,Mass_200mean,NextSubhalo`` 

usage: forest_sorter.py [-h] [-f FNAME_IN] [-o FNAME_OUT] [-s SORT_FIELDS]
                        [-i HALO_ID] [-p ID_FIELDS] [-x INDEX_MULT_FACTOR]

optional arguments:
  -h, --help            show this help message and exit
  -f FNAME_IN, --fname_in FNAME_IN
                        Path to the input HDF5 data file. Required.
  -o FNAME_OUT, --fname_out FNAME_OUT
                        Path to the output HDF5 data file. Required.
  -s SORT_FIELDS, --sort_fields SORT_FIELDS
                        Field names we will be sorted on. ORDER IS IMPORTANT.
                        Order using the outer-most sort to the inner-most.
                        Separate each field name with a comma. Default:
                        ForestID,Mass_200mean.
  -i HALO_ID, --HaloID HALO_ID
                        Field name for halo ID. Default: ID.
  -p ID_FIELDS, --ID_fields ID_FIELDS
                        Field names for those that contain IDs. Separate field
                        names with a comma. Default:
                        ID,Tail,Head,NextSubHalo,Dummy1,Dumm2).
  -x INDEX_MULT_FACTOR, --index_mult_factor INDEX_MULT_FACTOR
                        Conversion factor to go from a unique, per-snapshot
                        halo index to a temporally unique haloID. Default:
                        1e12.

usage: forest_sorter_test.py [-h] [-f FNAME_IN] [-o FNAME_OUT]
                             [-s SORT_FIELDS] [-i HALO_ID] [-p ID_FIELDS]
                             [-x INDEX_MULT_FACTOR] [-n NHALOS_TEST]
                             [-g GEN_DATA]

optional arguments:
  -h, --help            show this help message and exit
  -f FNAME_IN, --fname_in FNAME_IN
                        Path to test HDF5 data. Default: /test_data.hdf5
  -o FNAME_OUT, --fname_out FNAME_OUT
                        Path to sorted output HDF5 data file. Default:
                        /test_sorted.hdf5
  -s SORT_FIELDS, --sort_fields SORT_FIELDS
                        Field names we will be sorted on. ORDER IS IMPORTANT.
                        Order using the outer-most sort to the inner-most.
                        Separate each field name with a comma. Default:
                        ForestID,Mass_200mean.
  -i HALO_ID, --HaloID HALO_ID
                        Field name for halo ID. Default: ID.
  -p ID_FIELDS, --ID_fields ID_FIELDS
                        Field names for those that contain IDs. Separate field
                        names with a comma. Default:
                        ID,Tail,Head,NextSubHalo,Dummy1,Dumm2).
  -x INDEX_MULT_FACTOR, --index_mult_factor INDEX_MULT_FACTOR
                        Conversion factor to go from a unique, per-snapshot
                        halo index to a temporally unique haloID. Default:
                        1e12.
  -n NHALOS_TEST, --NHalos_test NHALOS_TEST
                        Minimum number of halos to test. Default: 10,000
  -g GEN_DATA, --gen_data GEN_DATA
                        Flag whether we want to generate data. If this is set
                        to 0, the tests will be run on the `fname_out` sorted
                        data that was created running on `fname_in`. Default:
                        1.

