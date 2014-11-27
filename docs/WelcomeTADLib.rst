.. TADLib documentation master file, created by
   sphinx-quickstart on Thu Oct 30 10:18:28 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TADLib: A Library to Explore Chromatin Interaction Patterns for Topologically Associating Domains
*************************************************************************************************

Introduction
============
Chromosome conformation capture (3C) derived techniques, especially Hi-C,
have revealed that topologically associating domain (TAD) is a structural
basis for both chromatin organization and regulation in three-dimensional
(3D) space. To systematically investigate the relationship between structure
and function, it is important to develop a quantitative parameter (or feature)
to measure the structural characteristics of TAD. TADLib is a package to explore
the chromatin interaction pattern of TAD.

Inspired by the observation that there exist great differences in chromatin
interaction pattern and gene expression level among TADs, a chromatin interaction
feature is developed to capture the aggregation degree of long-range chromatin
interactions. Application to human and mouse cell lines shows that there
exist heterogeneous structures among TADs and the structural rearrangement across
cell lines is significantly associated with transcription activity remodeling.

TADLib package is written in Python and provides a four-step pipeline:

- Identifying TAD from Hi-C data (optional)
- Selecting long-range chromatin interactions in each TAD
- Finding the aggregation patterns of selected interactions
- Calculating chromatin interaction feature of TAD

.. note:: By default, we suppose that the input Hi-C data are corrected appropriately.
   Otherwise, systematic biases in source data will negatively impact chromatin
   interaction selection and then parameter calculation. Several correction schemes
   are available online [1,2,3]_

Links
=====
- `PyPI <https://pypi.python.org/pypi/TADLib>`_ (Download and Installation)
- `Repository <https://github.com/XiaoTaoWang/TADLib>`_ (At GitHub)

Modules
=======
coniss
    Provide an interface to identify TAD from Hi-C data by using *CONISS* model,
    a constrained hierarchical clustering algorithm. However, TADs identified from
    other methods can also be used if they are preferred.

analyze
    Core module. The operations include chromatin interaction selection, interaction
    block generation and feature calculation.

polygon
    Complementary module. This module includes convex hull generation and other
    convex-based operations.

Requirements
============
TADLib is developed and tested on UNIX-like operating system, and following Python
packages are recommended:

- Python (2.x >= 2.6, not compatible with 3.x)
- Numpy (>= 1.6)
- Scipy library (>= 0.10)
- Matplotlib (>= 1.1)
- scikit-learn (>= 0.10)

For *coniss*, additional Python package *rpy2* (a powerful python package
providing simple and robust access to R within python) and two R packages
*rioja* and *vegan* are required.

.. note:: Tested systems: Red Hat Enterprise Linux Server release 6.4 (Santiago),
   CentOS release 6.4 (Final), Fedora release 20 (Heisenbug)

Installation
=============
Firstly, install all the required python packages:

There's an exhaustive instruction at http://www.scipy.org/install.html

We strongly recommend using `conda <http://conda.pydata.org/miniconda.html>`_,
an excellent Python package and environment manager.

Once Miniconda is installed, you can use the conda command to install any
other packages. Open a terminal and type::

    $ conda install numpy scipy matplotlib scikit-learn

More details about conda can be found at http://conda.pydata.org/docs/

Then you can install TADLib just as other packages stored in PyPI:

Use *easy_install*::

    $ conda install setuptools
    $ easy_install install TADLib

Or download the `source code <https://pypi.python.org/pypi/TADLib>`_ manually,
extract it and run the setup.py script::

    $ python setup.py install

Finally, run this command in a terminal prompt::

    $ python -m "tadlib.tests.testall"

If no exception occurs, congratulations, TADLib has been installed successfully!

Let's spend a few minutes dealing with *coniss* dependencies:

Download R source code here: http://cran.rstudio.com/, unpack it, change to
the extracted directory, get to a terminal prompt and type::

    $ ./configure --prefix=R_HOME --enable-R-shlib
    $ make
    $ make check
    $ make install

where R_HOME is the installation directory you choose for R.

You may need to reset your PATH and LD_LIBRARY_PATH environment variables
then. Use vi command to add this two lines to ``~``/.bashrc:

- export PATH=R_HOME/bin:$PATH
- export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:R_HOME/lib64/R/lib

To update the environment variables::

    $ source ~/.bashrc

Then it's time to install *rioja* and *vegan*. Open a R prompt::

    > install.packages("rioja")
    > install.packages("vegan")

rpy2 can be installed in this way::

    $ conda install pip python=2.7.5=2
    $ pip install rpy2 singledispatch

To test if rpy2 was installed correctly, run::

    $ python -m "rpy2.tests"


QuickStart
==========
A sample data is also distributed with the source code. First, change to the
sample directory: (we suppose you are still in the distribution root directory)::

    $ cd lib/data

Open a Python interperter and load Hi-C data:

>>> from tadlib import analyze
>>> # Data happen to lie in the working directory
>>> workInters = analyze.Inters(path = '.', template = 'h1hesc.chr%s_chr%s.txt')
>>> workInters.data['1'].shape # Chromosome 1
(10535,)

Load TAD regions:

>>> source = 'h1hesc.domain.txt'
>>> workTAD = analyze.TAD(source)
>>> workTAD.data[0]
('1', 5570000, 6090000)

Consider such a TAD region: chr1:5570000-6090000

Initialize a *Core* object using interaction matrix:

>>> left = 5570000 / 10000 # Resolution 10000
>>> right = 6090000 / 10000
>>> matrix = analyze.getmatrix(workInters.data['1'], left, right)
>>> workCore = analyze.Core(matrix)

Extract long-range interactions:

>>> workCore.longrange()
>>> len(workCore.pos) # Number
45

Aggregation pattern discovery:

>>> workCore.DBSCAN()
>>> workCore.Nc # Cluster Number
2

Feature calculation:

>>> workCore.gdensity()
>>> print workCore.gden # Subjects to small variations across platforms
0.463978520686

If you have installed *coniss* successfully, you can identify TADs as follows:

>>> from tadlib import coniss
>>> workQueue = coniss.Queue(template = 'h1hesc.chr%s_chr%s.txt')
>>> workQueue.TAD.dtype.names
('chr', 'start', 'end')

CALFEA
======
We have wrapped this pipeline into a user-friendly software called CALFEA
(CALculate FEAture for TADs).

Usage
-----
``calfea <-f filename> <--callTAD | --loadTAD> [options]``

Example and Explanation
-----------------------
Use the sample data again. Create a new directory under the distribution root::

    $ mkdir working

Change current directory to "working"::

    $ cd working

And type in command below::

    $ calfea --callTAD -f test1_feature.txt -p ../lib/data -F TXT -R 10000 -T h1hesc.chr%s_chr%s.txt -c 0 1 2 --immortal -S test1 -O test1 --verbose 3

As an example, we present most available parameters here.

- ``--callTAD``

At first, you should specify a flag ``--callTAD/--loadTAD`` to tell calfea whether
to identify TADs itself or just load from an existing file.

- ``-f/--filename`` FILENAME

Output file name. The output lines have 5 fields: *ChromID*, *TAD Start*,
*TAD End*, *Density* and *Gap Ratio*. We trace *Gap Ratio* for each TAD because
gap regions are always eliminated from original interaction matrix. By default,
results are printed to the standard output streams ("stdout").

- ``-p/--path`` PATH

The folder containing Hi-C data. Both absoulte and relative path are okay.
Chromatin are partitioned into bins under a certain resolution and read pairs
are assigned to these bins. You have to store Hi-C data chromosome by chromosome
following common naming pattern. For example, under pattern "chr%s_chr%s.int",
file names shall look like "chr1_chr1.int", "chr2_chr2.int", and so forth. Note
only intra-chromosome data are allowed, and don't place inter-chromosome ones
under the folder. (Default: .)

- ``-F/--format`` {TXT,NPZ}

Hi-C data format. Data should be provided in TXT format when you run *calfea*
for the first time.  If ``--immortal`` is specified, loaded data will be
saved to a new file in *.npz* format, which will speed up data loading process
greatly next time. (Default: TXT)

- ``-R/--resolution`` RESOLUTION

Resolution of the binned data. (Default: 10000)

- ``-T/--template`` TEMPLATE

Naming pattern for matching Hi-C data file names. Regular expression is used.
(Default: chr%s_chr%s.int)

- ``-c/--cols`` COLS COLS COLS

3 columns reading from TXT source data, with 0 being the first. Here,
``-c 0 1 2`` will extract the 1st, 2nd and 3rd columns.

- ``--immortal``

When specified, a Numpy .npz file will be generated under **PATH**.

- ``-S/saveto`` SAVETO

Of course, you need to give the prefix of output *.npz* file name with
``--immortal``.

- ``-O/--output`` OUTPUT

Prefix of the TAD file which will be generated under **PATH**. Required by
``--callTAD``.

- ``--verbose`` {0,1,2,3}

Logging level. 0: only show error messages, 1: also report warnings, 2:
show process information, 3: debug. We log all messages of all levels to a
disk file (calfea.log), while simultaneously logging process information or
above to the console. (Default: 2)

After this command, two files **test1_feature.txt** and **calfea.log** are
created under current working directory, and another two **test1.domain.txt**
and **test1.npz** are created under the sample data folder "../lib/data".

TAD identification is time/memory consuming, so you'd better load an existing
TAD file for the second time::

    $ calfea --loadTAD -s ../lib/data/test1.domain.txt -p ../lib/data -F NPZ -P test1 -f test2_feature.txt --verbose 3

Note that we use NPZ format this time and ``-T``, ``-c``, ``--immortal``,
``-S`` and ``-O`` are no longer needed. Flag ``-s`` and ``-P`` may be strange
for you:

- ``-s/--source`` SOURCE

Complete TAD source file name. The source file must contain 3 columns,
indicating *ChromID*, *TAD Start* and *TAD End*, respectively.

- ``-P/--prefix`` PREFIX

Prefix of input *.npz* file name. Path should not be included. 

This time, only **test2_feature.txt** is generated. We use a rotating file
for logging. According to our settings, when the size of **calfea.log**
gets about 15K, it's closed and renamed to **calfea.log.1**, and at the same
time, a new file **calfea.log** is silently opened for output. In a word,
the system saves old log files by appending the extensions ".1", ".2" etc.,
and the file being written to is always **calfea.log**.

Other Options
-------------

- ``--offset`` OFFSET

Offset from the main diagonal of TAD interaction matrix. All entries with
offset less than this value will be ignored during our analysis. (Default: 2)

- ``-C/--chroms`` [CHROMS [CHROMS ...]]

Chromosomes to read. Specially, "#" stands for chromosomes with
numerical labels. ``--chroms`` with zero argument will generate an empty
list, in which case all chromosome data will be loaded. (Default: ['#', 'X'])

- ``--chromname`` CHROMNAME

Template of chromosome names contained in external TAD file. If you specify
``--chromosome chr``, the content should like this::

    chr1    0   350000
    chr1    350000  800000
    chr1    800000  1450000

- ``-w/--window`` WINDOW

Window size used in TAD identification, with *RESOLUTION* as the unit. For
example, ``-w 2000`` and ``-R 10000`` generate a 20Mb window. Then this window
slides along each chromosome, and at each time, only a region of *WINDOW* size
is taken into account and subjected to clustering. (Default: 2000)

- ``-v/--version``

Print version number and exit.

- ``-h/--help``

Show help message and exit.

Next Steps
----------
That concludes the basic tutorial. It should be enough to get you up and
running our pipeline. However, if you want more details about the underlying
algorithms and the code, please carry on.

API Documentation
==================
.. toctree::
   :maxdepth: 3

   analyze
   polygon
   coniss

References
==========
.. [1] Imakaev M, Fudenberg G, McCord RP et al. Iterative correction of Hi-C data
   reveals hallmarks ofchromosome organization. Nat Methods, 2012, 9: 999-1003.

.. [2] Hu M, Deng K, Selvaraj S et al. HiCNorm: removing biases in Hi-C data via
   Poisson regression. Bioinformatics, 2012, 28: 3131-3.

.. [3] Yaffe E, Tanay A. Probabilistic modeling of Hi-C contact maps eliminates
   systematic biases to characterize global chromosomal architecture. Nat Genet,
   2011, 43: 1059-65.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

