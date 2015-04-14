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
and function, it is important to develop a quantitative parameter to measure
the structural characteristics of TAD. TADLib is such a package to explore
the chromatin interaction pattern of TAD.

Inspired by the observation that there exist great differences in chromatin
interaction pattern and gene expression level among TADs, a chromatin interaction
feature called Aggregation Preference (AP) is developed to capture the aggregation
degree of long-range chromatin interactions. Application to human and mouse cell
lines (including both traditional Hi-C and in situ Hi-C datasets) shows that there
exist heterogeneous structures among TADs and the structural rearrangement across
cell lines is significantly associated with transcription activity remodelling.

TADLib package is written in Python and provides a three-step pipeline:

- Selecting long-range chromatin interactions in each TAD
- Finding aggregation patterns of selected interactions
- Calculating chromatin interaction feature of TAD

.. note:: By default, we suppose that the input Hi-C data are corrected appropriately.
   Otherwise, systematic biases in source data will negatively impact chromatin
   interaction selection and then parameter calculation. Several correction schemes
   are available online [1]_, [2]_, [3]_.

Links
=====
- `PyPI <https://pypi.python.org/pypi/TADLib>`_ (Download and Installation)
- `Repository <https://github.com/XiaoTaoWang/TADLib>`_ (At GitHub)

Modules
=======
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
- scikit-learn (>= 0.10)

.. note:: Tested systems: Red Hat Enterprise Linux Server release 6.4 (Santiago),
   CentOS release 6.4 (Final), Fedora release 20 (Heisenbug), Ubuntu 14.04 LTS

Installation
=============
Firstly, install all the required python packages:

There's an exhaustive instruction at http://www.scipy.org/install.html

We strongly recommend using `conda <http://conda.pydata.org/miniconda.html>`_,
an excellent Python package and environment manager.

Once Miniconda is installed, you can use the conda command to install any
other packages. Open a terminal and type::

    $ conda install numpy scipy scikit-learn

More details about conda can be found at http://conda.pydata.org/docs/

Then you can install TADLib just as other packages stored in PyPI:

Use *easy_install*::

    $ conda install setuptools
    $ easy_install TADLib

Or download the `source code <https://pypi.python.org/pypi/TADLib>`_ manually,
extract it and run the setup.py script::

    $ python setup.py install

Finally, run this command in a terminal prompt::

    $ python -m "tadlib.tests.testall"

If no exception occurs, congratulations, TADLib has been installed successfully!

QuickStart
==========
A sample data is also distributed with the source code. First, change to the
sample directory: (we suppose you are still in the distribution root directory)::

    $ cd lib/data

Open a Python interperter and load Hi-C data:

>>> from tadlib import analyze
>>> # Data happen to lie in the working directory
>>> workInters = analyze.Inters(path = '.', template = 'gm12878.chr%s_chr%s.txt')
>>> workInters.data['1'].shape # Chromosome 1
(60016,)

Load TAD regions:

>>> source = 'gm12878.domain.txt'
>>> workTAD = analyze.TAD(source)
>>> workTAD.data[0]
('1', 228780000, 229375000)

Consider such a TAD region: chr1:229385000-229665000

Initialize a *Core* object using interaction matrix:

>>> left = 229385000 / 5000 # Resolution 5000
>>> right = 229665000 / 5000
>>> matrix = analyze.getmatrix(workInters.data['1'], left, right)
>>> workCore = analyze.Core(matrix)

Extract long-range interactions:

>>> workCore.longrange(pw = 4, ww = 7)
>>> len(workCore.pos) # Number
62

Aggregation pattern discovery:

>>> workCore.DBSCAN()
>>> workCore.Nc # Cluster Number
4

Feature calculation:

>>> workCore.gdensity()
>>> print workCore.gden # Subjects to small variations across platforms
0.627528491366


CALFEA
======
We have wrapped this pipeline into a user-friendly software called CALFEA
(CALculate FEAture for TADs).

Usage
-----
``calfea [options]``

Example and Explanation
-----------------------
Use the sample data again. Create a new directory under the distribution root::

    $ mkdir working

Change current directory to "working"::

    $ cd working

And type in command below::

    $ calfea -f test1_feature.txt -s ../lib/data/gm12878.domain.txt -p ../lib/data -F TXT -R 5000 -T gm12878.chr%s_chr%s.txt -c 0 1 2 --immortal -S test1 --pw 4 --ww 7 --verbose 3

As an example, we present most available parameters here.

- ``-f/--filename`` FILENAME

Output file name. The output lines have 5 fields: *ChromID*, *TAD Start*,
*TAD End*, *Aggregation Preference* and *Gap Ratio*. We trace *Gap Ratio* for
each TAD because gap regions are always eliminated from original interaction
matrix. By default, results are printed to the standard output streams ("stdout").

- ``-s/--source`` SOURCE

Complete TAD source file name. The source file must contain 3 columns,
indicating *ChromID*, *TAD Start* and *TAD End*, respectively.

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

- ``--pw`` PW

Width of the peak region. We use it in interaction selection and noise filtering
procedure. (Default: 2)

- ``--ww`` WW

Donut size. We use "donut" because the background of a peak looks like a donut
in 2D contact matrix. (Default: 5)

- ``--verbose`` {0,1,2,3}

Logging level. 0: only show error messages, 1: also report warnings, 2:
show process information, 3: debug. We log all messages of all levels to a
disk file (calfea.log), while simultaneously logging process information or
above to the console. (Default: 2)

After this command, two files **test1_feature.txt** and **calfea.log** are
created under current working directory, and binary file named **test1.npz**
are created under the sample data folder "../lib/data". We use a rotating file
for logging. According to our settings, when the size of **calfea.log**
gets about 30K, it's closed and renamed to **calfea.log.1**, and at the same
time, a new file **calfea.log** is silently opened for output. In a word,
the system saves old log files by appending the extensions ".1", ".2" etc.,
and the file being written to is always **calfea.log**.

Other Options
-------------
- ``-C/--chroms`` [CHROMS [CHROMS ...]]

Chromosomes to read. Specially, "#" stands for chromosomes with
numerical labels. ``--chroms`` with zero argument will generate an empty
list, in which case all chromosome data will be loaded. (Default: ['#', 'X'])

- ``-P/--prefix`` PREFIX

Prefix of input *.npz* file name. Path should not be included. 

- ``--chromname`` CHROMNAME

Template of chromosome names contained in external TAD file. If you specify
``--chromosome chr``, the content should like this::

    chr1    0   350000
    chr1    350000  800000
    chr1    800000  1450000

- ``--top`` TOP

Parameter for noisy interaction filtering. By default, 30% noisy interactions
will be eliminated. (Default: 0.7)

- ``--ratio`` RATIO

Specifies the sample ratio of significant interactions for TAD. (Default: 0.05)

- ``--gap`` GAP

Maximum gap ratio of a TAD. (Default: 0.2)

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

