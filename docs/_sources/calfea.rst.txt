Aggregation Preference
======================

Introduction
------------
Inspired by the observation that there exist great differences in chromatin interaction
pattern among TADs, we proposed an empirical parameter called Aggregation Preference(AP)
to measure the overall aggregation degree of significant chromatin interactions inside
TAD. Application to human and mouse cell types (including both traditional Hi-C and in situ
Hi-C data sets) shows that there exist heterogeneous structures among TADs and the structural
rearrangement across cell types is significantly associated with transcriptional remodelling.

Generally, it takes 3 steps to calculate the AP value:

1. Select long-range significant chromatin interactions in each TAD
2. Find aggregation patterns of selected interactions by using a density-based clustering
   algorithm called DBSCAN
3. Calculate the AP value of TAD


QuickStart
----------
A sample data is also distributed with the source code. First, change to the
sample directory: (we suppose you are still in the distribution root directory)::

    $ cd lib/calfea/data

Open a Python interperter and load Hi-C data:

>>> from tadlib.calfea import analyze
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

AP calculation:

>>> workCore.gdensity()
>>> print workCore.gden # Subjects to small variations across platforms
0.627528491366


Guide for CALFEA
----------------
We have wrapped this pipeline into a user-friendly software called CALFEA
(CALculate FEAture for TADs).

Usage
^^^^^
``calfea [options]``

Example and Explanation
^^^^^^^^^^^^^^^^^^^^^^^
Use the sample data again. Create a new directory under the distribution root::

    $ mkdir working

Change current directory to "working"::

    $ cd working

And type in command below::

    $ calfea -f test1_feature.txt -s ../lib/calfea/data/gm12878.domain.txt -p ../lib/calfea/data -F TXT -R 5000 -T gm12878.chr%s_chr%s.txt -c 0 1 2 --immortal -S test1 --pw 4 --ww 7

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

  3 columns reading from TXT Hi-C data, with 0 being the first. Here,
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

After this command, two files **test1_feature.txt** and **calfea.log** are
created under current working directory, and binary file named **test1.npz**
are created under the sample data folder "../lib/calfea/data". We use a rotating file
for logging. According to our settings, when the size of **calfea.log**
gets about 30K, it's closed and renamed to **calfea.log.1**, and at the same
time, a new file **calfea.log** is silently opened for output. In a word,
the system saves old log files by appending the extensions ".1", ".2" etc.,
and the file being written to is always **calfea.log**.

Other Options
^^^^^^^^^^^^^
- ``-C/--chroms`` [CHROMS [CHROMS ...]]

  Chromosomes to read. Specially, "#" stands for chromosomes with
  numerical labels. ``--chroms`` with zero argument will generate an empty
  list, in which case all chromosome data will be loaded. (Default: ['#', 'X'])

- ``-P/--prefix`` PREFIX

  Prefix of input *.npz* file name. Path should not be included.

- ``--chromname`` CHROMNAME

  Leading string of chromosome names contained in TAD source file. If you specify
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
^^^^^^^^^^
That concludes the basic tutorial. It should be enough to get you up and
running our pipeline. However, if you want more details about the underlying
algorithms and the code, please carry on.

API Documentation
-----------------
API reference of our defined classes and functions for Aggregation Preference(AP)
calculation.

.. toctree::
   :maxdepth: 2

   calfea_api
