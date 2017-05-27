Hierarchical TAD
****************

Introduction
============
HiTAD is a method to detect hierarchical TADs, including TADs, sub-TADs and smaller domains.
Except local insulations, HiTAD further constrains TADs as the optimal domains to globally
separate intra-chromosomal interactions. Under objective functions derived from chromatin
interactions, HiTAD adopts an iterative optimization procedure to detect hierarchical TADs.
HiTAD performs well in domain sensitivity, replicate reproducibility and inter cell-type
conservation. Application to human and mouse cell types (including both traditional Hi-C
and in situ Hi-C data sets) reveals that there exist common change types for hierarchical
TADs, which are involved in shaping higher-order compartment, replication timing and
transcriptional regulation.

Snapshot
========
This tutorial will cover the two most commonly used tools/submodules in HiTAD:

- hitad

  A command-line tool streamlining our 5-step identification algorithm:

  1. Calculate adaptive directionality index (DI) for each bin.
  2. Detect bottom boundaries by 5-state Gaussian mixture Hidden Markov Model using adaptive
     DIs as input.
  3. Identify TAD from bottom domains by optimally separating intra-chromosomal interactions
     under chromosome-level objective functions.
  4. Recursively identify inner domains under domain-level objective functions to optimally
     separate intra-domain interactions.
  5. (Optional but recommended) Perform domain alignment between two replicates and only
     the reproducible domains are maintained to generate the final hierarchical TAD list.

  .. note:: If the two replicates are comparable in sequencing depth, then step 5 can improve
     the accuracy and reliability while guaranteeing the sensitivity of called domains; otherwise,
     you'd better merge the replicates into one dataset and step 5 will be skipped in this case.

- aligner

  A submodule containing classes and functions for our proposed domain alignment strategy. In
  our work, this strategy is used in reproducible domain detecting and change type defining
  at both boundary and domain level.

Tutorial
========

Command-line tool
-----------------

Usage
^^^^^
``hitad [options]``

Data Preparation
^^^^^^^^^^^^^^^^^
You need to perform a standard Hi-C data processing pipeline (mapping, filtering, binning, correcting)
and get the bin-level contacts for HiTAD.

Only intra-chromosomal contacts are allowed and they should be stored chromosome by chromosome in
text files. These files should be named following the pattern ``chr*.txt`` and be placed under the
same folder, where "*" indicates the chromosome labels (1,2,...,X,Y for human). Each file should
contain 3 columns (separated by whitespace)::

    <bin1> <bin2> <IF>

- bin1, bin2: bin index on the corresponding chromosome
- IF : normalized/corrected interaction frequency

Learn HiTAD Using the TXT Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First, change to the example directory and see what the inputs should look like (suppose you are still
in the distribution root directory)::

    $ cd lib/hitad/data
    $ ls -lh

There are two sub-directories named "NPZ" and "TXT", corresponding to two supported data formats
respectively::

    total 16K
    drwxr-x--- 2 xtwang CPeng 4.0K May  3 13:36 NPZ
    drwxr-x--- 4 xtwang CPeng 4.0K May  3 13:36 TXT

"NPZ" is a binary format we define to accelerate contact data loading. You'll know it in more details
soon. For now, just pay attention to the "TXT"::

    $ cd TXT
    $ ls -lh *

There are two sub-directories and two files::

    -rw-r----- 1 xtwang CPeng   37 May  3 13:36 datasets-onerep
    -rw-r----- 1 xtwang CPeng   65 May  3 13:36 datasets-tworeps

    IMR90-HindIII-rep1:
    total 2.7M
    -rw-r----- 1 xtwang CPeng 2.7M May  3 13:36 chr22.txt
    
    IMR90-HindIII-rep2:
    total 2.8M
    -rw-r----- 1 xtwang CPeng 2.8M May  3 13:36 chr22.txt

The two sub-directories correspond to contact data of two replicates, respectively. And the two
files are metadata files specifying the resolution (in base-pair unit), replicate label and data path
(in TXT case, it should point to the folders) information in a structured way::

    $ cat datasets-tworeps

    res:40000
      rep1:../lib/hitad/data/TXT/IMR90-HindIII-rep1
      rep2:../lib/hitad/data/TXT/IMR90-HindIII-rep2

    $ cat datasets-onerep

    res:40000
      rep1:../lib/hitad/data/TXT/IMR90-HindIII-rep1

If you want to use two replicates to improve reliability of identified domains (as we did in our work),
you should provide hitad with a metadata file like "datasets-tworeps"; if you have only one replicate
or the sequencing depths between two replicates are quite different, your metadata file should look
like "datasets-onerep".

With contact data and metadata file ready, running hitad is quite straightforward.

Go back to the distribution root and create a new directory there for working::

    $ cd ../../../..
    $ mkdir working

Change to "working" directory::

    $ cd working

Then type in the command below::

    $ hitad -O IMR90-HindIII-40K-tworeps.txt -d ../lib/hitad/data/TXT/datasets-tworeps --logFile hitad-40K.log

- ``-O/--output`` OUTPUT

  Output file name.

- ``-d/--datasets`` DATASETS

  Metadata file path describing your Hi-C data

- ``--logFile`` LOGFILE

  Log file name. (Default: hitad.log)

After this command(it takes about several seconds), two files "IMR90-HindIII-40K-tworeps.txt"
and "hitad-40K.log" will be generated under current working directory. We use a rotating file
for logging. According to our settings, when the size of "hitad-40K.log" gets about 200K, it's
closed and renamed to "hitad-40K.log.1". At the same time, a new file "hitad-40K.log" is silently
opened for output. In a word, the system saves old log files by appending the extensions
".1", ".2" etc., and the current log is always written to "hitad-40K.log"::

    $ head -6 IMR90-HindIII-40K-tworeps.txt
    
    22      18240000        18560000        0
    22      19000000        19640000        0
    22      19920000        20120000        0
    22      20800000        21480000        0
    22      22040000        22320000        0
    22      22320000        22720000        1

The output file contains 4 columns indicating chromosome label, domain start (bp), domain end (bp),
and hierarchical label, respectively. In our notation, TAD is denoted as level 0, sub-TAD is denoted
as level 1, and subsequent domain level is denoted as level 2, etc.

Transform TXT into NPZ
^^^^^^^^^^^^^^^^^^^^^^^
First, let's take a look at a NPZ file within Python environment::

    $ cd ../lib/hitad/data/NPZ

Open a Python interperter and type:

>>> import numpy as np
>>> lib = np.load('GM12878-HindIII.40K.rep1.npz')
>>> lib.keys()
['22']
>>> type(lib['22'])
<type 'numpy.ndarray'>
>>> lib['22'][:5]
array([(401, 686, 0.149032), (401, 693, 0.175786), (401, 695, 0.122059),
       (401, 868, 0.156891), (401, 1218, 0.211272)], 
      dtype=[('bin1', '<i8'), ('bin2', '<i8'), ('IF', '<f8')])

A NPZ file is a pool of arrays using chromosome labels as the keys. Each array is structured and has
3 fields: "bin1", "bin2" and "IF", just as the information contained in a TXT file mentioned
above.

NPZ files can be generated in two ways:

1. Run hitad with ``--npzpre`` specified.
2. Use `runHiC <https://github.com/XiaoTaoWang/HiC_pipeline>`_ pipeline. runHiC is a user-friendly
   command-line software developed by our lab for Hi-C data processing.

Let's rerun hitad with ``--npzpre`` and see what will happen::

    $ cd ../../../../working/
    $ hitad -O IMR90-HindIII-40K-tworeps.txt -d ../lib/hitad/data/TXT/datasets-tworeps --npzpre IMR90-HindIII --logFile hitad-40K.log

- ``--npzpre`` NPZPRE

  Prefix of output .npz file name. If no path is contained, NPZ files will be placed under
  current working directory.

Two additional files "IMR90-HindIII.40K.rep1.npz" and "IMR90-HindIII.40K.rep2.npz" will
be generated under current working directory for two replicates respectively.

Run HiTAD with NPZ Data
^^^^^^^^^^^^^^^^^^^^^^^
A metadata file is also required in the case of NPZ data::

    $ cd ../lib/hitad/data/NPZ
    $ cat datasets-tworeps
    
    res:40000
      rep1:../lib/hitad/data/NPZ/GM12878-HindIII.40K.rep1.npz
      rep2:../lib/hitad/data/NPZ/GM12878-HindIII.40K.rep2.npz

.. note:: Data path should point to the npz file.

Execute following commands and try to understand the result by yourself::

    $ cd ../../../../working/
    $ hitad -O GM12878-HindIII-40K-tworeps.txt -d ../lib/hitad/data/NPZ/datasets-tworeps --logFile hitad-40K.log

Other Options
^^^^^^^^^^^^^
- ``-C/--chroms`` [CHROMS [CHROMS ...]]
  
  List of chromosome labels. Only Hi-C data within the specified chromosomes will be included.
  Specially, "#" stands for chromosomes with numerical labels. The default setting is ``--chrom "#" X``,
  which will include 1, 2, ..., 22 and X chromosome for human genome. (Default: ['#', 'X'])

- ``--maxsize`` MAXSIZE

  Maximum domain size in base-pair unit. (Default: 4000000)

- ``--cache`` CACHE

  Cache folder path. HiTAD will try to create the folder and all intermediate
  ones if they don't exist. (Default: /tmp)

- ``--removeCache``

  Remove cached data before existing. (Default: False)

- ``-p/--cpu-core`` CPU_CORE

  Number of processes to launch. (default: 1)

- ``-v/--version``

  Print version number and exit.

- ``-h/--help``

  Show help message and exit.


Domain Alignment
----------------
Traditionally, domains and boundaries are aligned by matching boundaries with nearest ones
between two data sets. This strategy generally assign a threshold in advance to determine
whether a boundary of one data set exists in the other data set. However, the selection
of the threshold is quite empirical and artificial. To deal with this problem, we propose
a parameter-free alignment strategy by maximizing the overlap ratio between matched domains
and considering all domains at the same time. We also generalize the strategy on hierarchical
TADs and further define change types at both boundary and domain levels systematically.

All related classes and functions are defined in :py:mod:`tadlib.hitad.aligner`.

In this section, I will show you how to perform alignment/comparison between two
domain sets.

During the previous sections, we have generated 2 domain sets, "GM12878-HindIII-40K-tworeps.txt"
and "IMR90-HindIII-40K-tworeps.txt", under the "working" directory.

Open a Python interperter and load these domain lists by
:py:func:`tadlib.hitad.aligner.readHierDomain`:

>>> from tadlib.hitad.aligner import *
>>> ilist = readHierDomain('IMR90-HindIII-40K-tworeps.txt')
>>> glist = readHierDomain('GM12878-HindIII-40K-tworeps.txt')
>>> print len(ilist), len(glist)
80 75
>>> print ilist[:3]
[['22', 18240000, 18560000, 0], ['22', 19000000, 19640000, 0], ['22', 19920000, 20120000, 0]]

We also provide a function :py:func:`tadlib.hitad.aligner.hierFormat` to transform
[chrom,start,end] domains to [chrom,start,end,level] domains. That will be helpful if
you want to align domain sets generated by other domain callers (such as `Arrowhead <https://github.com/theaidenlab/juicer/wiki/Arrowhead>`_
and `TADtree <http://compbio.cs.brown.edu/projects/tadtree/>`_).

>>> pseudo = [d[:3] for d in ilist]
>>> print pseudo[:3]
[['22', 18240000, 18560000], ['22', 19000000, 19640000], ['22', 19920000, 20120000]]
>>> transformed = hierFormat(pseudo)
>>> transformed == ilist
True

Hold domain lists with :py:class:`tadlib.hitad.aligner.DomainSet`:

>>> iset = DomainSet('IMR90', ilist, 40000)
>>> gset = DomainSet('GM12878', glist, 40000)

To construct hierarchical domain mapping between two domain sets, we need
to initialize a :py:class:`tadlib.hitad.aligner.DomainAligner` and invoke the
:py:meth:`tadlib.hitad.aligner.DomainAligner.align` method:

>>> work = DomainAligner(iset, gset)
>>> work.align('IMR90','GM12878')

:py:class:`tadlib.hitad.aligner.DomainAligner` also defines domain-level change types, including
conserved TADs, semi-conserved TADs, merged TADs, split TADs, and unaligned TADs:

>>> conserved = work.conserved('IMR90','GM12878') # Conserved TADs
>>> print len(conserved)
20
>>> print sorted(conserved)[0] # The first conserved TAD pair
(('22', 19920000, 20120000), ('22', 19920000, 20120000))

>>> semi = work.inner_changed('IMR90','GM12878') # Semi-Conserved TADs
>>> print len(semi)
4
>>> print sorted(semi)[0] # The first semi-conserved TAD pair
(('22', 27080000, 28240000), ('22', 27080000, 28200000))

>>> merged = work.merged('IMR90','GM12878') # Merged TADs
>>> print sorted(merged)[0] # The first merged region
(('22', 31720000, 32360000), ('22', 31720000, 32400000))
>>> merged[(('22', 31720000, 32360000), ('22', 31720000, 32400000))] # Merged details
[[['22', 31720000, 32040000, 0], ['22', 32040000, 32360000, 0]], [['22', 31720000, 32400000, 0]]]

>>> split = work.split('IMR90','GM12878') # Split TADs
>>> print sorted(split)[0] # The first split region
(('22', 19000000, 19640000), ('22', 19160000, 19720000))
>>> split[(('22', 19000000, 19640000), ('22', 19160000, 19720000))] # Split details
[[['22', 19000000, 19640000, 0]], [['22', 19160000, 19440000, 0], ['22', 19440000, 19720000, 0]]]

Boundary-level change types are defined in :py:class:`tadlib.hitad.aligner.BoundAligner` which
is also based on our domain-based alignment algorithm:

>>> boundview = BoundAligner(iset, gset)
>>> # Pairs of conserved TAD boundaries
>>> conserved_bounds = boundview.conserved_tad_bounds('IMR90','GM12878')
>>> print conserved_bounds.items()[:3]
[(('22', 43560000), ('22', 43560000)), (('22', 45600000), ('22', 45520000)), (('22', 46720000), ('22', 46920000))]

>>> # Pairs of conserved sub-TAD boundaries.
>>> conserved_subs = boundview.conserved_sub_bounds('IMR90','GM12878')
>>> print conserved_subs.items()[:3]
[(('22', 36000000), ('22', 36000000)), (('22', 34320000), ('22', 34320000)), (('22', 30480000), ('22', 30480000))]

>>> # TAD to sub-TAD switch cases
>>> tad_sub = boundview.tad2sub('IMR90','GM12878')
>>> len(tad_sub)
1

>>> # sub-TAD to TAD switch cases.
>>> sub_tad = boundview.sub2tad('IMR90','GM12878')
>>> len(sub_tad)
6

>>> # TAD boundaries that exist in IMR90, but disappear in GM12878.
>>> disappear_TAD = boundview.disappeared_tad('IMR90','GM12878')
>>> sorted(disappear_TAD)[0] # The first disappeared TAD boundary
('22', 18240000)

>>> # Sub-TAD boundaries that exist in IMR90, but disappear in GM12878.
>>> disappear_sub = boundview.disappeared_sub('IMR90','GM12878')
>>> sorted(disappear_sub)[0] # The first disappeared sub-TAD boundary
('22', 22720000)

API Documentation
=================
API reference of our defined classes and functions for hierarchical TAD
identification and domain alignment.

.. toctree::
   :maxdepth: 2

   hitad_api
