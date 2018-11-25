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
We have wrapped this pipeline into a user-friendly software called CALFEA
(CALculate FEAture for TADs).

Usage
^^^^^
``calfea [options]``

To run *calfea*, you need to provide a Hi-C matrix in `cool <https://github.com/mirnylab/cooler>`_
format and corresponding TAD list (a TXT file with 3 columns: *chrom*, *start* and *end*).

Depending on what data you already have, there are different tools you can choose to generate
*cool*:

- If you are starting from the beginning (FASTQ/SRA), I recommend using `runHiC <https://github.com/XiaoTaoWang/HiC_pipeline>`_,
  a user-friendly and efficient Hi-C data processing tool developed by our lab.

- If you are an old user of TADLib and have NPZ/TXT Hi-C matrix at hand, you can use the *toCooler*
  script distributed with another software of mine `hicpeaks <https://github.com/XiaoTaoWang/HiCPeaks>`_.

- In other case, try `cooler official tools <https://cooler.readthedocs.io/en/latest/cli.html>`_. 

The command looks like this::

    $ calfea -O test.txt -t tad_file.txt -p cool_uri --pw 2 --ww 5

As an example, we present most available parameters here.

- ``-O/--output`` OUTPUT

  Output file name. The output lines have 5 fields: *ChromID*, *TAD Start*,
  *TAD End*, *Aggregation Preference* and *Gap Ratio*. We trace *Gap Ratio* for
  each TAD because gap regions are always eliminated from original interaction
  matrix.

- ``-t/--tad-file`` TAD_FILE

  TAD source file name. The file must contain 3 columns,
  indicating *ChromID*, *TAD Start* and *TAD End*, respectively.

- ``-p/--path`` PATH

  Path to the *cool* URI. Note that URI is not equal to file path. Refer to
  the `cool schema <https://cooler.readthedocs.io/en/latest/schema.html>`_ for more details.

- ``--pw`` PW

  Width of the peak region. We use it in interaction selection and noise filtering
  procedure. (Default: 2)

- ``--ww`` WW

  Donut width. We use "donut" because the background of a peak looks like a donut
  in 2D contact matrix. (Default: 5)

After this command, two files **test.txt** and **calfea.log** are
created under current working directory. We use a rotating file
for logging. According to our settings, when the size of **calfea.log**
gets about 30K, it's closed and renamed to **calfea.log.1**, and at the same
time, a new file **calfea.log** is silently opened for output. In a word,
the system saves old log files by appending the extensions ".1", ".2" etc.,
and the file being written to is always **calfea.log**.

Other Options
^^^^^^^^^^^^^
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
