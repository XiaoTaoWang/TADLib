Introduction
------------
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

Installation
------------
Please check the file "INSTALL.rst" in the distribution.

Links
-----
- `Detailed Documentation <http://pythonhosted.org/TADLib/>`_
- `Repository <https://github.com/XiaoTaoWang/TADLib>`_ (At GitHub, Track the package issue)
- `PyPI <https://pypi.python.org/pypi/TADLib>`_ (Download and Installation)

Notes
-----
By default, we suppose that the input Hi-C data are corrected appropriately.
Otherwise, systematic biases in source data will negatively impact chromatin
interaction selection and then parameter calculation. Several correction schemes
are available online:

.. [1] Yaffe E, Tanay A. Probabilistic modeling of Hi-C contact maps eliminates
   systematic biases to characterize global chromosomal architecture. Nat Genet,
   2011, 43: 1059-65.

.. [2] Imakaev M, Fudenberg G, McCord RP et al. Iterative correction of Hi-C data
   reveals hallmarks ofchromosome organization. Nat Methods, 2012, 9: 999-1003.

.. [3] Hu M, Deng K, Selvaraj S et al. HiCNorm: removing biases in Hi-C data via
   Poisson regression. Bioinformatics, 2012, 28: 3131-3.
