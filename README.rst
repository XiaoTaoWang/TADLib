.. note:: Since version 0.4.0, the default data format has changed to `cool <https://github.com/mirnylab/cooler>`_,
   to comply with 4DN standards.


Introduction
============
Chromosome conformation capture (3C) derived techniques, especially Hi-C, have
revealed that topologically associating domain (TAD) is a structural basis for
both chromatin organization and biological functions in three-dimensional (3D)
space. TAD is also hierarchically organized by smaller structural units, which
are relevant to biological functions. To systematically investigate the relationship
between structure and function, it is important to develop quantitative methods to
identify and measure the organization of TAD. TADLib is such a library to explore
the chromatin interaction patterns inside TAD from Hi-C chromatin interactions.

Currently, TADLib consists of three methods:

- Aggregation Preference (AP)
    AP is a quantitative parameter to measure the overall density of significant
    chromatin interactions inside TAD. Inspired by the observation that there exist
    great differences in chromatin interaction pattern among TADs, an empirical
    parameter called Aggregation Preference (AP) can be used to capture these
    aggregation degree of significant chromatin interactions. Application to human
    and mouse cell types (including both traditional Hi-C and in situ Hi-C data sets)
    shows that there exist heterogeneous structures among TADs and the structural
    rearrangement across cell types is significantly associated with transcriptional
    remodelling. [1]_
- Hierarchical TAD (HiTAD)
    HiTAD is a method to detect hierarchical TADs, including TADs, sub-TADs and
    smaller domains. Except local insulations, HiTAD further constrains TADs as the
    optimal domains to globally separate intra-chromosomal interactions. Under
    objective functions derived from chromatin interactions, HiTAD adopts an iterative
    optimization procedure to detect hierarchical TADs. HiTAD performs well in domain
    sensitivity, replicate reproducibility and inter cell-type conservation. Application
    to human and mouse cell types (including both traditional Hi-C and in situ Hi-C data
    sets) reveals that there exist common change types for hierarchical TADs, which are
    involved in shaping higher-order compartment, replication timing and transcriptional
    regulation. [2]_
- Domain Caller (domaincaller)
    Domaincaller is an implementation of the original Directionality Index (DI) based
    TAD caller proposed by Dixon et al. [3]_ Instead of original separate scripts for
    each stage of the caller, this module provides a convenient command line interface
    integrating the whole pipeline, including calculating the DI track, performing HMM
    and post-processing. It supports the `.cool <https://github.com/mirnylab/cooler>`_
    matrix format, so has low memory requirements when dealing with high resolution data.
    The TAD calls on the `Human ES dataset <http://chromosome.sdsc.edu/mouse/hi-c/hESC.tar.gz>`_
    are pretty similar to the published list, both in boundary locations and CTCF enrichment.
    

Links
=====
- `Detailed Documentation <https://xiaotaowang.github.io/TADLib/>`_
    - `Installation <https://xiaotaowang.github.io/TADLib/install.html>`_
    - `Aggregation Preference <https://xiaotaowang.github.io/TADLib/calfea.html>`_
    - `Hierarchical TAD <https://xiaotaowang.github.io/TADLib/hitad.html>`_
- `Code Repository <https://github.com/XiaoTaoWang/TADLib>`_ (At GitHub, Track the package issue)
- `PyPI <https://pypi.python.org/pypi/TADLib>`_ (Download and Installation)
	
Citation
========
.. [1] Wang XT, Dong PF, Zhang HY, Peng C. Structural heterogeneity and functional diversity
   of topologically associating domains in mammalian genomes. Nucleic Acids Research, 2015,
   doi: 10.1093/nar/gkv684

.. [2] Wang XT, Cui W, Peng C. HiTAD: detecting the structural and functional hierarchies of
   topologically associating domains from chromatin interactions. Nucleic Acids Research, 2017,
   doi: 10.1093/nar/gkx735

.. [3] Dixon JR, Selvaraj S, Yue F, Kim A, Li Y, Shen Y, Hu M, Liu JS, Ren B. Topological domains
   in mammalian genomes identified by analysis of chromatin interactions. Nature, 2012,
   doi: 10.1038/nature11082
