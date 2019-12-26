TAD visulization
****************
Since 0.4.2, I've included a visulization module, mainly designed to customize triangular rotated
heatmap and plot TADs/loops on it.

Main features:

1. `cool <https://github.com/mirnylab/cooler>`_ file support
2. extract/manipulate coordinates of each contact in triangle shape


Example
^^^^^^^
Here are some screenshots from [1]_ by using this module:

.. image:: ./_static/HiTAD.supple-fig14.png
        :align: center


Basic Tutorial
^^^^^^^^^^^^^^
First, let's change to the ``data`` folder distributed with ``TADLib``::

    $ pwd
    
    /Users/xtwang/workspace/MyPackages/TADLib/data

Then follow the commands below within a Python interpreter::

    In [1]: from tadlib.visualize.heatmaps import *
    In [2]: vis = Triangle('GM12878-MboI.hg19.cool::10000', '21', 45000000, 46500000)
    In [3]: vis.matrix_plot()
    In [4]: vis.plot_loops('GM12878-MboI.10K.peakachu.bedpe')
    In [5]: vis.plot_TAD('GM12878-MboI.20K.HiTAD.txt', linewidth=1.5)
    In [6]: vis.show()

.. image:: ./_static/HiTAD.visualize.demo.png
        :align: center



Reference
^^^^^^^^^
.. [1] Wang XT, Cui W, Peng C. HiTAD: detecting the structural and functional hierarchies of
   topologically associating domains from chromatin interactions. Nucleic Acids Research, 2017,
   doi: 10.1093/nar/gkx735






