Install TADLib
==============

Requirements
------------
TADLib is developed and tested on UNIX-like operating system, and following Python
packages are required:

- Python 3.5+
- numpy
- scipy
- scikit-learn
- cooler
- pomegranate 0.10.0+
- networkx 1.x (not compatible with networkx 2.0 yet)

Install Conda
-------------
We strongly recommend using the conda package manager.

.. note:: If you have the Anaconda Distribution installed, you already have it.

Choose an appropriate `Miniconda installer <https://conda.io/miniconda.html>`_ for your system,
then in your terminal window type the following and follow the prompts on the installer screens::

    $ bash Miniconda3-latest-Linux-x86_64.sh

After that, update the environment variables to finish the Conda installation::

    $ source ~/.bashrc

Install the Required Packages
-----------------------------
First set up the channels to make all packages listed above accessible (note that the order is
important to guarantee the correct priority)::
    
    $ conda config --add channels conda-forge
    $ conda config --add channels defaults
    $ conda config --add channels r
    $ conda config --add channels bioconda

Then just type and execute the following command::
    
    $ conda install setuptools numpy scipy scikit-learn cooler pomegranate=0.10.0 networkx=1

Install TADLib
--------------
Finally, *TADLib* can be installed from PyPI by pip::

    $ pip install TADLib

TADLib has been installed successfully if no exception occurs in the above process.
