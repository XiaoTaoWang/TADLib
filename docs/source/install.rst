Install TADLib
==============

Requirements
------------
TADLib is developed and tested on UNIX-like operating system, and following Python
packages are required:

- Python (2.7, not compatible with 3.x for now)
- Numpy (>= 1.11)
- Scipy (>= 0.18)
- Scikit-Learn (>= 0.18)
- GHMM (>= 0.9)

Install Conda
-------------
We strongly recommend using the conda package manager.

.. note:: If you have the Anaconda Distribution installed, you already have it.

Download the latest `Linux Miniconda installer for Python 2.7 <https://conda.io/miniconda.html>`_,
then in your terminal window type the following and follow the prompts on the installer screens::

    $ bash Miniconda2-latest-Linux-x86_64.sh

After that, update the environment variables to finish the Conda installation::

    $ source ~/.bashrc

Install the Required Packages
-----------------------------
Conda allows separation of packages into separate repositories, or channels. The main *defaults*
channel has a large amount of common packages including *Numpy*, *Scipy*, and *Scikit-Learn* listed
above. To install these packages, just type and execute the following command::

    $ conda install setuptools numpy scipy scikit-learn


Then it's straightforward to install all the required packages through the following one-line command::

    $ conda install setuptools numpy scipy scikit-learn ghmm 

Set up Channels
---------------
*GHMM* is not available in the *defaults* channel but included in the *bioconda* channel, and
to make it accessible, you need to add the *bioconda* channel as well as the other channels bioconda
depends on(note that the order is important to guarantee the correct priority)::

    $ conda config --add channels conda-forge
    $ conda config --add channels defaults
    $ conda config --add channels r
    $ conda config --add channels bioconda

Then you can install *GHMM* by::

    $ conda install ghmm

Install TADLib
--------------
Download the `source code <https://pypi.python.org/pypi/TADLib>`_ of TADLib, extract it and run
the setup.py script::

    $ python setup.py install

TADLib has been installed successfully if no exception occurs in the above process.
