Installation Guide for TADLib
==============================

Requirements
============
TADLib is developed and tested on UNIX-like operating system, and following Python
packages are recommended:

- Python (2.7, not compatible with 3.x)
- Numpy (>= 1.6)
- Scipy library (>= 0.13)
- scikit-learn (>= 0.11)

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
