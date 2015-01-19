# Created on Wed Nov 12 14:01:48 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

"""
Setup script for TADLib (A Library to Explore Chromatin Interaction Patterns
for Topologically Associating Domains).

This is a free software under GPLv3. Therefore, you can modify, redistribute
or even mix it with other GPL-compatible codes. See the file LICENSE
included with the distribution for more details.

"""
import os, sys, lib
from distutils.core import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major != 2) or (sys.version_info.minor < 6):
    print 'PYTHON VERSION MUST BE 2.6 or 2.7. YOU ARE CURRENTLY USING PYTHON ' + sys.version
    sys.exit(2)

# Guarantee Unix Format
text = open('scripts/calfea', 'rb').read().replace('\r\n', '\n')
open('scripts/calfea', 'wb').write(text)

setup(
    name = 'TADLib',
    version = lib.__version__,
    author = lib.__author__,
    author_email = 'wangxiaotao868@gmail.com',
    url = 'https://github.com/XiaoTaoWang/TADLib/',
    description = 'A Library to Explore Chromatin Interaction Patterns for Topologically Associating Domains',
    keywords = 'chromosome structure feature Hi-C TAD CONISS polygon',
    package_dir = {'tadlib':'lib'},
    packages = ['tadlib', 'tadlib.tests'],
    scripts = ['scripts/calfea'],
    package_data = {'tadlib':['data/*']},
    long_description = read('README.rst'),
    classifiers = [
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Mathematics'
        ]
    )
