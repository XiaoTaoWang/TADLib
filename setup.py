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

if (sys.version_info.major != 2) or (sys.version_info.minor != 7):
    print 'PYTHON 2.7 IS REQUIRED. YOU ARE CURRENTLY USING PYTHON ' + sys.version
    sys.exit(2)

# Guarantee Unix Format
for src in ['scripts/calfea','scripts/hitad']:
    text = open(src, 'rb').read().replace('\r\n', '\n')
    open(src, 'wb').write(text)

setup(
    name = 'TADLib',
    version = lib.__version__,
    author = lib.__author__,
    author_email = 'wangxiaotao868@163.com',
    url = 'https://github.com/XiaoTaoWang/TADLib/',
    description = 'A Library to Explore Chromatin Interaction Patterns for Topologically Associating Domains',
    keywords = 'TAD Aggregation Preference AP sub-TAD hierarchy Hi-C',
    package_dir = {'tadlib':'lib'},
    packages = ['tadlib', 'tadlib.calfea', 'tadlib.hitad'],
    scripts = ['scripts/calfea', 'scripts/hitad'],
    package_data = {'tadlib.calfea':['data/*'],
                    'tadlib.hitad':['data/NPZ/*',
                                    'data/TXT/datasets*',
                                    'data/TXT/IMR90-HindIII-rep1/*',
                                    'data/TXT/IMR90-HindIII-rep2/*']},
    long_description = read('README.rst'),
    classifiers = [
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Mathematics'
        ]
    )
