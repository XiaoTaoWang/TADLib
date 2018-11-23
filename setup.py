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
import os, sys, tadlib, glob
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if ((sys.version_info.major==2) and (sys.version_info.minor!=7)) or \
   ((sys.version_info.major==3) and (sys.version_info.minor<5)):
    print('PYTHON 2.7/3.5+ IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)

setuptools.setup(
    name = 'TADLib',
    version = tadlib.__version__,
    author = tadlib.__author__,
    author_email = 'wangxiaotao686@gmail.com',
    url = 'https://github.com/XiaoTaoWang/TADLib/',
    description = 'A Library to Explore Chromatin Interaction Patterns for Topologically Associating Domains',
    keywords = 'TAD Aggregation Preference AP sub-TAD hierarchy Hi-C cooler',
    long_description = read('README.rst'),
    long_description_content_type='text/x-rst',
    scripts = glob.glob('scripts/*'),
    packages = setuptools.find_packages(),
    classifiers = [
        'Programming Language :: Python',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
    )
