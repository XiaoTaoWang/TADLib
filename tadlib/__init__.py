# Created on Wed Oct 29 16:55:39 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

import xmlrpc.client
from pkg_resources import parse_version as V

__author__ = 'XiaoTao Wang'
__version__ = '0.4.2'
__license__ = 'GPLv3+'

## Check for update
try:
    pypi = xmlrpc.client.ServerProxy('http://pypi.python.org/pypi')
    available = pypi.package_releases('TADLib')
    if V(__version__) < V(available[0]):
        print('*'*75)
        print('Version {0} is out of date, Version {1} is available.'.format(__version__, available[0]))
        print()
        print('*'*75)
except:
    pass

