# Created on Thu Nov 13 15:47:14 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

"""
In this docstring, we will test if all modules are installed successfully.

Import necessary modules in advance:

>>> import os
>>> import numpy as np
>>> from tadlib import analyze, polygon

:mod:`analyze` first:

>>> path = os.path.join(os.path.dirname(analyze.__file__), 'data')
>>> source = os.path.join(path, 'gm12878.domain.txt')

>>> template = 'gm12878.chr%s_chr%s.txt'
>>> workInters = analyze.Inters(path = path, template = template)
>>> workInters.data.keys()
['1']
>>> len(workInters.data['1'])
60016

>>> workTAD = analyze.TAD(source)
>>> workTAD.data.shape
(3,)

>>> left = 229385000 / 5000 # Resolution: 5000
>>> right = 229665000 / 5000
>>> matrix = analyze.getmatrix(workInters.data['1'], left, right) # Chrom 1
>>> workCore = analyze.Core(matrix)
>>> len(matrix) == len(workCore.newM)
True
>>> workCore.longrange(pw = 4, ww = 7)
>>> len(workCore.pos)
61
>>> workCore.DBSCAN()
>>> workCore.Nc
4
>>> workCore.clusters['areas'].max()
83.5
>>> workCore.gdensity()
>>> print round(workCore.gden, 3)
0.638

>>> M = np.loadtxt(os.path.join(path, 'matrix.txt'))
>>> M[1,:] = 0; M[:,1] = 0
>>> newM, convert = analyze.manipulation(M)
>>> print newM # doctest: +NORMALIZE_WHITESPACE
[[ 0.3313  0.4975  0.0617]
 [ 0.8783  0.0286  0.0314]
 [ 0.0705  0.6353  0.3036]]
    
>>> convert == {0: 0, 1: 2, 2: 3}
True

Next, :mod:`polygon`:

Convex hull of the sample data set:

>>> points = np.loadtxt(os.path.join(path, 'points.txt'), dtype = int)
>>> P = polygon.Polygon(points)
>>> len(P.anchors)
9

>>> check = np.array([[9, 10], [6, 3], [10, 18]]) # The Queries
>>> P.isinside(check).round(4) # doctest: +NORMALIZE_WHITESPACE
array([ -4.2426, -10.    ,   0.5423])

Vertices of a square:

Clockwise:    

>>> sq = [(0,0), (0,1), (1,1), (1,0)]
>>> polygon.shoelace(sq)
-2.0
    
Counter-clockwise:

>>> sq = [(0,0), (1,0), (1,1), (0,1)]
>>> polygon.shoelace(sq)
2.0

>>> line = [(2, 0.4), (2, 0.8), (2, 4), (2, 100)]
>>> polygon.collinear(line)
True

"""

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose = True,
                    name = 'testall',
                    optionflags = doctest.REPORT_ONLY_FIRST_FAILURE)
