# Created on Wed Oct 01 17:36:03 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

"""
An supplementary module for analyze.

To perform analysis, both interaction data and TADs are required.
:mod:`tadlib.analyze` module defines two classes :class:`tadlib.analyze.Inters`
and :class:`tadlib.analyze.TAD`, in which case you not only need to prepare
interaction data but also the matching TADs. This module provides a single
unified interface :class:`Queue` and identifies TADs according to your
interaction data automatically.

*CONISS*, a constrained hierarchical clustering algorithm, is used for TAD
identification. The cluster number is determined by a broken stick model.

Contents
--------
:class:`Queue`
    A customized **Inters** object. Identify TADs after loading interaction
    data.

Notes
-----
The core algorithms are implemented via two *R* packages *rioja* and *vegan*.
So *rpy2*, a python package, is utilized to access and embed *R* within our
code.

Examples
--------
Say your interaction data (resolution: 10000) file "h1hesc.chr1_chr1.txt"
is located in '~/data', **Queue** object can be created as follows:

>>> from tadlib import coniss
>>> path = '~/data'
>>> template = 'h1hesc.chr%s_chr%s.txt'
>>> workQueue = coniss.Queue(path = path, template = template)
>>> workQueue.TAD.dtype.names
('chr', 'start', 'end')

"""

import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, r
from analyze import Inters, getmatrix, manipulation
import logging, os

log = logging.getLogger(__name__)

class Queue(Inters):
    """Load interaction data and identify TADs.
    
    A customized **Inters** object. Pair-end reads from Hi-C should be
    processed into interactions between bins at a certain resolution. TADs
    can be identified by the new method :meth:`calTAD` using a constrained
    hierarchical clustering algorithm.
    
    TADs will be exported to a formatted text file, in which 3 columns
    indicating chromosome labels, start and end positions are contained.
    
    Parameters
    ----------
    path : str
        The path to the folder with source interaction data. (Default: '.',
        current working directory)
        
    Format : {'TXT', 'NPZ'}
        Two choices of the source format. (Default: 'TXT')
        
    template : str
        Template of the source file names. Only required for TXT sources.
        (Default: 'chr%s_chr%s.int', in which '%s' indicates chromosome
        ID)
        
    resolution : int
        Resolution of the interaction data. (Default: 10000, i.e., the bin
        size equals 10kb)
        
    chroms : list
        The list with the strings indicating which chromosome data will be
        read. Specially, '#' stands for chromosomes with numerical labels.
        If an empty list is provided, all chromosome data will be loaded.
        (Default: ['#', 'X'])
        
    cols : None or list of int (length: 3)
        Which columns to read, with 0 being the first. For example,
        ``cols = [1, 3, 4]`` will extract the 2nd, 4th and 5th columns.
        (Default: None, automatic recognization will be triggered)
        
    prefix : str or None
        Prefix of input .npz file name. For example, ``prefix = 'h1hesc'``
        will load a file named "h1hesc.npz". (Default: None)
        
    immortal : bool
        If True, save loaded data to a new Numpy .npz file.
        (Default: False)
    
    saveto : None or str
        Prefix of the output .npz file name. (Default: None)
    
    window : int
        Window size used in TAD identification, with **resolution** as the
        unit. ``window = 2000`` and ``resolution = 10000`` generate a
        20Mb window. It seems like a window sliding along the genome, and
        at each time, only a region of **window** size is taken into account
        and subjected to clustering. (Default: 2000)
    
    output : str or None
        Prefix of the generated TAD file name. (Default: None)
        
    Attributes
    ----------
    Some important attributes from **Inters**:
    
    location : str
        Location of both interaction data and corresponding TADs.
    
    resolution : int
        Resolution of the data.
    
    data : dict
        Loaded interaction data, with chromosome ID as the key and a
        customized Numpy Structured Array as the value. The structured
        array, which is defined through numpy.dtype object, has 3 fields
        named "bin1", "bin2" and "IF", respectively.
    
    labels : list of str
        List of sorted chromosome labels.
    
    Own attributes:
    
    window : int
        Window size for TAD identification.
    
    TAD : Structured Array
        Defined by numpy.dtype object. Contain 3 fields: "chr", "start"
        and "end".
    
    See Also
    --------
    numpy.dtype, tadlib.analyze.Inters
    
    """
    
    def __init__(self, path = '.', Format = 'TXT', resolution = 10000,
                 template = 'chr%s_chr%s.int', chroms = ['#', 'X'], cols = [],
                 prefix = None, immortal = False, saveto = None, window = 2000,
                 output = None):
        # Inherited from Inters
        Inters.__init__(self, path, Format = Format, template = template, 
                        resolution = resolution, chroms = chroms, cols = cols,
                        prefix = prefix, immortal = immortal, saveto = saveto)
        
        # Queue's own attributes
        self.window = window
        
        # Identify TADs
        self.calTAD()
        
        if output:
            log.debug('Export TADs ...')
            tadfile = '.'.join([str(output), 'domain', 'txt'])
            np.savetxt(os.path.join(self.location, tadfile), self.TAD,
                       fmt=['%s', '%d', '%d'], delimiter='\t')
            log.debug('TADs are saved in %s', 
                      os.path.join(self.location, tadfile))
    
    def calTAD(self):
        """Identify TADs using constrained hierarchical clustering.
        
        TADs, also known as self-associating domains, are continuous
        self-similar regions from the matrix perspective. Therefore,
        TAD identification essentially belongs to the clustering problem,
        in which clusters should be constrained by sample order.
        
        An algorithm called *CONISS* [1]_ is used for clustering. And the
        cluster number is determined with a broken stick model. [2]_
        
        See Also
        --------
        tadlib.analyze.getmatrix
        
        Notes
        -----
        *CONISS* has been implemented by a *R* package *rioja*, while broken
        stick model can be constructed by *bstick* function defined in
        another *R* package *vegan*.
        
        *CONISS* reads a whole interaction matrix into the memory, which is
        almost impossible for a large chromatin. So selecting a proper
        **window** size is subtle but important. Large **window** always
        contributes to better clusters but serious memory/time consumption.
        Relative small **window** uses less memory/time while losing accuracy
        to some extent.
        
        References
        ----------
        .. [1] Grimm EC. CONISS: A FORTRAN 77 program for stratigraphically
           constrained cluster analysis by the method of incremental sum of
           squares. Computers & Geosciences, 1987, 13: 13-35.
               
        .. [2] Bennett K. Determination of the number of zones in a
           biostratigraphic sequence. New Phytologist, 1996, 132: 155-170.
        
        """
        log.debug('Identifying TADs ...')
        # Create TAD export format
        strL = max([len(i) for i in self.labels])
        dtype = np.dtype({'names':['chr', 'start', 'end'],
                          'formats':['S'+str(strL), np.int, np.int]})
        # Temporary list for TAD repository
        pool = []
        # Import R packages
        # R elementary packages
        base = importr('base')
        stats = importr('stats')
        # Required by constrained hierarchical clustering
        vegan = importr('vegan')
        rioja = importr('rioja')
        
        rmatrix = r('matrix')
        # Loop over each chromosomes
        for i in self.labels:
            log.debug('Chromosome %s ...', i)
            current = self.data[i]
            # Chromosome Size
            cS = max(current['bin1'].max(), current['bin2'].max())
            # Window size
            step = min(self.window, cS-1)
            interval = range(0, cS, step)
            if interval[-1]!=cS:
                interval[-1] = cS
            B = [0] # Boundaries
            L = 0 # Left-side of the window
            # Window Shift
            for j in range(len(interval)-1):
                # Pre-defined window
                left = L; right = interval[j+1]
                matrix = getmatrix(current, left, right)
                mask = matrix.sum(axis=0)!=0
                # Consider those pool regions
                if mask.sum()>10:
                    # Remove vacant bins of the window
                    new, D = manipulation(matrix, left)
                    # Correlation Matrix for clustering
                    cor = np.corrcoef(new)
                    # Convert numpy array to R matrix object
                    rcor = rmatrix(numpy2ri.numpy2ri(cor), nrow=new.shape[0])
                    # CONISS-Constrained Clustering
                    clust = rioja.chclust(stats.dist(rcor, method='maximum'),
                                          method='coniss')
                    # Broken-Stick Model
                    model = vegan.bstick(clust, ng=new.shape[0]-1, plot=False)
                    # Choose a reasonable cluster number
                    num = np.array(model[0])
                    dispersion = np.array(model[1])
                    bstick = np.array(model[2])
                    idx = np.where(bstick > dispersion)[0][0]
                    N = num[idx] - 1
                    # Cut the tree
                    temp = list(base.cumsum(base.table(stats.cutree(clust,
                                                                    k=N))))
                    P = [D[m] for m in temp[:-1]]
                    P.append(D[temp[-1] - 1] + 1)
                    # Update boundaries
                    if P[-1] == right:
                        if j == (len(interval) - 2):
                            B.extend(P)
                        else:
                            B.extend(P[:-1])
                        if N > 1:
                            L = P[-2]
                        else:
                            L = left
                    else:
                        B.extend(P)
                        L = P[-1]
                else:
                    L = left
            # Update pool
            for j in xrange(len(B)-1):
                pool.append((i, B[j]*self.resolution, B[j+1]*self.resolution))
        self.TAD = np.array(pool, dtype = dtype)
