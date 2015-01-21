# Created on Sat Sep 27 16:18:27 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

"""
Interaction selection, pattern discovery and feature calculation.

Contents
--------
:class:`Inters`
    Class for loading interaction data.

:class:`TAD`
    Class for loading TAD (Topologically Associating Domain) data.

:class:`Core`
    Define and calculate interaction feature for TAD.

:func:`getmatrix`
    Function for creating interaction matrix from interaction data.

:func:`manipulation`
    Deal with vacant rows and columns in interaction matrix.

:func:`extract_kth_diag`
    Return coordinates of a specified diagonal. A complementary function
    for numpy matrix operations.

:func:`axis_control`
    A customized axis control operations for matplotlib.

:func:`isnumber`
    Judge if a string is numerical or not.

Architecture
------------
The core algorithm is implemented by **Core**, whose only parameter, an
interaction matrix can be created by **getmatrix**, inputs of which are
always from **Inters** and **TAD**.

Above relations are shown in schema as follows:

.. image:: _static/ObjectRelationship.png

Examples
--------
Say your interaction data file "h1hesc.chr1_chr1.txt" and TAD file
"h1hesc.domain.txt" at both located in "~/data/", then we have:

>>> from tadlib import analyze
>>> path = '~/data/'
>>> source = '~/data/h1hesc.domain.txt'
>>> template = 'h1hesc.chr%s_chr%s.txt'
>>> workInters = analyze.Inters(path = path, template = template)
>>> workInters.data.keys()
['1']
>>> len(workInters.data['1'])
10535

>>> workTAD = analyze.TAD(source)
>>> workTAD.data.shape
(7,)

Suppose a TAD, [4670000, 5570000) on chromosome 1, then you can create a
Core object like this:

>>> left = 6490000 / 10000 # Resolution: 10000
>>> right = 6770000 / 10000
>>> matrix = analyze.getmatrix(workInters.data['1'], left, right) # Chrom 1
>>> workCore = analyze.Core(matrix)
>>> len(matrix) == len(workCore.newM)
False
>>> workCore.longrange()
>>> len(workCore.pos)
9
>>> workCore.DBSCAN()
>>> workCore.Nc
1
>>> workCore.clusters['areas']
array([ 2.5])
>>> workCore.gdensity()
>>> print round(workCore.gden, 4)
0.3563

"""

from __future__ import division
import glob, re, os, sys
import logging
import warnings
import polygon
import numpy as np
from scipy.stats import itemfreq, poisson
from scipy.spatial import distance
from scipy.interpolate import interp1d, splrep, splev
import scipy.special as special
from sklearn import cluster

import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
# Other basic settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])

log = logging.getLogger(__name__)

class Inters(object):
    """
    Load interaction data from TXT or Numpy .npz file.
    
    Here, interaction data are produced by Hi-C [1]_, which is a high-
    throughput variant of 3C [2]_ and is widely used to study chromatin
    interactions across an entire genome.
    
    You need to perform some conventional processes, such as read alignment
    and filtering in advance. Technically, our pipeline accepts binned data,
    and only within-chromosome ones are allowed.
    
    Although not required, some correction procedures [3]_ are recommended
    to eliminate systematic biases.
    
    If source data are provided in TXT format, we suppose they are split
    and saved chromosome by chromosome, and the files' name should contain
    chromosome information. 3 columns are required in these files. In order,
    they are bin1, bin2, and corresponding frequency.
    
    Data can also be stored in Numpy .npz format, with chromosome ID as
    the key. Such file type makes I/O extremely fast.
    
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
        
    saveto : str or None
        Prefix of the output .npz file name. (Default: None)
        
    Attributes
    ----------
    location : str
        Absolute path to the data folder. Output .npz file will also be
        exported here.
        
    chroms : set
        Set of the parameter chroms.
        
    template : str
        Template of the source file names.
        
    resolution : int
        Resolution of the data.
        
    cols : None or list of int (length: 3)
        Which columns in TXT source file to read.
        
    data : dict
        Loaded interaction data, with chromosome ID as the key and a
        customized Numpy Structured Array as the value. The structured
        array, which is defined through numpy.dtype object, has 3 fields
        named "bin1", "bin2" and "IF", respectively.
        
    interFiles : list of str
        List of TXT source files. Sorted according to corresponding
        chromosome labels. Only created if ``Format = 'TXT'``.
        
    labels : list of str
        List of sorted chromosome labels.
    
    label2idx : dict
        Map from chromosomes labels to zero-based indices.
    
    idx2label : dict
        Map from zero-based indices to chromosome labels.
        
    map : dict
        Mapping from labels to interFiles. Only available when
        ``Format = 'TXT'``.
    
    See Also
    --------
    numpy.load : Load an array(s) from .npy or .npz files.
    numpy.savez : Save several arrays into a single .npz file.
    numpy.dtype : Create a data type object.
    
    Notes
    -----
    Our pipeline is more suitable for high-resolution data.
    (``resolution <= 50000``)
    
    References
    ----------
    .. [1] Lieberman-Aiden E, van Berkum NL, Williams L et al. Comprehensive
       mapping of long-range interactions reveals folding principles of the
       human genome. Science, 2009, 326: 289-293.
           
    .. [2] Dekker J, Rippe K, Dekker M et al. Capturing chromosome
       conformation. Science, 2002, 295: 1306-1311.
           
    .. [3] Imakaev M, Fudenberg G, McCord RP et al. Iterative correction of
       Hi-C data reveals hallmarks ofchromosome organization. Nat Methods,
       2012, 9: 999-1003.
    
    """
    
    def __init__(self, path = '.', Format = 'TXT', resolution = 10000,
                 template = 'chr%s_chr%s.int', chroms = ['#', 'X'],
                 cols = None, prefix = None, immortal = False, saveto = None):
        # Set main attributes
        self.location = os.path.abspath(os.path.expanduser(path))
        self.chroms = set(chroms)
        self.template = template
        self.resolution = resolution
        self.cols = cols
        
        log.debug('Initializing an Inters object ...')
        
        # Argument types and dependency relationships
        if (immortal == True) and (not saveto):
            log.error('Setting Error: saveto is required if immortal'
                      ' is True')
            sys.exit(1)
        
        if (Format == 'TXT') and (not self.template):
            log.error('Setting Error: template is required if Format '
                      'is set to be TXT')
            sys.exit(1)
        
        if (Format == 'NPZ') and (not prefix):
            log.error('Setting Error: prefix is required if Format '
                      'is set to be NPZ')
            sys.exit(1)

        if not Format in ['TXT', 'NPZ']:
            log.error('Setting Error: Unknown source format')
            sys.exit(1)
        
        # Summary of actions
        log.debug('Following tasks will be done:')
        if (Format == 'NPZ') and immortal:
            log.debug('Argument template and cols will be ignored')
            log.debug('Load interaction data from a Numpy .npz file')
            log.debug('Save interaction data to another Numpy .npz file')
        elif (Format == 'NPZ') and (not immortal):
            log.debug('Argument template and cols will be ignored')
            log.debug('Load interaction data from a Numpy .npz file')
        elif (Format == 'TXT') and (immortal):
            log.debug('Load interaction data from TXT files')
            log.debug('Save interaction data to a binary file in Numpy .npz'
                      ' format')
        else:
            log.debug('Load interaction data from TXT files')
        
        # Interaction data will be stored in a dict, with chromosome labels
        # as the keys
        self.data = {}
        
        if Format == 'NPZ':
            npzfile = '.'.join([prefix, 'npz'])
            log.debug('Data source: %s', npzfile)
            log.debug('Loading ...')
            if not os.path.exists(os.path.join(self.location, npzfile)):
                log.error('There is no file named %s at %s ', 
                          (npzfile, self.location))
                sys.exit(2)
            else:
                tempfile = np.load(os.path.join(self.location, npzfile))
                self.labels = []
                for i in tempfile.files:
                    # Filter using chroms argument
                    if ((not self.chroms) or (i.isdigit() and '#' in self.chroms)
                        or (i in self.chroms)):
                        self.data[i] = tempfile[i]
                        self.labels.append(i)
                # Sort labels, compatible with 'TXT' case
                self._sortlabels()
                self.labels = zip(*sorted(self.idx2label.items(),
                                          key = lambda x: x[0]))[1]
                
                log.debug('Done!')
                tempfile.close()
        else:
            log.debug('Scanning folder: %s' % self.location)
            # Scan the folder and obtain a list of chromosome labels
            self._scanFolder()
            log.debug('Loading interaction data from TXT files, this '
                      'process may be time/space consuming')
            self._readInters()
        
        if immortal:
            npzfile = '.'.join([saveto, 'npz'])
            tempfile = os.path.join(self.location, npzfile)
            log.debug('Saving data into %s', tempfile)
            np.savez(tempfile, **self.data)
            log.debug('Done!')
            
    
    def _extractChrLabel(self, filename):
        """
        Extract chromosome label from a file name according to the template.
        
        """
        # Full filename including path prefix
        _, interName = os.path.split(filename)
        regexp = self.template % (('(.*)',) * self.template.count('%s'))
        search_results = re.search(regexp, filename)
        label = search_results.group(1)
        
        # Remove leading zeroes.
        if label.isdigit():
            label = str(int(label))
    
        return label

    def _scanFolder(self):
        """
        Create a map from chromosome labels to file names under the folder.
        
        """
        
        if not os.path.isdir(self.location):
            log.error('%s is not a folder', self.location)
            sys.exit(2)
        
        # Interaction file names
        self.interFiles = [os.path.join(self.location, i)
                           for i in glob.glob(os.path.join(self.location,
                           self.template % (('*',) * self.template.count('%s'))))]
        
        if len(self.interFiles) == 0:
            log.error('No files found at %s', self.location)
            sys.exit(2)
        
        log.debug('%d interaction files are found in total', 
                  len(self.interFiles))
        
        # Read chromosome labels
        log.debug('Filter interaction files according to chroms ...')
        self.labels = []
        filtered = [] # Depend on user's selection
        for i in self.interFiles:
            label = self._extractChrLabel(i)
            if ((not self.chroms) or (label.isdigit() and '#' in self.chroms)
                or (label in self.chroms)):
                self.labels.append(label)
                filtered.append(i)
        
        self.interFiles = filtered
        if len(self.interFiles) == 0:
            log.error('No interaction files left after filtering!')
            sys.exit(2)
        
        log.debug('%d interaction files are left after filtering', 
                  len(self.interFiles))
        
        self._sortlabels()
        
        # Sort interaction file names and labels according to the
        # indices
        self.labels = zip(*sorted(self.idx2label.items(),
                                  key = lambda x: x[0]))[1]
        self.interFiles.sort(key = lambda path: 
                             self.label2idx[self._extractChrLabel(path)])
        
        # Map from labels to files
        self.map = dict(zip(self.labels, self.interFiles))
        log.debug('Following interaction files will be considered for'
                  ' further analysis:')
        for i in self.interFiles:
            log.debug(i)
    
    def _sortlabels(self):
        """Sort chromosome labels and construct map between chromosome
        labels and zero-based indices.
        
        Numerical labels are sorted naturally. For non-numerical labels,
        give the priority to XYM over the rest.
        
        """
        if ('#' in self.chroms) or (not self.chroms):
            # Convert labels to indices:
            # For numerical labels
            num_labels = [i for i in self.labels if i.isdigit()]
            # Sort labels naturally, i.e. place '2' before '10'.
            num_labels.sort(key=lambda x: int(re.findall(r'\d+$', x)[0]))
            count = len(num_labels)
            self.label2idx = dict(
                [(num_labels[i], int(i)) for i in xrange(len(num_labels))])
            self.idx2label = dict(
                [(int(i), num_labels[i]) for i in xrange(len(num_labels))])
            # For non-numerical labels. Give the priority to XYM over
            # the rest
            nonnum_labels = [i for i in self.labels if not i.isdigit()]
            for i in ['M', 'Y', 'X']:
                if i in nonnum_labels:
                    nonnum_labels.pop(nonnum_labels.index(i))
                    nonnum_labels.insert(0, i)
            for i in nonnum_labels:
                self.label2idx[i] = count
                self.idx2label[count] = i
                count += 1
        else:
            self.label2idx = dict(
                [(self.labels[i], int(i)) for i in xrange(len(self.labels))])
            self.idx2label = dict(
                [(int(i), self.labels[i]) for i in xrange(len(self.labels))])
        
    def _readInters(self):
        """
        Read interaction data chromosome by chromosome.
        
        """
        
        itype = np.dtype({'names':['bin1', 'bin2', 'IF'],
                          'formats':[np.int, np.int, np.float]})
        
        # Try to parse the interaction file format
        # We assume the columns are separated by comma or any space
        log.debug('Parsing interaction files ...')
        if not self.cols:
            log.warning('You haven\'t specify the desired columns! An '
                        'automatic recognization will be performed, but '
                        'unexpected consequences may happen.')
        # Loop
        for i in self.labels:
            log.debug('Current file: %s', self.map[i])
            tempfile = open(self.map[i])
            skip = -1
            for line in tempfile:
                # Skip all comment lines
                skip += 1
                if not line.startswith('#'):
                    delimiter = None
                    parse = line.strip().split()
                    if len(parse) < 3:
                        delimiter = ','
                        parse = line.strip().split(',')
                    break
            if line.isspace():
                log.error('Empty! Please check your file!')
                sys.exit(2)
            
            if len(parse) < 3:
                log.error('This format "%s" is not allowed!', line.strip())
                sys.exit(2)
            
            cols = [idx for idx in range(len(parse)) if isnumber(parse[idx])]
            if len(cols) < 3:
                # Header? Try next line!
                line = tempfile.readline()
                skip += 1
                parse = line.strip().split(delimiter)
                cols = [idx for idx in range(len(parse)) if isnumber(parse[idx])]
                if len(cols) < 3:
                    # Format Error
                    log.error('Illegal Format: "%s"! Three columns '
                              'indicating two bins and corresponding'
                              ' interaction frequency are required!',
                              line.strip())
                    sys.exit(2)
            
            # When label is contained in the main text
            if len(cols) > 3:
                for line in tempfile:
                    parse = line.strip().split(delimiter)
                    cols = [idx for idx in range(len(parse)) if 
                            ((isnumber(parse[idx])) and (parse[idx] != i))]
                    if len(cols) == 3:
                        break
            if len(cols) != 3:
                log.error('Format "%s" cannot be resolved, more than '
                          '3 numerical columns may exist! Computer '
                          'is puzzled!', line.strip())
                sys.exit(2)
            else:
                if self.cols:
                    if (set(self.cols) & set(cols)) != set(self.cols):
                        warnings.warn('Your cols setting may have a trouble!'
                                      ' Note that the index of Python is '
                                      '0-based.')
                        
                    log.debug('Column %s, %s, and %s will be  parsed as '
                              '"bin1", "bin2", and "Interaction frequency", '
                              'respectively.', *tuple(self.cols))
                    
                    self.data[i] = np.loadtxt(self.map[i], dtype = itype,
                                              delimiter = delimiter,
                                              skiprows = skip,
                                              usecols = self.cols)
                else:
                    log.debug('Column %s, %s, and %s will be  parsed as '
                              '"bin1", "bin2", and "Interaction frequency", '
                              'respectively.', *tuple(cols))
                    
                    self.data[i] = np.loadtxt(self.map[i], dtype = itype,
                                              delimiter = delimiter,
                                              skiprows = skip,
                                              usecols = cols)
        log.debug('Interaction files are read successfully!')

class TAD(object):
    """
    Load TAD data from a TXT file.
    
    TAD -- Topologically Associating Domain.
    
    TADs reflect megabase-scale substructures of chromatin in mammals [1]_,
    even in some lower animals(Drosophila embryos) [2]_. In short, TADs are
    identified naturally from interaction data, and chromatin interactions
    occur preferentially among the bins within the same domain.
    
    Parameters
    ----------
    source : str
        Complete source file name.
        Suppose a file named "h1hesc.domain" is located in
        "/home/xtwang/data/TAD", then ``source = '~/data/TAD/h1hesc.domain'``
        or ``source = '/home/xtwang/data/TAD/h1hesc.domain'``.
    
    chromname : None or str
        Template of chromosome names.
        Suppose ``chromname = 'chr'``, then the source data should be as
        follows::
        
            chr1    0   350000
            chr1    350000  800000
            chr1    800000  1450000
        
        Default: None
    
    cols : list of int (length: 3)
        Which columns to read, with 0 being the first. For example,
        ``cols = [0, 1, 2]`` will extract the 1st, 2nd and 3rd columns.
        (Default: [0, 1, 2])
    
    Attributes
    ----------
    data : Structured Array
        3 fields are defined by numpy.dtype object, and they are "chr",
        "start" and "end".
    
    See Also
    --------
    numpy.dtype : Create a data type object
    
    References
    ----------
    .. [1] Dixon JR, Selvaraj S, Yue F et al. Topological domains in
       mammalian genomes identified by analysis of chromatin interactions.
       Nature, 2012, 485: 376-380.
       
    .. [2] Sexton T, Yaffe E, Kenigsberg E et al. Three-dimensional folding
       and functional organization principles of the Drosophila genome.
       Cell, 2012, 148: 458-472.
    
    """
    def __init__(self, source, chromname=None, cols=[0, 1, 2]):
        
        source = os.path.abspath(os.path.expanduser(source))
        
        log.debug('Initializing a TAD object ...')
        
        # Check settings
        if len(cols) != 3:
            log.error('Setting Error: the length of cols must be 3')
            sys.exit(1)
        
        if not os.path.exists(source):
            log.error('We cannot find %s' % source)
            sys.exit(1)
        
        log.debug('Loading TADs ...')
        # Read TAD data from source
        # Skip any comment lines
        pool = [line.strip().split() for line in open(source) if
                not line.startswith('#')]
        
        if not pool:
            log.error('Empty source file!')
            sys.exit(1)
        
        # Skip header lines
        pool = [line for line in pool if isnumber(line[cols[1]])]
        
        # Manipulate chromosome labels
        maxL = 0 # Maximum label length        
        if not chromname:
            chromname = ''
        for i in xrange(len(pool)):
            pool[i][cols[0]] = pool[i][cols[0]][len(chromname):]
            if len(pool[i][cols[0]]) > maxL:
                maxL = len(pool[i][cols[0]])
            # Standardize API
            pool[i] = (pool[i][cols[0]], int(pool[i][cols[1]]),
                       int(pool[i][cols[2]]))
        
        # Create a structured array
        dtype = np.dtype({'names':['chr', 'start', 'end'],
                          'formats':['S'+str(maxL), np.int, np.int]})
        
        self.data = np.array(pool, dtype = dtype)
        
        log.debug('External TADs are loaded successfully!')

class Core(object):
    """
    Interaction analysis at TAD level.
    
    The analysis uses an interaction matrix as the input, which contains
    all interaction information of a TAD. 
    
    High IFs off the diagonal region can be identified using
    :meth:`longrange`. :meth:`DBSCAN` performs a density-based clustering
    algorithm to detect aggregation patterns in those IFs. Furthermore,
    :meth:`gdensity` defines an interaction feature, which is proved
    closely correlating with active/inactive level.
    
    Parameters
    ----------
    matrix : numpy.ndarray, (ndim = 2)
        Interaction matrix of a TAD. Can be extracted by :func:`getmatrix`
        giving interaction data (under certain resolution). Each entry
        indicates interaction frequency between corresponding two bins.
    
    start : int
        Off-diagonal level, i.e., only entries that are at least **start**
        away from the diagonal will be taken into account.
    
    left : int
        Starting point of TAD, in **resolution** unit. For example, if bin
        size is 10kb, ``left = 50`` means position 500000(bp) on the genome.
    
    Attributes
    ----------
    matrix : numpy.ndarray, (ndim = 2)
        Interaction matrix.
    
    start : int
        Off-diagonal level.
    
    left : int
        Starting point of TAD.
    
    newM : numpy.ndarray, (ndim = 2)
        A modified **matrix**. All vacant rows and columns in original
        **matrix** are removed.
    
    convert : dict
        A coordinate map from **newM** to **matrix**.
    
    After calling :meth:`longrange`:
    
    dim : int
        Size of **newM**.
    
    diag : int
        The best "k" (the kth diagonal) for IF extraction.
    
    pos : numpy.ndarray, (shape = (N, 2))
        Coordinates of selected IFs in **newM**.
    
    sig : float
        Significance of selected IFs under poisson assumption.
    
    Np : int
        Number of selected IFs.
    
    aberrant : numpy.ndarray, (ndim = 1)
        Aberrant level, i.e., ratio of low IFs near the diagonal region,
        used to determine the best "k" value.
    
    x : numpy.ndarray, (ndim = 1)
        Indices of **aberrant** components.
    
    xs : numpy.ndarray, (ndim = 1)
        x coordinates of interpolating fitting curve of **aberrant**.
    
    ys : numpy.ndarray, (ndim = 1)
        y coordinates of the fitting curve.
    
    dy1 : numpy.ndarray, (ndim = 1)
        The first-order differences of **ys**.
    
    dy2 : numpy.ndarray, (ndim = 1)
        Second-order differences of **ys**.
    
    After calling :meth:`DBSCAN`:
    
    cluster_id : numpy.ndarray, (shape = (N,))
        Cluster labels for each point in **pos**. -1 indicates noisy points.
    
    Nc : int
        Cluster number.
    
    clusters : dict
        Details of each cluster. "Average density", "radius",
        "area (polygon)", and "object number" are all recorded.
    
    Hulls : list of :func:`tadlib.polygon.Polygon` object
        Record clusters which can be enclosed by a convex polygon.
    
    eps : float
        One of the DBSCAN algorithm inputs, calculated in an analytical
        way.
    
    MinPts : int
        Another parameter of DBSCAN algorithm.
    
    ptrace : list of int
        Labels of clusters in which all points are collinear. These
        clusters cannot be enclosed by a convex polygon.
    
    After :meth:`gdensity`:
    
    gden : float, [0, 1]
        Weighted average density. Calculated using cluster density
        information.
    
    label : {0, 1}
        1 if ``Nc > 0``.
    
    See Also
    --------
    getmatrix
    
    """
    
    def __init__(self, matrix, start = 2, left = 0):
        
        ## Set basic attributes
        # Initial Interaction matrix
        self.matrix = matrix
        # Off-diagonal level
        self.start = start
        # Starting point of the region
        self.left = left
        # Manipulation, remove vacant rows and columns
        self.newM, self.convert = manipulation(self.matrix, self.left)
    
    
    def longrange(self):
        """Extract high IFs away from the diagonal region of an interaction
        matrix.
        
        It's known that Hi-C data have revealed a widespread rule for
        chromatin folding -- "power-law". Fragments on chromatin always
        interact with genomic proximal ones and chromain are driven into
        higher order structure. Under the rule, all high IFs should locate
        in diagonal region for a TAD.
        
        However, there does exist high IFs which are far away from the
        diagonal region breaking the rule. These interactions may occur
        between promoters and enhancers, which have potential to regulate
        gene expression.
        
        Extracting those IFs may bring some new insight into the
        relationship between structures and regulations.
        
        Notes
        -----
        The interaction matrix is expressed as :math:`M_{ij}` with
        :math:`\\leq i \\leq n, i \\leq j \\leq n`. Here, :math:`M_{ij}`
        is an upper triangular matrix.
        
        For each diagonal :math:`k`, coordinates (i, j) are divided into
        two groups: :math:`R_{kc}`, where :math:`\\left|i-j\\right|\\leq k`,
        and :math:`R_{ko}`, where :math:`\\left|i - j\\right| > k`.
        
        IFs within :math:`R_{ko}` having equal strength with those in
        :math:`R_{kc}` are extracted. In this way, each :math:`k` has its
        own such IF population. However, different :math:`k` has different
        "aberrant level", which we use for determining the best :math:`k`:
        
        .. math:: \omega_k = \\frac{L_k}{N_k}
        
        where :math:`\omega_k` is so-called "aberrant level",
        :math:`L_k` is the number of IFs within :math:`R_{kc}` but having
        lower strength than selected IFs, and :math:`N_k` is the total
        number of :math:`R_{kc}`.
        
        :math:`\omega_k` always decreases sharply along with the first few
        :math:`k`. Therefore, we calculate first and second order differences
        of it and the first inflection point is regarded as the best
        :math:`k`.
        
        """
        
        self.dim = self.newM.shape[0]
        
        # Upper triangular matrix
        upper = np.triu(self.newM, k = self.start)
        nonzero = upper[np.nonzero(upper)]
        
        ## The ideological source -- Pareto distribution
        # A statistic table
        freq_base = itemfreq(nonzero)
        # Cumulative Population
        cumsum = np.cumsum(freq_base[:,1]) / nonzero.size
        # Cumulative wealth (IF)
        wealth = np.cumsum(freq_base[:,0] * freq_base[:,1]) / nonzero.sum()
        
        ## Preliminary screening
        ## For loop termination
        ## The maximum off-diagonal level
        # The richest 50%
        obj = np.ceil(0.5 * nonzero.size)
        # Coefficients of polynomial
        coeff = [-1, 2 * (self.dim - self.start) - 1,
                 2 * (self.dim - self.start) - 2 * obj]
        # Polynomial roots
        root = np.roots(coeff)
        # Constraint condition
        root = root[(root > 0) & (root < (self.dim - self.start))]
        n = int(np.ceil(root[0])) if len(root) > 0 else 0
        for i in range(n, self.dim - self.start):
            temp = upper.copy()
            temp[np.triu_indices(self.dim, i + self.start)] = 0
            residual = np.nonzero(temp)[0].size - obj
            if residual>=0:
                break
        # Calibrate the root
        n = i
        
        ## Assumed Poisson Model
        ## Trace the strength of selected IFs
        G = self.newM[np.triu_indices(self.dim, k = self.start)]
        Poiss = poisson(G.mean())
        
        ## Aberrant Power-law for IF selection
        # Selected IF positions and corresponding level
        candidate = []
        # Work in concert with Poiss
        thres = []
        # Aberrant level, i.e., ratio of low IFs near the diagonal, used to
        # determine the best off-diagonal level
        aberrant = []
        
        # Loop start
        for i in range(self.start, self.start + n + 1):
            top_p = upper.copy()
            ## Two parts
            # Off diagonal
            off = np.triu(self.newM, k = i + 1)
            out = off[np.nonzero(off)]
            # Close to diagonal
            top_p[np.triu_indices(self.dim, i + 1)] = 0
            nearest = top_p[np.nonzero(top_p)]
            
            # Find the percent on cumulative population curve which
            # reflects the wealth (IF) ranking
            arg = np.abs(cumsum - (1 - nearest.size / nonzero.size)).argmin()
            # The top x% wealth
            top = 1 - wealth[arg]
            total = nonzero.sum() * top
            # Cash flow away from the diagonal
            free =  total - top_p.sum()
            
            ## Compute threshold
            # Statistics for the out part
            out_base = itemfreq(out)
            # Arrange from large to small
            out_base = np.flipud(out_base)
            got = np.abs(np.cumsum(out_base[:,0] * out_base[:,1]) - free)
            thre = out_base[:,0][got.argmin()]
            
            ## Extract long-range interaction positions
            mask = off > thre
            temp = np.where(mask)
            extract = np.array([(temp[0][j], temp[1][j]) for j in 
                                range(temp[0].size)])
            
            if len(extract) > 0:
                candidate.append((extract, i))
                aberrant.append((nearest <= thre).sum() / nearest.size)
                thres.append(thre)
        
        ## Different cases
        if candidate:
            self.aberrant = np.array(aberrant)
            idx = self.aberrant.argmin()
            self.x = np.arange(self.aberrant.size)
            # Constraint for gradient calculation
            if self.x.size >= 3:
                # Constraint of B-spline (order = 3)
                if self.x.size > 4:
                    ## Linear Spline
                    ip = interp1d(self.x, self.aberrant)
                    # Downsample the data evenly
                    times = np.arange(2, 4)
                    scheck = self.x.size / times
                    snum = scheck[scheck > 6][-1] if (scheck > 6).sum() > 0 else self.x.size
                    xi = np.linspace(0, self.x.size - 1, snum)
                    yi = ip(xi)
                    
                    ## B-spline
                    tcl = splrep(xi, yi)
                    self.xs = np.linspace(0, self.x.size - 1,
                                          self.x.size * 5)
                    self.ys = splev(self.xs, tcl)
                    # Finite differences
                    self.dy1 = np.gradient(self.ys)
                    self.dy2 = np.gradient(self.dy1)
                else:
                    self.xs = self.x
                    self.ys = self.aberrant
                    self.dy1 = np.gradient(self.aberrant)
                    self.dy2 = np.gradient(self.dy1)
                
                ## Inflection point
                m = (self.dy2[1:] <= 0) & (self.dy2[:-1] >= 0)
                for i in np.where(m)[0]:
                    tid = np.int(np.ceil(self.xs[i]))
                    if self.aberrant[tid] < self.aberrant[0]:
                        idx = tid
                        break
            
            # Selected interaction positions
            self.pos = candidate[idx][0]
            # Off-diagonal level
            self.diag = candidate[idx][1]
            # Significance of selected IFs
            self.sig = 1 - Poiss.cdf(thres[idx])
        else:
            ## Dummy assignment
            self.diag = 0; self.pos = []
            self.sig = np.nan
        
        # The number
        self.Np = len(self.pos)
        self.candidate = candidate
    
    def DBSCAN(self):
        """Detect natural patterns in selected IFs using DBSCAN.
        
        DBSCAN is a dennsity-based clustering algorithm. [1]_ Two input
        parameters **eps** and **MinPts** are calculated in an analytical
        way. [2]_ Before further analysis, some basic features are
        calculated for each cluster, including "density", "radius" and the
        "area", among which, "area" stands for corresponding convex polygon
        area.
        
        See Also
        --------
        sklearn.cluster.DBSCAN : an implementation of DBSCAN
        epsilon : calculate **eps** and **MinPts**
        tadlib.polygon.Polygon : calculations based on polygon.
        
        Notes
        -----
        Both "radius" and "density" are defined based on core objects of a
        cluster. "radius" is the average distance from the core object to its
        MinPts-nearest neighbors while "density" is the inverse of it.
        
        References
        ----------
        .. [1] Ester M, Kriegel H, Sander J, Xu X. A density-based
           algorithm for discovering clusters in large spatial databases
           with noise. Proc. 2nd Int. Conf. on Knowledge Discovery and Data
           Mining, Portland, OR, AAAI, Press, 1996, pp. 226-231
        
        .. [2] Daszykowski M, Walczak B, Massart DL. Looking for Natural
           Patterns in Data, Part 1: Density Based Approach. Chemom. Intell.
           Lab. Syst., 2001, 56: 83-92.
               
        """
        
        ## Generalize
        # Pseudo cluster number
        self.Nc = 0
        # Trace for scatter plot
        self.ptrace = []
        
        # Lower bound for input
        if self.Np >= 3:
            self.epsilon()
            # A simple but prerequisite condition
            if self.eps > 0:
                # Density-based cluster identification
                db = cluster.DBSCAN(eps = self.eps,
                                    min_samples = self.MinPts).fit(self.pos)
                # Cluster Label, -1 means noise
                self.cluster_id = db.labels_.astype(int)
                
                # Correction, there may be a bug in DBSCAN
                table = itemfreq(self.cluster_id)
                mask = (table[:,0] != -1) & (table[:,1] == 1)
                ridx = table[mask][:,0]
                mask = np.ones(self.Np, dtype = bool)
                for i in ridx:
                    mask &= (self.cluster_id != i)
                self.pos = self.pos[mask]
                self.cluster_id = self.cluster_id[mask]
                
                # Number of clusters
                self.Nc = len(set(self.cluster_id)) - \
                              (1 if -1 in self.cluster_id else 0)
                
                if self.Nc > 0:
                    ## Cluster-based attributes
                    ## Calculate average density / radius of each cluster
                    ## Only core objects are taken into account
                    ## Area, objects and object number are also recorded
                    self.clusters = {}
                    self.clusters['density'] = np.zeros(self.Nc)
                    self.clusters['radius'] = np.zeros(self.Nc)
                    self.clusters['areas'] = np.zeros(self.Nc)
                    self.clusters['obj'] = []
                    self.clusters['Cn'] = np.zeros(self.Nc, dtype = int)
                    
                    # Convex Hulls
                    self.Hulls = []
                    
                    ## Core points may also change after correction
                    core_indices = db.core_sample_indices_
                    cores_mask = np.zeros(self.Np, dtype = bool)
                    cores_mask[core_indices] = True
                    cores_mask = cores_mask[mask]
                    
                    # For each cluster
                    for i in xrange(self.Nc):
                        extract = (self.cluster_id == i)
                        t_C = self.pos[extract]
                        # Objects
                        self.clusters['obj'].append(t_C)
                        # Total object number
                        self.clusters['Cn'][i] = extract.sum()
                        
                        ## Average radius / density
                        # Core points of current cluster
                        choose = extract & cores_mask
                        cores = self.pos[choose]
                        # Distances from core points to any other points in
                        # current cluster
                        dists = distance.cdist(cores, t_C)
                        dists.sort()
                        # Average radius
                        self.clusters['radius'][i] = np.mean(dists[:, 1:self.MinPts].sum(axis = -1) / \
                                                             (self.MinPts - 1))
                        # Inverse of average radius, i.e., density
                        self.clusters['density'][i] = np.mean((self.MinPts - 1) / \
                                    dists[:, 1:self.MinPts].sum(axis = -1))
                                    
                        # Collinear test
                        judge = polygon.collinear(t_C)
                        if not judge:
                            # Represented as a convex polygon
                            P = polygon.Polygon(t_C)
                            # Area of the polygon
                            P.calarea()
                            self.clusters['areas'][i] = P.area
                            self.Hulls.append(P)
                            
                        else:
                            self.ptrace.append(i)
                            self.clusters['areas'][i] = 0
                
                # Change the number of points
                self.Np = len(self.pos)
    
    def gdensity(self):
        """Weighted density calculation.
        
        :meth:`longrange` and :meth:`DBSCAN` have to be called in advance.
        
        Density of a TAD is the weighted average density of each cluster.
        Weight is the ratio of object number of a cluster to :attr:`Np`.
        
        """
        if self.Nc > 0:
            W = self.clusters['Cn'] / self.Np
            self.gden = np.sum(self.clusters['density'] * W)
            self.label = 1
        else:
            self.gden = 0
            self.label = 0
        
            
    def epsilon(self):
        """Analytical way of estimating input parameters for DBSCAN.
        
        """
        ## Neighborhood of a point
        if len(self.pos.shape) > 1:
            # m objects, n variables
            m, n = self.pos.shape
        else:
            m = self.pos.shape[0]
            n = 1
        
        # Minimum number of points considered as a cluster
        self.MinPts = np.int(np.ceil(m / 25)) + 1
        # Enclose the total space
        prod = np.prod(self.pos.max(axis = 0) - self.pos.min(axis = 0))
        # Interpolation
        gamma = special.gamma(0.5 * n + 1)
        denom = (m * np.sqrt(np.pi ** n))
        self.eps = ((prod * self.MinPts * gamma) / denom) ** (1.0 / n)
    
    def plot(self, prefix = 'image', F = 'png', dpi = 300, s = 35,
             original = False, cmap = my_cmap, ec = '#4C4C4C'):
        """Graphical representation of every analysis stage.
        
        3 figures are created: Heat map of interaction matrix, scatter plot
        of selected IFs, scatter plot of clusters (enclosed by a convex
        polygon).
        
        Parameters
        ----------
        prefix : str
            Prefix of the output figure names. Default: 'image'
        
        F : str
            Format of the figures. Default: 'png'
        
        dpi : int
            The resolution in dots per inch. Directly delivered to
            plt.savefig. Default: 300
        
        s : int
            Size of points in scatter plot. Default: 35
        
        original : bool
            True if use the original interaction matrix. Default: False
        
        cmap : matplotlib.colors.Colormap
            A matplotlib.colors.Colormap instance. Default: :data:`my_cmap`
        
        ec : str
            Line color. Default: '#4C4C4C' (Dark grey)
        
        """
        ## Heatmap
        # Use the original matrix
        if original:
            matrix = self.matrix
        # Revised matrix
        else:
            matrix = self.newM
        
        nonzero = matrix[np.nonzero(matrix)]
        # 95-th percentile, to clarify the heatmap
        p = np.percentile(nonzero, 95)
        plt.imshow(matrix, origin = 'lower', interpolation = 'none',
                   cmap = cmap, vmax = p)
        plt.colorbar(format = '%d')
        plt.savefig('.'.join([prefix, 'heatmap', F]), dpi = dpi)
        plt.close()
        
        if self.Np > 0:
            if original:
                # Coordinates-conversion
                x = np.r_[[self.convert[i[0]] for i in self.pos]]
                y = np.r_[[self.convert[i[1]] for i in self.pos]]
                orid = self.convert[self.diag]
                bx = np.arange(len(matrix) - orid + 1)
                by = bx + orid
            else:
                x = self.pos[:,0]
                y = self.pos[:,1]
                bx = np.arange(self.dim - self.diag + 1)
                by = bx + self.diag
        
        ## Selected interactions
        if self.Np > 0:
            fig = plt.figure(figsize = (8, 8))
            plt.scatter(x, y, c = cmap(0.9), edgecolor = 'none', s = s)
            # Boundary (k-th diagonal)
            plt.plot(bx, by, color = ec, linestyle = '--', linewidth = 0.8)
            # Axis Control
            plt.xlim(0, len(matrix))
            plt.ylim(0, len(matrix))
            plt.axis('equal')
            plt.title('k = ' + str(self.diag))
            plt.savefig('.'.join([prefix, 'selected', F]), dpi = dpi)
            plt.close(fig)
        
        ## Clusters
        if self.Np > 0:
            # Color Map
            spectral = plt.cm.rainbow
            fig = plt.figure(figsize = (8, 8))
            if self.Nc > 0:
                noise = self.cluster_id == -1
                plt.scatter(x[noise], y[noise], c = 'k', edgecolor = 'none',
                            s = s)
                # In cluster but not in hulls
                for i in self.ptrace:
                    mask = self.cluster_id == i
                    # Random colors
                    seed = spectral(np.random.rand())
                    cs = [seed for i in xrange(len(x))]
                    plt.scatter(x[mask], y[mask], edgecolor = 'none', s = s,
                                c = cs)
                    plt.plot(x[mask], y[mask], color = ec, linestyle = '-')
                    
                # In hulls
                for h in self.Hulls:
                    if original:
                        points = np.r_[[(self.convert[i[0]],
                                         self.convert[i[1]])
                                       for i in h.points]]
                    else:
                        points = h.points
                    seed = spectral(np.random.rand())
                    cs = [seed for i in xrange(len(points))]
                    plt.scatter(points[:,0], points[:,1], edgecolor = 'none',
                                s = s, c = cs)
                    # Plot the Polygon
                    for simplex in h.simplices:
                        plt.plot(points[simplex, 0], points[simplex, 1],
                                 color = ec, linestyle = '-')
            else:
                # All points are noise
                plt.scatter(x, y, c = 'k', edgecolor = 'none', s = s)
            # Boundary
            plt.plot(bx, by, color = ec, linestyle = '--', linewidth = 0.8)
            # Axis Control
            plt.xlim(0, len(matrix))
            plt.ylim(0, len(matrix))
            plt.axis('equal')
            plt.title('k = ' + str(self.diag))
            plt.savefig('.'.join([prefix, 'clusters', F]), dpi = dpi)
            plt.close(fig)
    
    def kchoose(self, prefix = 'Aberrant', F = 'png', dpi = 300,
                ec = '#4C4C4C'):
        """
        Graphical representation of principles in :attr:`diag` determination.
        
        :attr:`aberrant` and its finite differences are plotted, and
        :attr:`diag` is labeled on the figures.
        
        Parameters
        ----------
        prefix : str
            Prefix of the output figure names. Default: 'Aberrant'
        
        F : str
            Format of the figures. Default: 'png'
        
        dpi : int
            The resolution in dots per inch. Directly delivered to
            plt.savefig. Default: 300
        
        ec : str
            Line color. Default: '#4C4C4C' (Dark grey)
        
        """
        if self.Np > 0:
            if self.x.size > 2:
                ## The trend of aberrant level
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.x, self.aberrant, 'ro', self.xs, self.ys)
                # Marker
                idx = self.diag - self.start
                ax.vlines(idx, ymin = self.aberrant.min(),
                          ymax = self.aberrant.max(),
                          linestyles = 'dashed', colors = ec)
                axis_control(ax)
                ax.set_xlabel('k')
                ax.set_ylabel('Percentage')
                ax.set_xticks(self.x)
                ax.set_xticklabels([i + self.start for i in self.x])
                plt.savefig('.'.join([prefix, 'trend', F]), dpi = dpi)
                plt.close()
            
                ## First-order difference
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.xs, self.dy1)
                ax.vlines(idx, ymin = self.dy1.min(), ymax = self.dy1.max(),
                          linestyles = 'dashed', colors = ec)
                axis_control(ax)
                ax.set_xlabel('k')
                ax.set_xticks(self.x)
                ax.set_xticklabels([i + self.start for i in self.x])
                plt.savefig('.'.join([prefix, 'diff-1', F]), dpi = dpi)
                plt.close()
            
                ## Second-order difference
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.xs, self.dy2)
                ax.vlines(idx, ymin = self.dy2.min(), ymax = self.dy2.max(),
                          linestyles = 'dashed', colors = ec)
                axis_control(ax)
                ax.set_xlabel('k')
                ax.set_xticks(self.x)
                ax.set_xticklabels([i + self.start for i in self.x])
                plt.savefig('.'.join([prefix, 'diff-2', F]), dpi = dpi)
                plt.close()
    
def extract_kth_diag(M, k = 0):
    """Coordinates of specified diagonal.
    
    Parameters
    ----------
    M : numpy.ndarray, (shape = (N, N))
        Array from which the diagonal is taken.
    
    k : int
        Offset of the diagonal from the main diagonal. Positive or negative.
        Default: 0
    
    Returns
    -------
    pos : numpy.ndarray, (shape = (x, 2))
        Coordinates of the specified diagonal.
    
    See Also
    --------
    numpy.diagonal : Return sepcified diagonals
    
    Examples
    --------
    >>> import numpy as np
    >>> from tadlib.analyze import extract_kth_diag
    >>> M = np.random.rand(4, 4) # A Random Matrix
    >>> pos = extract_kth_diag(M, k = 1)
    >>> print pos
    [[0 1]
     [1 2]
     [2 3]]
    >>> np.diagonal(M, offset = 1)
    array([ 0.78321734,  0.63147834,  0.64473559])
    >>> M[pos[:,0], pos[:,1]]
    array([ 0.78321734,  0.63147834,  0.64473559])
        
    """
    ridx, cidx = np.diag_indices_from(M)
    # ridx and cidx share the same buffer
    cidx = cidx.copy()
    
    if k > 0:
        cidx += k
    else:
        cidx -= k
    
    k = np.abs(k)
    
    x = ridx[:-k]; y = cidx[:-k]
    pos = np.array([(x[i], y[i]) for i in xrange(x.size)])
    
    return pos

def axis_control(ax):
    """A customized axis control operations.
    
    Parameters
    ----------
    ax : AxesSubplot
        An AxesSubplot instance.
    
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y',labelsize=9,direction='out')
    ax.tick_params(axis='x',labelsize=9,direction='out')
        

def isnumber(value):
    """A numerical string or not.
    
    A string is numerical if it can be converted using **float** function.
    
    Parameters
    ----------
    value : str
    
    Returns
    -------
    True if value is numerical.
    
    """
    try:
        float(value)
    except:
        return False
    else:
        return True
                    
def getmatrix(inter, l_bin, r_bin):
    """Extract regional interaction data and place it into a matrix.
    
    Parameters
    ----------
    inter : numpy structured array
        Three fields are required, "bin1", "bin2" and "IF", data types of
        which are int, int and float, respectively.
    
    l_bin : int
        Left bin index of the region.
        
    r_bin : int
        Right bin index of the region.
        
    Returns
    -------
    inter_matrix : numpy.ndarray
        The value of each entry is the interaction frequency between
        corresponding two bins.
    
    Notes
    -----
    Original interaction data is always binned under some resolution (10 kb
    in our research). **inter** is used to store such binned data, which can
    be seen as a sparse matrix in a numpy-style way.
    
    Sparse matrix is always required for high-resolution analysis of large
    genome, such as human and mouse.
        
    """
    # Construct a matrix
    inter_matrix = np.zeros((r_bin - l_bin, r_bin - l_bin), dtype = float)
    # Extract the regional data
    mask = (inter['bin1'] >= l_bin) & (inter['bin1'] < r_bin) & \
           (inter['bin2'] >= l_bin) & (inter['bin2'] < r_bin)
    inter_extract = inter[mask]
    
    # Fill the matrix
    for i in inter_extract:
        # Off-diagonal parts
        if i['bin1'] != i['bin2']:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
            inter_matrix[i['bin2'] - l_bin][i['bin1'] - l_bin] += i['IF']
        else:
            # Diagonal part
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
    
    return inter_matrix

def manipulation(matrix, start = 0):
    """Remove vacant rows and columns of a matrix.
    
    Parameters
    ----------
    matrix : numpy.ndarray, (ndim = 2)
        Interaction matrix generated by **getmatrix**.
    
    start : int
        The begining of the region. (Default: 0)
    
    Returns
    -------
    newM : numpy.ndarray, (ndim = 2)
        A modified **matrix** in which all vacant rows and columns are
        removed.
    
    convert : dict
        Index map from **newM** to **matrix**.
    
    See Also
    --------
    getmatrix
    
    Examples
    --------
    >>> import numpy as np
    >>> from tadlib.analyze import manipulation
    >>> matrix = np.random.rand(4, 4)
    >>> matrix[1,:] = 0; matrix[:,1] = 0
    >>> print matrix
    [[ 0.24822414  0.          0.07782508  0.01812965]
     [ 0.          0.          0.          0.        ]
     [ 0.93870151  0.          0.21986474  0.20462965]
     [ 0.13022712  0.          0.78674168  0.77068304]]
    
    >>> newM, convert = manipulation(matrix)
    >>> print newM
    [[ 0.24822414  0.07782508  0.01812965]
     [ 0.93870151  0.21986474  0.20462965]
     [ 0.13022712  0.78674168  0.77068304]]
    
    >>> print convert
    {0: 0, 1: 2, 2: 3}

    """
    mask = matrix.sum(axis = 0) == 0
    index = list(np.where(mask)[0])
    # Remove vacant rows
    temp = np.delete(matrix, index, 0)
    # Vacant columns
    newM = np.delete(temp, index, 1)
    convert = dict(zip(np.arange(len(newM)),
                       np.where(np.logical_not(mask))[0] + start))
    
    return newM, convert
