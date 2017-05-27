# Created on Sat Sep 27 16:18:27 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

from __future__ import division
import glob, re, os, sys
import logging
import warnings
import polygon
import numpy as np
from scipy.stats import poisson, itemfreq
from scipy.spatial import distance
from scipy.interpolate import interp1d, splrep, splev
import scipy.special as special
from sklearn import cluster

## Customize the logger
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
    (``resolution <= 20000``)
    
    References
    ----------
    .. [1] Lieberman-Aiden E, van Berkum NL, Williams L et al. Comprehensive
       mapping of long-range interactions reveals folding principles of the
       human genome. Science, 2009, 326: 289-293.
           
    .. [2] Dekker J, Rippe K, Dekker M et al. Capturing chromosome
       conformation. Science, 2002, 295: 1306-1311.
           
    .. [3] Imakaev M, Fudenberg G, McCord RP et al. Iterative correction of
       Hi-C data reveals hallmarks of chromosome organization. Nat Methods,
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
            log.debug('    Argument template and cols will be ignored')
            log.debug('    Load interaction data from a Numpy .npz file')
            log.debug('    Save interaction data to another Numpy .npz file')
        elif (Format == 'NPZ') and (not immortal):
            log.debug('    Argument template and cols will be ignored')
            log.debug('    Load interaction data from a Numpy .npz file')
        elif (Format == 'TXT') and (immortal):
            log.debug('    Load interaction data from TXT files')
            log.debug('''    Save interaction data to a binary file in Numpy .npz format''')
        else:
            log.debug('    Load interaction data from TXT files')
        
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
        pool = [line for line in pool if ((isnumber(line[cols[1]])) and (isnumber(line[cols[2]])))]
        
        # Manipulate chromosome labels
        maxL = 0 # Maximum label length        
        if not chromname:
            chromname = ''
        for i in xrange(len(pool)):
            pool[i][cols[0]] = pool[i][cols[0]][len(chromname):]
            if len(pool[i][cols[0]]) > maxL:
                maxL = len(pool[i][cols[0]])
            # Standardize API
            pool[i] = (pool[i][cols[0]], float(pool[i][cols[1]]),
                       float(pool[i][cols[2]]))
        
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
    :meth:`tadlib.calfea.analyze.Core.longrange`. :meth:`tadlib.calfea.analyze.Core.DBSCAN`
    performs a density-based clustering algorithm to detect aggregation patterns
    in those IFs. Furthermore, two structural features, called AP
    (Aggregation Preference) and Coverage in our original research, can be
    calculated by :meth:`tadlib.calfea.analyze.Core.gdensity` and
    :meth:`tadlib.calfea.analyze.Core.totalCover` respectively.
    
    Parameters
    ----------
    matrix : numpy.ndarray, (ndim = 2)
        Interaction matrix of a TAD. Can be extracted by
        :func:`tadlib.calfea.analyze.getmatrix`
        giving interaction data (under certain resolution). Each entry
        indicates interaction frequency between corresponding two bins.
    
    left : int
        Starting point of TAD, in **resolution** unit. For example, if bin
        size is 10kb, ``left = 50`` means position 500000(bp) on the genome.
    
    Attributes
    ----------
    newM : numpy.ndarray, (ndim = 2)
        A modified **matrix**. All vacant rows and columns in original
        **matrix** are removed.
    
    convert : dict
        A coordinate map from **newM** to **matrix**.
    
    After :meth:`tadlib.calfea.analyze.Core.longrange`:
    
    cEM : numpy.ndarray, (ndim = 2)
        Expected interaction matrix. An upper triangular matrix. Value in each
        entry will be used to construct a Poisson Model for statistical
        significance calculation.
    
    Ps : numpy.ndarray, (ndim = 2)
        An upper triangular matrix. Value in each entry indicates the p-value
        under corresponding Poisson Model.
    
    pos : numpy.ndarray, (shape = (N, 2))
        Coordinates of selected IFs in **newM**.
    
    Np : int
        Number of selected IFs.
    
    After :meth:`tadlib.calfea.analyze.Core.DBSCAN`:
    
    cluster_id : numpy.ndarray, (shape = (N,))
        Cluster labels for each point in **pos**. -1 indicates noisy points.
    
    Nc : int
        Cluster number.
    
    clusters : dict
        Details of each cluster. "Average density", "radius",
        "area (polygon)", "point coordinates", and "object number" are all
        recorded.
    
    Hulls : dict
        Details of each convex hull (clusters which can be enclosed by a
        convex polygon).
    
    ptrace : list of int
        Labels of clusters in which all points are collinear. These
        clusters cannot be enclosed by a convex polygon.
    
    After :meth:`tadlib.calfea.analyze.Core.gdensity`:
    
    gden : float, [0, 1]
        Weighted average density. Calculated using cluster density
        information.
    
    After :meth:`tadlib.calfea.analyze.Core.totalCover`:
    
    coverage : float, [0, 1]
        Total coverage of clusters.
    
    See Also
    --------
    tadlib.calfea.analyze.getmatrix
    
    """
    def __init__(self, matrix, left = 0):
        
        # Manipulation, remove vacant rows and columns
        self.newM, self.convert = manipulation(matrix, left)
        
        ## Determine proper off-diagonal level
        Len = self.newM.shape[0]
        idiag = np.arange(0, Len)
        iIFs = []
        for i in idiag:
            temp = np.diagonal(self.newM, offset = i)
            iIFs.append(temp.mean())
        iIFs = np.array(iIFs)
        
        idx = np.where(iIFs > 0)[0][0]
        
        if idx > 0:
            log.warning('''There is no data in diagonal region.''')
        
        self._start = idx
        IFs = iIFs[idx:]
        diag = idiag[idx:]
        
        self._Ed = _fitting(diag, IFs)
        
    def longrange(self, pw = 2, ww = 5, top = 0.7, ratio = 0.05):
        """
        Select statistically significant interactions of the TAD. Both
        genomic distance and local interaction background are taken into
        account.
        
        Parameters
        ----------
        pw : int
            Width of the interaction region surrounding the peak. Default: 2
        
        ww : int
            The size of the donut sampled. Default: 5
        
        top : float, [0.5, 1]
            Parameter for noisy interaction filtering. Default: 0.7
        
        ratio : float, [0.01, 0.1]
            Specifies the sample size of significant interactions.
            Default: 0.05
        
        Notes
        -----
        *pw* and *ww* are sensitive to data resolution. It performs well
        when we set *pw* to 4 and *ww* to 7 at 5 kb, and (2, 5) at 10 kb. [1]_
        
        References
        ----------
        .. [1] Rao, S.S., Huntley, M.H., Durand, N.C. et al. A 3D map of the
           human genome at kilobase resolution reveals principles of chromatin
           looping. Cell, 2014, 159: 1665-1680.
        
        """
        dim = self.newM.shape[0]
        
        ps = 2 * pw + 1 # Peak Size
        ws = 2 * ww + 1 # Initial window size
        bs = 2 * pw + 1 # B -- Blurry
        
        start = ww if (ww > self._start) else self._start
        # Upper triangular matrix
        upper = np.triu(self.newM, k = start)
        bUpper = np.triu(self.newM, k = 0)
        
        # Expanded Matrix
        expM = np.zeros((dim + ww*2, dim + ww*2))
        expBM = np.zeros((dim + ww*2, dim + ww*2))
        expM[ww:-ww, ww:-ww] = upper
        expBM[ww:-ww, ww:-ww] = bUpper
        
        tm = np.all((expBM == 0), axis = 0)
        Mask = np.zeros((dim + ww*2, dim + ww*2), dtype = bool)
        Mask[:,tm] = True
        Mask[tm,:] = True
        expCM = np.ones_like(expM, dtype = int)
        expCM[Mask] = 0
        
        ## Expected matrix
        EM_idx = np.triu_indices(dim, k = start)
        EM_value = self._Ed[EM_idx[1] - EM_idx[0] - self._start]
        EM = np.zeros((dim, dim))
        EM[EM_idx] = EM_value
        ## Expanded Expected Matrix
        expEM = np.zeros((dim + ww*2, dim + ww*2))
        expEM[ww:-ww, ww:-ww] = EM
        
        ## Construct pool of matrices for speed
        # Window
        OPool_w = {}
        EPool_w = {}
        ss = range(ws)
        for i in ss:
            for j in ss:
                OPool_w[(i,j)] = expM[i:(dim+i), j:(dim+j)]
                EPool_w[(i,j)] = expEM[i:(dim+i), j:(dim+j)]
        # Peak
        OPool_p = {}
        EPool_p = {}
        ss = range(ww-pw, ps+ww-pw)
        for i in ss:
            for j in ss:
                OPool_p[(i,j)] = expM[i:(dim+i), j:(dim+j)]
                EPool_p[(i,j)] = expEM[i:(dim+i), j:(dim+j)]
        
        # For Blurry Matrix
        OPool_b = {}
        OPool_bc = {}
        ss = range(ww-pw, bs+ww-pw)
        for i in ss:
            for j in ss:
                OPool_b[(i,j)] = expBM[i:(dim+i), j:(dim+j)]
                OPool_bc[(i,j)] = expCM[i:(dim+i), j:(dim+j)]
        
        ## Background Strength  --> Background Ratio
        bS = np.zeros((dim, dim))
        bE = np.zeros((dim, dim))
        for w in OPool_w:
            if (w[0] != ww) and (w[1] != ww):
                bS += OPool_w[w]
                bE += EPool_w[w]
        for p in OPool_p:
            if (p[0] != ww) and (p[1] != ww):
                bS -= OPool_p[p]
                bE -= EPool_p[p]
                
        bE[bE==0] = 1
        bR = bS / bE
        
        ## Corrected Expected Matrix
        cEM = EM * bR
        self.cEM = cEM
        
        ## Contruct the Blurry Matrix
        BM = np.zeros((dim, dim))
        CM = np.zeros((dim, dim), dtype = int)
        
        for b in OPool_b:
            BM += OPool_b[b]
            CM += OPool_bc[b]
            
        mBM = BM / CM
        Mask = np.isnan(mBM)
        mBM[Mask] = 0
        
        ## Poisson Models
        Poisses = poisson(cEM)
        Ps = 1 - Poisses.cdf(upper)
        self.Ps = Ps
        rmBM = mBM[EM_idx] # Raveled array
        # Only consider the top x%
        top_idx = np.argsort(rmBM)[np.int(np.floor((1-top)*rmBM.size)):]
        # The most significant ones
        rPs = Ps[EM_idx][top_idx]
        Rnan = np.logical_not(np.isnan(rPs)) # Remove any invalid entry
        RrPs = rPs[Rnan]
        sig_idx = np.argsort(RrPs)[:np.int(np.ceil(ratio/top*RrPs.size))]
        if sig_idx.size > 0:
            self.pos = np.r_['1,2,0', EM_idx[0][top_idx][Rnan][sig_idx], EM_idx[1][top_idx][Rnan][sig_idx]]
        else:
            self.pos = np.array([])
            
        self.Np = len(self.pos)
        self._area = EM_idx[0].size
        
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
        tadlib.calfea.polygon.Polygon : calculations based on polygon.
        
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
        
        self.Nc = 0
        # Trace for scatter plot
        self.ptrace = []
        
        # Lower bound for input
        if self.Np >= 5:
            self._epsilon()
            # Minimum epsilon
            if self._eps < np.sqrt(2):
                self._eps = np.sqrt(2)
            # A simple but prerequisite condition
            if self._eps > 0:
                # Density-based cluster identification
                db = cluster.DBSCAN(eps = self._eps,
                                    min_samples = self._MinPts).fit(self.pos)
                # Cluster Label, -1 means noise
                self.cluster_id = db.labels_.astype(int)
                table = itemfreq(self.cluster_id)
                mask = (table[:,0] != -1) & (table[:,1] == 1)
                ridx = table[mask][:,0]
                mask = np.ones(self.Np, dtype = bool)
                if len(ridx) > 0:
                    ridx.sort()
                    for i in ridx[::-1]:
                        mask &= (self.cluster_id != i)
                    for i in ridx[::-1]:
                        self.cluster_id[self.cluster_id > i] -= 1
                        
                self.pos = self.pos[mask]
                self.cluster_id = self.cluster_id[mask]
                
                # Number of clusters
                self.Nc = len(set(self.cluster_id)) - \
                              (1 if -1 in self.cluster_id else 0)
                
                if self.Nc > 0:
                    ## Cluster-based attributes
                    self.clusters = {}
                    self.clusters['density'] = np.zeros(self.Nc)
                    self.clusters['radius'] = np.zeros(self.Nc)
                    self.clusters['areas'] = np.zeros(self.Nc)
                    self.clusters['obj'] = []
                    self.clusters['Cn'] = np.zeros(self.Nc, dtype = int)
                    
                    # Hull-based attributes
                    self.Hulls = {}
                    self.Hulls['polygons'] = []
                    self.Hulls['density'] = []
                    
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
                        self.clusters['radius'][i] = np.mean(dists[:, 1:self._MinPts].sum(axis = -1) / \
                                                             (self._MinPts - 1))
                        # Inverse of average radius, i.e., density
                        self.clusters['density'][i] = np.mean((self._MinPts - 1) / \
                                    dists[:, 1:self._MinPts].sum(axis = -1))
                                    
                        # Collinear test
                        judge = polygon.collinear(t_C)
                        if not judge:
                            # Represented as a convex polygon
                            P = polygon.Polygon(t_C)
                            # Area of the polygon
                            P.calarea()
                            self.clusters['areas'][i] = P.area
                            self.Hulls['polygons'].append(P)
                            self.Hulls['density'].append(self.clusters['density'][i])
                            
                        else:
                            self.ptrace.append(i)
                            self.clusters['areas'][i] = 0
                
                # Change the number of points
                self.Np = len(self.pos)
    
    def gdensity(self):
        """Weighted density calculation.
        
        :meth:`tadlib.calfea.analyze.Core.longrange` and
        :meth:`tadlib.calfea.analyze.Core.DBSCAN` have to be called in advance.
        
        Density of a TAD is the weighted average density of each cluster.
        Weight is the ratio of object number of a cluster to :attr:`Np`.
        
        """
        if self.Nc > 0:
            Num = len(self.pos[self.cluster_id != -1])
            W = self.clusters['Cn'] / Num
            self.gden = np.sum(self.clusters['density'] * W)
        else:
            self.gden = 0
    
    def totalCover(self):
        """
        Total coverage of clusters.
        
        :meth:`tadlib.calfea.analyze.Core.longrange` and
        :meth:`tadlib.calfea.analyze.Core.DBSCAN` have to be called in advance.
        """
        if self.Nc > 0:
            csum = self.clusters['areas'].sum()
            self.coverage = csum / self._area
        else:
            self.coverage = 0
        
            
    def _epsilon(self):
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
        self._MinPts = np.int(np.ceil(m / 25)) + 1
        # Enclose the total space
        prod = np.prod(self.pos.max(axis = 0) - self.pos.min(axis = 0))
        # Interpolation
        gamma = special.gamma(0.5 * n + 1)
        denom = (m * np.sqrt(np.pi ** n))
        self._eps = ((prod * self._MinPts * gamma) / denom) ** (1.0 / n)

def _fitting(x, y):
    
    ## Linear Spline
    ip = interp1d(x, y)
    # Downsample the data evenly
    times = np.arange(2, 4)
    scheck = x.size / times
    snum = scheck[scheck > 6][-1] if (scheck > 6).sum() > 0 else x.size
    xi = np.linspace(x.min(), x.max(), snum)
    yi = ip(xi)
    
    ## B-spline
    tcl = splrep(xi, yi)
    ys = splev(x, tcl)
    
    # Finite differences
    dy1 = np.gradient(ys)
    
    ## Unstable region
    m = (dy1[1:] >= 0) & (dy1[:-1] <= 0)
    if len(np.where(m)[0]) != 0:
        i = np.where(m)[0][0]
        ys[x > x[i]] = ys[i]
    
    return ys

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
        Interaction matrix generated by :func:`tadlib.calfea.analyze.getmatrix`.
    
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
    tadlib.calfea.analyze.getmatrix
    
    Examples
    --------
    >>> import numpy as np
    >>> from tadlib.calfea.analyze import manipulation
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
