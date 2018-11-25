# Created on Sat Sep 27 16:18:27 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

from __future__ import division
import glob, re, os, sys
import logging
from tadlib.calfea import polygon
import numpy as np
from scipy.stats import poisson, itemfreq
from scipy.spatial import distance
from scipy.interpolate import interp1d, splrep, splev
import scipy.special as special
from sklearn import cluster

## Customize the logger
log = logging.getLogger(__name__)

def load_TAD(source_fil, chromname=None, cols=[0, 1, 2]):
    """
    Load TAD data from a TXT file.

    Parameters
    ----------
    source : str
        Path to the TAD file.
    
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
    
    Returns
    -------
    data :  Numpy Structured Array
        The parsed TAD intervals are contained in a numpy structured array
        containing 3 fields: "chr", "start" and "end".

    """
    source = os.path.abspath(os.path.expanduser(source_fil))

    # Read TAD data from source
    # Skip any comment lines
    pool = [line.strip().split() for line in open(source, 'r') if
            not line.startswith('#')]
    # Skip header lines
    pool = [line for line in pool if ((isnumber(line[cols[1]])) and (isnumber(line[cols[2]])))]

    # Manipulate chromosome labels
    maxL = 0 # Maximum label length
    if not chromname:
        chromname = ''
    for i in range(len(pool)):
        pool[i][cols[0]] = pool[i][cols[0]][len(chromname):]
        if len(pool[i][cols[0]]) > maxL:
            maxL = len(pool[i][cols[0]])
        # Standardize API
        pool[i] = (pool[i][cols[0]], float(pool[i][cols[1]]),
                float(pool[i][cols[2]]))
        
    # Create a structured array
    dtype = np.dtype({'names':['chr', 'start', 'end'],
                      'formats':['U'+str(maxL), np.int, np.int]})  
    
    data = np.array(pool, dtype = dtype)

    return data

class Core(object):
    """
    Interaction analysis at TAD level.
    
    High IFs off the diagonal region can be identified using
    :py:meth:`tadlib.calfea.analyze.Core.longrange`. :py:meth:`tadlib.calfea.analyze.Core.DBSCAN`
    performs a density-based clustering algorithm to detect aggregation patterns
    in those IFs. Furthermore, two structural features, called AP
    (Aggregation Preference) and Coverage in our original research, can be
    calculated by :py:meth:`tadlib.calfea.analyze.Core.gdensity` and
    :py:meth:`tadlib.calfea.analyze.Core.totalCover` respectively.
    
    Parameters
    ----------
    matrix : numpy.ndarray, (ndim = 2)
        Interaction matrix of a TAD.
    
    left : int
        Starting point of TAD. For example, if the bin size is 10kb,
        ``left = 50`` means position 500000(bp) on the genome.
    
    Attributes
    ----------
    newM : numpy.ndarray, (ndim = 2)
        Gap-free interaction matrix.
    
    convert : list
        Information required for converting *newM* to *matrix*.
    
    cEM : numpy.ndarray, (ndim = 2)
        Expected interaction matrix. An upper triangular matrix. Value in each
        entry will be used to construct a Poisson Model for statistical
        significance calculation.
    
    fE : numpy.ndarray, (ndim = 2)
        An upper triangular matrix. Each entry represents the fold enrichment
        of corresponding observed interaction frequency.
    
    Ps : numpy.ndarray, (ndim = 2)
        An upper triangular matrix. Value in each entry indicates the p-value
        under corresponding Poisson Model.
    
    pos : numpy.ndarray, (shape = (N, 2))
        Coordinates of the selected IFs in *newM*.
    
    Np : int
        Number of the selected IFs.
    
    cluster_id : numpy.ndarray, (shape = (N,))
        Cluster labels for each point in *pos*. -1 indicates noisy points.
    
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
    
    gden : float, [0, 1]
        Weighted average density.
    
    coverage : float, [0, 1]
        Total coverage of clusters.
    
    """
    def __init__(self, matrix, left = 0):

        # rescale matrix
        nonzero = matrix[matrix.nonzero()]
        if np.median(nonzero) < 1:
            min_nonzero = nonzero.min()
            scale = 1 / min_nonzero
            matrix = matrix * scale
        
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
            Width of the peak region. Default: 2
        
        ww : int
            Width of the donut. Default: 5
        
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
        
        mBM = np.zeros_like(BM)
        Mask = CM != 0
        mBM[Mask] = BM[Mask] / CM[Mask]
        
        ## Fold Enrichments
        self.fE = np.zeros_like(self.cEM)
        mask = self.cEM != 0
        self.fE[mask] = upper[mask] / self.cEM[mask]
        
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
        """Detect natural patterns in selected interactions using DBSCAN.
        
        DBSCAN is a dennsity-based clustering algorithm. [1]_ Two input
        parameters *eps* and *MinPts* are calculated in an analytical
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
                    for i in range(self.Nc):
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
        
        :py:meth:`tadlib.calfea.analyze.Core.longrange` and
        :py:meth:`tadlib.calfea.analyze.Core.DBSCAN` have to be called in advance.
        
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
        
        :py:meth:`tadlib.calfea.analyze.Core.longrange` and
        :py:meth:`tadlib.calfea.analyze.Core.DBSCAN` have to be called in advance.
        """
        if self.Nc > 0:
            csum = self.clusters['areas'].sum()
            self.coverage = csum / self._area
        else:
            self.coverage = 0
    
    def convertMatrix(self, M):
        """
        Convert an internal gap-free matrix(e.g., newM, cEM, fE, and Ps)
        into a new matrix with the same shape as the original interaction
        matrix by using the recorded index map(see the *convert* attribute).
        """
        idx = sorted(self.convert[0].values())
        newM = np.zeros((self.convert[1], self.convert[1]), dtype=M.dtype)
        y,x = np.meshgrid(idx, idx)
        newM[x,y] = M
            
        return newM
            
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
    """Remove gaps of the original interaction matrix.
    
    Parameters
    ----------
    matrix : numpy.ndarray, (ndim = 2)
        Interaction matrix.
    
    start : int
        The begining of the region. (Default: 0)
    
    Returns
    -------
    newM : numpy.ndarray, (ndim = 2)
        The gap-removed matrix.
    
    convert : list
        The first element is the index map from *newM* to *matrix*, and
        the second element records the length of *matrix*.
    
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
    [{0: 0, 1: 2, 2: 3}, 4]

    """
    mask = matrix.sum(axis = 0) == 0
    index = list(np.where(mask)[0])
    # Remove vacant rows
    temp = np.delete(matrix, index, 0)
    # Vacant columns
    newM = np.delete(temp, index, 1)
    mapidx = dict(zip(np.arange(len(newM)),
                      np.where(np.logical_not(mask))[0] + start))
    convert = [mapidx, matrix.shape[0]]
    
    return newM, convert
