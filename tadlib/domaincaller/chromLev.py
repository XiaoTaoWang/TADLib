# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:47:34 2016

@author: wxt
"""

from __future__ import division
import copy, collections
import numpy as np
from scipy import sparse
import matplotlib
matplotlib.use('agg')
from tadlib.hitad.aligner import BoundSet, DomainSet, DomainAligner, hierFormat, Container
from tadlib.calfea import analyze
from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        if self.vmin == self.vmax:
            return np.ma.masked_array(np.interp(value, [self.vmin], [0.5]))
        
        if self.vmin < self.midpoint < self.vmax:
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        elif self.vmin >= self.midpoint:
            x, y = [self.vmin, self.vmax], [0.5, 1]
        elif self.vmax <= self.midpoint:
            x, y = [self.vmin, self.vmax], [0, 0.5]
            
        return np.ma.masked_array(np.interp(value, x, y))

np.seterr(divide = "ignore")

class Chrom(object):
    """
    *Chrom* is defined to:
    
    - Hold Hi-C data within a certain chromosome
    - Identify hierarchical domains in 4 steps: 1.Calculate adaptive DIs.
      2.Identify original candidate bounds by 5-state Gaussian mixture Hidden
      Markov Model using adaptive DIs as input. 3.Select TAD bounds from
      candidate bounds. 4.Recursively identify inner domain structures of each TAD.
    - Visualize any region of the chromosome. Hierarchical domains will be
      plotted as boxes along with the diagonal of the heatmap, and adaptive
      DI track will be placed on top of the heatmap.
    
    Parameters
    ----------
    chrom : str
        Chromosome label.
    
    res : int
        Resolution of the Hi-C data in base-pair unit.
    
    hicdata : CSR sparse matrix
        Bin-level Hi-C matrix of the specified chromosome.
    
    Attributes
    ----------
    chrom : str
        Chromosome label.
    
    res : int
        Resolution in base-pair unit.
    
    chromLen : int
        Total bin number of the chromosome.
    
    rawMatrix : sparse matrix in Compressed Sparse Row format
        CSR sparse matrix is used to extract Hi-C data by slicing conveniently
        while guarantee low memory overhead.
    """
    defaultwindow = 2000000
    minsize = 5

    def __init__(self, chrom, res, hicdata):

        self.chrom = chrom
        self.res = res
        self._rm = 1
        self._dw = self.defaultwindow // res
        self.chromLen = hicdata.shape[0]
        self.hmm = None

        x, y = hicdata.nonzero()
        mat_ = hicdata[x, y]
        if isinstance(mat_, np.matrix):
            IF = np.array(mat_).ravel()
        else:
            # mat_ is a sparse matrix
            IF = np.array(mat_.todense()).ravel()

        IF[np.isnan(IF)] = 0
        self.rawMatrix = self._genSparseMatrix(x, y, IF)

        del x, y, IF, hicdata

        self._state = 'Submitted'

    def _genSparseMatrix(self, x, y, IF):

        extendLen = self.chromLen + 2*self.chromLen
        rawMatrix = sparse.csr_matrix((IF, (x + self.chromLen, y + self.chromLen)),
                                      shape = (extendLen, extendLen))

        return rawMatrix

    def detectPeaks(self, trends, mph=0, mpd=5):
        """
        Detect peaks (local maxima) in a 1-D array intuitively (a peak must
        be greater than its immediate neighbors).
        
        Parameters
        ----------
        trends : 1-D numpy ndarray
            Data.
        
        mph : float
            Only peaks that are greater than this value will be detected.
            (Default: 0)
        
        mpd : positive integer
            Only peaks whose indices are at least separated by this value will
            be reported. (Default: 5)
        
        Returns
        -------
        ind : 1-D numpy ndarray
            Indices of peaks detected in *trends*.
        """
        dx = trends[1:] - trends[:-1]
        ind = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]

        sp = np.where(trends==1)[0]
        if sp.size:
            ind = np.r_[sp[-1], ind]

        if dx[-1] > 0:
            ind = np.r_[ind, trends.size-1]

        # Filter peaks by mph
        if ind.size:
            ind = ind[trends[ind] > mph]

        # Remove small peaks closer than mpd
        if ind.size and mpd > 1:
            ind = ind[np.argsort(trends[ind])][::-1]
            idel = np.zeros(ind.size, dtype = bool)
            for i in range(ind.size):
                if not idel[i]:
                    idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                           & (trends[ind[i]] > trends[ind])
                    idel[i] = 0

            ind = np.sort(ind[~idel])

        return ind

    def randomCheck(self, seq, pthre = 0.05):
        """
        We use chi square test to test the randomness of a sequence by
        looking at the conversion frequency between neighbors in the sequence.
        
        Parameters
        ----------
        seq : str
            A string containing only '1' or '0'. (e.g. '101000101')
        
        pthre : float, 0-1
            Significance level of the hypothesis tests.
        
        Returns
        -------
        reject : bool
            True if we should reject the null hypothesis (the sequence is
            generated randomly) under the selected significance level.
        """
        from scipy.stats import chisquare

        pairwise = zip(seq[1:], seq[:-1])
        d = collections.defaultdict(int)
        for k in pairwise:
            d[k] += 1
        for k in [('0','0'), ('0','1'), ('1','0'), ('1','1')]:
            if not k in d:
                d[k] = 0
        obs = np.array(list(d.values()))
        exp = np.ones_like(obs) * obs.mean()

        _, pval =  chisquare(obs, exp)
        reject = pval<=pthre

        return reject

    def oriWindow(self, P):
        """
        Estimate the most appropriate window size for current bin to best
        capture the local interaction bias direction.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.detectPeaks : detect peaks given a 1-D array
        tadlib.hitad.chromLev.Chrom.randomCheck : randomness test for a two-valued
                                                  (0-1) sequence
        """
        noise = P == 0
        check = noise[:20]
        noiselevel = check.sum() / check.size
        if noiselevel > 0.6:
            return 0

        indices = np.arange(1, P.size+1)

        m = [P < 0, # Downstream bias
             P > 0] # Upstream bias

        trends_1 = m[0].cumsum().astype(float) / np.arange(1, m[0].size + 1)
        trends_2 = m[1].cumsum().astype(float) / np.arange(1, m[1].size + 1)

        inds = [self.detectPeaks(trends_1, 0.5, 5),
                self.detectPeaks(trends_2, 0.5, 5)]
        pool = {}
        for i in [0, 1]:
            for p in inds[i]:
                pool[p] = i

        for p in sorted(pool):
            seq = ''.join([str(int(i)) for i in m[pool[p]][:(p+1)]])
            tmp = indices[p]+self._rm+1
            if tmp >= self.minsize: # hasn't been systematically tested
                if self.randomCheck(seq):
                    return tmp
        
        return self._dw

    def minWindows(self, start, end, maxw):
        """
        Estimate best window size for each bin of a given range.
        
        Parameters
        ----------
        start, end : int
            Specify range of the bin indices.
        
        maxw : int
            Maximum allowable window size.
        
        Attributes
        ----------
        windows : 1-D numpy.ndarray, int32
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.oriWindow : Window size estimation for
                                                a single bin.
        """

        start += self.chromLen; end += self.chromLen
        self.windows = np.zeros(end - start, dtype = np.int32)
        for i in range(start, end):
            down = self.rawMatrix[i, i:(i+maxw)].toarray().ravel()
            up = self.rawMatrix[(i-maxw+1):(i+1), i].toarray().ravel()[::-1]
            down[:self._rm+1] = 0; up[:self._rm+1] = 0
            diff = up - down
            ws = self.oriWindow(diff[self._rm+1:])
            self.windows[i-start] = ws
    
    def calIS(self, idx, w=10):

        idx = idx + self.chromLen
        sub = self.rawMatrix[idx-w:idx, idx+1:idx+w+1]

        return sub.mean()
    
    def preciseBound(self, byregion):

        for r in byregion:
            tmp = byregion[r]
            for d in tmp:
                si = d[0]//self.res
                ini = np.inf
                for i in range(max(si-1,0), min(si+2,self.chromLen)):
                    IS = self.calIS(i)
                    if IS < ini:
                        d[0] = i * self.res
                        ini = IS

                ei = d[1]//self.res
                ini = np.inf
                for i in range(max(ei-1,0), min(ei+2,self.chromLen)):
                    IS = self.calIS(i)
                    if IS < ini:
                        d[1] = i * self.res
                        ini = IS     

    def calDI(self, windows, start):
        """
        Calculate Directionality Index (DI) for each bin with adaptive
        window size.
        
        Parameters
        ----------
        windows : 1-D numpy.ndarray, int32
            Returned by :py:meth:`tadlib.hitad.chromLev.Chrom.minWindows`.
        start : int
            Starting bin index, the window size of which is taken from the
            1st place of *windows*.
    
        Attributes
        ----------
        DIs : 1-D numpy ndarray, float
            Calculated adaptive DI array, which has the same size as the
            input *windows*.
        """

        start = start + self.chromLen

        self.DIs = np.zeros(len(windows))
        for i in range(start, start + len(windows)):
            w = windows[i-start]
            if w == 0:
                w = self._dw
            down = self.rawMatrix[i, i:(i+w)].toarray().ravel()
            up = self.rawMatrix[(i-w+1):(i+1), i].toarray().ravel()[::-1]
            down = down[self._rm+1:]; up = up[self._rm+1:]
            tmp = self._binbias(up, down)
            if tmp != 0:
                self.DIs[i-start] = tmp
            else:
                if w < self._dw:
                    w = self._dw
                    down = self.rawMatrix[i, i:(i+w)].toarray().ravel()
                    up = self.rawMatrix[(i-w+1):(i+1), i].toarray().ravel()[::-1]
                    down = down[self._rm+1:]; up = up[self._rm+1:]
                    self.DIs[i-start] = self._binbias(up, down)

        # trim outliers
        lthre = np.percentile(self.DIs[self.DIs<0], 0.1)
        hthre = np.percentile(self.DIs[self.DIs>0], 99.9)
        self.DIs[self.DIs<lthre] = lthre
        self.DIs[self.DIs>hthre] = hthre

    def _binbias(self, up, down):

        bias = 0

        zeromask = (up != 0) & (down != 0)
        if zeromask.sum() >= 5:
            up = up[zeromask]; down = down[zeromask]

        if up.size <= 1:
            return bias
        
        upmean = up.mean(); downmean = down.mean()
        SD_1 = np.sum((up - upmean) ** 2) / (up.size * (up.size - 1))
        SD_2 = np.sum((down - downmean) ** 2) / (down.size * (down.size - 1))
        SD_pool = np.sqrt(SD_1 + SD_2)
        if SD_pool != 0:
            bias = (upmean - downmean) / SD_pool

        return bias

    def splitChrom(self, DIs):
        """
        Split a chromosome into gap-free regions. HMM learning and domain
        identification procedures will be performed on these regions separately.
        
        Parameters
        ----------
        DIs : 1-D numpy ndarray, float
            Adaptive DI array of the whole chromosome. Generally, we detect
            runs of zeros in the array as gaps, which will be cut off the
            chromosome, making entire chromosome pieces of gap-free regions.
        
        Attributes
        ----------
        chromRegions : dict, {(start,end):DIs[start:end]}
            The keys are gap-free regions, and the values are corresponding
            adaptive DI pieces.
        gapbins : set
            Set of bins (in base-pair unit) located in gap regions.
        """
        # minregion and maxgaplen are set intuitively
        maxgaplen = max(100000 // self.res, 5)
        minregion = maxgaplen * 2

        valid_pos = np.where(DIs != 0)[0]
        gapsizes = valid_pos[1:] - valid_pos[:-1]
        endsIdx = np.where(gapsizes > (maxgaplen + 1))[0]
        startsIdx = endsIdx + 1

        chromRegions = {}
        for i in range(startsIdx.size - 1):
            start = valid_pos[startsIdx[i]]
            end = valid_pos[endsIdx[i + 1]] + 1
            if end - start > minregion:
                chromRegions[(start, end)] = DIs[start:end]

        if startsIdx.size > 0:
            start = valid_pos[startsIdx[-1]]
            end = valid_pos[-1] + 1
            if end - start > minregion:
                chromRegions[(start, end)] = DIs[start:end]
            start = valid_pos[0]
            end = valid_pos[endsIdx[0]] + 1
            if end - start > minregion:
                chromRegions[(start, end)] = DIs[start:end]

        if not startsIdx.size:
            if valid_pos.size > 0:
                start = valid_pos[0]
                end = valid_pos[-1]
                if end - start > minregion:
                    chromRegions[(start, end)] = DIs[start:end]

        gapmask = np.ones(DIs.size, bool)
        for r in chromRegions:
            gapmask[r[0]:r[1]] = 0
        gapbins = set(np.where(gapmask)[0]*self.res)
        
        self.regionDIs, self.gapbins = chromRegions, gapbins

    def viterbi(self, seq):
        """
        Find the most likely hidden state series using the viterbi algorithm.
        
        Parameters
        ----------
        seq : 1-D numbpy ndarray, float
            Adaptive DI array for any region.
        
        Returns
        -------
        path : list
            List of hidden state labels. Has the same length as the input
            *seq*.
        
        """
        path = [int(s.name) for i, s in self.hmm.viterbi(seq)[1][1:-1]]

        return path

    def _getBounds(self, path, junctions=['30']):
        """
        Call boundary sites from hidden state series. By default, these
        transition modes will be detected as boundaries: "no bias(2) >
        domain start(0)", "domain end(4) > domain start(0)", and "domain end(4)
        > no bias(2)".
        
        Parameters
        ----------
        path : list
            Hidden state series returned by :py:meth:`tadlib.hitad.chromLev.Chrom.viterbi`.
        junctions : list of strings
            Boundary definitions by using state transition modes.
            (Default: ['20','40','42'])
        
        Returns
        -------
        bounds : 1-D numpy ndarray, int
            Detected boundary positions in ascending order. 0 and len(*path*)
            will always be included.
        """
        pathseq = ''.join(map(str, path))
        pieces = [pathseq]
        for junc in junctions:
            gen = []
            for seq in pieces:
                tmp = seq.split(junc)
                if len(tmp) == 1:
                    gen.extend(tmp)
                else:
                    gen.append(tmp[0]+junc[0])
                    for s in tmp[1:-1]:
                        gen.append(junc[1]+s+junc[0])
                    gen.append(junc[1]+tmp[-1])

            pieces = gen
            
        bounds = np.r_[0, np.r_[list(map(len, pieces))]].cumsum()

        return bounds
        
    def pipe(self, seq, start):
        """
        Transform an observed sequence into a list of domains.
        
        Parameters
        ----------
        seq : 1-D numbpy ndarray, float
            Adaptive DI array for any region.
        start : int
            Chromosome bin index of the *seq* start.
        
        Returns
        -------
        domains : list
            List of domains in the format ``[start bin, end bin, noise level,
            hierarchical level]``.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.refNoise : Calculate the noise level of
                                               a given domain
        tadlib.hitad.aligner.BoundSet : where the meanings of the hierarchical
                                        level labels are explained in detail.
        """
        # bin-level domain (not base-pair-level domain!)
        bounds = self._getBounds(self.viterbi(seq), junctions=['30'])
        pairs = [[bounds[i], bounds[i+1]] for i in range(len(bounds)-1)]
        domains = []
        for b in pairs:
            # start, end, noise level, hierarchical level
            tmp = [b[0]+start, b[1]+start, 0, 0]
            domains.append(tmp)

        return domains

    def minCore(self, regionDIs):
        """
        Output domain list for each gap-free region.
        
        Parameters
        ----------
        regionDIs : dict
            Gap-free regions and corresponding adaptive DI arrays.
        
        Returns
        -------
        minDomains : dict
            Gap-free regions and corresponding identified bottom domain list.
            Different from :py:meth:`tadlib.hitad.chromLev.Chrom.pipe`, the
            start and the end of a domain are in base-pair unit.

        """
        tmpDomains = {}
        for region in sorted(regionDIs):
            seq = regionDIs[region]
            domains = self.pipe(seq, region[0])
            cr = (region[0]*self.res, region[1]*self.res)
            tmpDomains[cr] = []
            for domain in domains:
                domain[0] = domain[0] * self.res
                domain[1] = domain[1] * self.res
                domain[2] = self.refNoise(domain)
                tmpDomains[cr].append(domain)
        
        minDomains = self._orifilter(tmpDomains)

        return minDomains

    def getDomainList(self, byregion):
        """
        Combine by-region domains into a single list.
        
        Parameters
        ----------
        byregion : dict
            The keys are tuples representing gap-free regions of the chromosome,
            and the values are corresponding identified domain lists.
        
        Returns
        -------
        DomainList : list
            A merged domain list of all regions
        """
        DomainList = []
        for region in sorted(byregion):
            DomainList.extend(byregion[region])

        return DomainList

    def _orifilter(self, oriDomains):
        """
        Perform size filtering on the input domain lists.
        
        Parameters
        ----------
        oriDomains : dict
            The keys are tuples representing gap-free regions of the chromosome,
            and the values are corresponding identified domain lists. Start
            and end of the domain should be in base-pair unit.
        
        Returns
        -------
        filtered : dict
            Pairs of gap-free regions and corresponding filtered domain lists.
        """
        filtered = {}
        for region in oriDomains:
            tmplist = []
            for d in oriDomains[region]:
                if d[1] - d[0] >= (self.minsize*self.res):
                    tmplist.append(d)
            if len(tmplist):
                filtered[region] = tmplist
        return filtered

    def iterCore(self, minDomains, tmpDomains):
        """
        Calculate the mismatch ratio for the input two domain lists. Return
        1 if *minDomains* is empty.
        
        :py:meth:`tadlib.hitad.chromLev.Chrom.oriIter` uses this method to
        determine whether to break the iteration.
        
        Parameters
        ----------
        minDomains : dict
            Target domains calculated by the last loop.
        tmpDomains : dict
            Query domains returned by the current loop.
        
        Returns
        -------
        tol : float
            Mismatch ratio.
        """
        tmplist = self.getDomainList(copy.deepcopy(minDomains))
        reflist = []
        for refd in tmplist:
            reflist.append([self.chrom, refd[0], refd[1], 0])
        if not len(reflist):
            tol = 1
        else:
            tmplist = self.getDomainList(copy.deepcopy(tmpDomains))
            alignlist = []
            for alignd in tmplist:
                alignlist.append([self.chrom, alignd[0], alignd[1], 0])
            Ref = DomainSet('ref', reflist, self.res)
            Align = DomainSet('align', alignlist, self.res)
            worker = DomainAligner(Ref, Align)
            worker.align('ref','align')
            count = len(worker.conserved('ref','align'))
            tol = 1 - count / len(Ref.Domains)

        return tol

    def oriIter(self, minDomains):
        """
        Iteratvely approximate adaptive window sizes and return the final
        bottom domain list which will be used in subsequent procedures. For
        each loop, window sizes are updated according to the latest bottom
        domains and next loop will re-run the identification pipeline using
        new window sizes. The iteration terminates if domains between two
        consecutive loops are very similar (estimated by
        :py:meth:`tadlib.hitad.chromLev.Chrom.iterCore`).
        
        Parameters
        ----------
        minDomains : dict
            Initial domains served as the target domain list for
            :py:meth:`tadlib.hitad.chromLev.Chrom.iterCore` at the first
            iteration. We set it empty in our pipeline.
        
        Attributes
        ----------
        minDomains : dict
            The keys are tuples representing gap-free regions of the chromosome,
            and the values are corresponding identified bottom domain lists.
            Start and end of the domain are in base-pair unit.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.calDI : calculate DI values according to
                                            input window sizes
        tadlib.hitad.chromLev.Chrom.iterCore : estimate the degree of divergence
                                               between two domain lists
        """
        for n_i in range(5):
            tmpDomains = self.minCore(self.regionDIs)
            tol = self.iterCore(minDomains, tmpDomains)
            minDomains = tmpDomains
            for region in minDomains:
                for d in minDomains[region]:
                    ds = d[0]//self.res; de = d[1]//self.res
                    Len = de - ds
                    ws = np.max(np.r_['0,2,1', np.arange(1,Len+1), np.arange(Len,0,-1)],
                                axis=0)
                    self.windows[ds:de] = ws
            self.calDI(self.windows, 0)
            self.splitChrom(self.DIs)
            if tol < 0.05:
                break

        self.minDomains = minDomains
    
    def callDomains(self):
        """
        Direct API for our hierarchical domain identification pipeline:
        
        - Adaptively estimate window size for each bin.
          (:py:meth:`tadlib.hitad.chromLev.Chrom.minWindows`)
        - Calculate adaptive DIs. (:py:meth:`tadlib.hitad.chromLev.Chrom.calDI`)
        - Iteratively correct adaptive window size and bottom boundary positions.
          (:py:meth:`tadlib.hitad.chromLev.Chrom.oriIter`)
        """
        self.minWindows(0, self.chromLen, self._dw)
        self.calDI(self.windows, 0)
        self.splitChrom(self.DIs)
        self.oriIter({})
        #self.preciseBound(self.minDomains)

        self.domains = self.getDomainList(self.minDomains)

        self._state = 'Completed'
        

    def getSelfMatrix(self, start, end):
        """
        Return the contact matrix of any given region.
        
        Parameters
        ----------
        start, end : int
            The region interval in base-pair unit.
        
        Returns
        -------
        Matrix : 2-D numpy ndarray, float
            Sub contact matrix.
        """
        startidx = start // self.res + self.chromLen
        endidx = end // self.res + self.chromLen

        Matrix = self.rawMatrix[startidx:endidx, startidx:endidx].toarray()

        x, y = np.nonzero(Matrix)
        Matrix[y, x] = Matrix[x, y]

        return Matrix

    def refNoise(self, domain):
        """
        Return noise level of a domain, which is simply defined as the zero
        entry ratio of the contact matrix.
        """
        if domain[1] - domain[0] < self.res*self.minsize:
            return 1

        matrix = self.getSelfMatrix(domain[0], domain[1])

        total = matrix.size - np.arange(len(matrix), len(matrix)-self._rm-1, -1).sum()*2 +\
                len(matrix)

        if total < 5:
            return 1
        else:
            nx, ny = np.nonzero(matrix)
            mask = np.abs(ny-nx) > self._rm
            sigNum = mask.sum() # Number of nonzero entries apart from diag
            return 1-sigNum/total


