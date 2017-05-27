# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:47:34 2016

@author: wxt
"""

from __future__ import division
import logging, copy, collections
import numpy as np
from scipy import sparse
from aligner import BoundSet, DomainSet, DomainAligner, hierFormat, Container

np.seterr(divide = "ignore")

## Customize the logger
log = logging.getLogger(__name__)

class Chrom(object):
    """
    *Chrom* is defined to:
    
    - Hold Hi-C data within a certain chromosome
    - Identify hierarchical domains in 4 steps: 1.Calculate adaptive DIs for
      all bins of the considered chromosome. 2.Identify original candidate
      bounds by 5-state Gaussian mixture Hidden Markov Model using adaptive
      DIs as input. 3.Select TAD bounds from candidate bounds. 4.Recursively
      identify inner domain structures of each TAD.
    - Visualize any region of the chromosome. Hierarchical domains will be
      plotted as boxes along with the diagonal of the heatmap, and adaptive
      DI track will be placed on top of the heatmap.
    
    Parameters
    ----------
    chrom : str
        Chromosome label.
    
    res : int
        Resolution of the Hi-C data in base-pair unit.
    
    hicdata : numpy.ndarray
        Hi-C data stored in our customized Numpy Structured Array. The
        structured array is defined in :py:class:`tadlib.hitad.genomeLev.Genome`,
        and it has 3 fields named "bin1", "bin2" and "IF", respectively.
    
    replabel : str
        Biological replicate label.
    
    maxapart : int
        Maximum allowable TAD size in base-pair unit. (Default: 4000000)
    
    Attributes
    ----------
    chrom : str
        Chromosome label.
    
    res : int
        Resolution in base-pair unit.
    
    chromLen : int
        Total bin number of the chromosome. Obtained from Hi-C data.
    
    rawMatrix : sparse matrix in Compressed Sparse Row format
        CSR sparse matrix is used to extract Hi-C data by slicing conveniently
        while guarantee low memory overhead.
        
    """

    defaultwindow = 2000000
    minsize = 5

    def __init__(self, chrom, res, hicdata, replabel, maxapart=4000000):

        self.chrom = chrom
        self.res = res
        self.maxapart = maxapart
        self.replabel = replabel
        self._rm = 1
        self._dw = self.defaultwindow // res
        self._mw = self.maxapart // res

        x = hicdata['bin1']; y = hicdata['bin2']; IF = hicdata['IF']
        self.chromLen, self.rawMatrix = self._genSparseMatrix(x, y, IF)

        del x, y, IF, hicdata

        self._state = 'Submitted'

    def _genSparseMatrix(self, x, y, IF):

        indices = np.r_['0,2,1', x, y]
        indices.sort(axis = 0)
        x, y = indices
        chromLen = y.max() + 1
        extendLen = chromLen + 2*chromLen
        rawMatrix = sparse.csr_matrix((IF, (x + chromLen, y + chromLen)),
                                      shape = (extendLen, extendLen))

        return chromLen, rawMatrix

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
        from itertools import izip
        from scipy.stats import chisquare

        pairwise = izip(seq[1:], seq[:-1])
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
        
        Parameters
        ----------
        P : 1-D numpy.ndarray
            Calculated by the formula P\ :sub:`i`\ (k) = M\ :sub:`i-k,k` -
            M\ :sub:`i,i+k`, where *M* represents the contact matrix, *i*
            represents current bin index, and *k* indicates bin-level genomic
            distances. This is an intuitive definition of interaction bias
            for each genomic distance: if the sign of an entry is positive,
            then bin *i* is more likely to contact with the upstream bin
            *i-k*; otherwise it's more inclined to contact with the downstream
            bin *i+k*.
        
        Design
        ------
        Because of the presence of local domains, we can see runs of values
        in *P* with consistent signs. The object of this method is to determine
        the cut of the first run. For the beginning of a domain, the values
        in the first run tend to be negative, and for the end of a domain,
        they are almost positive.
        
        The cut sites of the runs can be approximated by detecting peaks on
        the trend of positive/negative value ratio under increasing *k*. And
        the existence of runs under certain *k* can be determined by checking
        the randomness of the sequence.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.detectPeaks : detect peaks given a 1-D array
        tadlib.hitad.chromLev.Chrom.randomCheck : randomness test for a two-valued
                                                  (0-1) sequence
        """
        noise = P == 0
        check = noise[:50]
        noiselevel = check.sum() / check.size
        if noiselevel > 0.2:
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
        Estimate best window size for every bin of a given range.
        
        Parameters
        ----------
        start, end : int
            Specify range of the bins.
        
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
        for i in xrange(start, end):
            down = self.rawMatrix[i, i:(i+maxw)].toarray().ravel()
            up = self.rawMatrix[(i-maxw+1):(i+1), i].toarray().ravel()[::-1]
            down[:self._rm+1] = 0; up[:self._rm+1] = 0
            diff = up - down
            ws = self.oriWindow(diff[self._rm+1:])
            self.windows[i-start] = ws

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
    
        Returns
        -------
        DIs : 1-D numpy ndarray, float
            Calculated adaptive DI array, which has the same size as the
            input *windows*.
        """

        start = start + self.chromLen

        self.DIs = np.zeros(len(windows))
        for i in xrange(start, start + len(windows)):
            w = windows[i-start]
            if w > 0:
                down = self.rawMatrix[i, i:(i+w)].toarray().ravel()
                up = self.rawMatrix[(i-w+1):(i+1), i].toarray().ravel()[::-1]
                down = down[self._rm+1:]; up = up[self._rm+1:]
                self.DIs[i-start] = self._binbias(up, down)

        # trim outliers
        lthre = np.percentile(self.DIs, 0.1)
        hthre = np.percentile(self.DIs, 99.9)
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

    def oriHMMParams(self):
        """
        Set initial parameters for the Hidden Markov Model (HMM).
        
        Attributes
        ----------
        HMMParams : dict
            Has 3 keys: "A", state transition matrix, "B" (emission probabilities),
            specifying parameters (Means, Variances, Weights) of the mixture
            Gaussian distributions for each hidden state, and "pi", indicating
            the hidden state weights. This dict will be updated after learning
            procedure.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.paramFromModel : Learn HMM parameters
                                                     from sequence data.
        """

        if not 'HMMParams' in self.__dict__:
            # Five Hidden States:
            # 0--start, 1--downstream, 2--no bias, 3--upstream, 4--end
            A = [[0., 1., 0., 0., 0.],
                 [0., 0.4, 0.3, 0.3, 0.],
                 [0.05, 0., 0.5, 0.45, 0.],
                 [0., 0., 0., 0.5, 0.5],
                 [0.99, 0., 0.01, 0., 0.]]
            pi = [0.05, 0.3, 0.3, 0.3, 0.05]
            numdists = 3 # Three-distribution Gaussian Mixtures
            W = 1. / numdists
            var = 7.5 / (numdists - 1)
            means = [[], [], [], [], []]
            for i in range(numdists):
                means[4].append(i * 7.5 / ( numdists - 1 ) + 2.5)
                means[3].append(i * 7.5 / ( numdists - 1 ))
                means[2].append((i - (numdists-1)/2) * 7.5 / (numdists - 1))
                means[1].append(-i * 7.5 / ( numdists - 1 ))
                means[0].append(-i * 7.5 / ( numdists - 1 ) - 2.5)

            B = [[means[i], [var for j in range(numdists)], [W for k in range(numdists)]] for i in range(len(pi))]

            self.HMMParams = {'A': A, 'B': B, 'pi': pi}

    def splitChrom(self, DIs):
        """
        Split a chromosome into gap-free regions. HMM learning and domain
        identification procedure will be performed on these region separately.
        
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
        maxgaplen = 100000 // self.res
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
            start = valid_pos[0]
            end = valid_pos[-1]
            if end - start > minregion:
                chromRegions[(start, end)] = DIs[start:end]

        if not len(chromRegions):
            raise ValueError('Empty DI sequences for HMM training')

        gapmask = np.ones(DIs.size, bool)
        for r in chromRegions:
            gapmask[r[0]:r[1]] = 0
        gapbins = set(np.where(gapmask)[0]*self.res)
        
        self.regionDIs, self.gapbins = chromRegions, gapbins

    def paramFromModel(self, seqs):
        """
        Train HMM model on *seqs* using `GHMM <http://ghmm.sourceforge.net/index.html>`_.
        Update *HMMParams* attribute after training.
        
        Parameters
        ----------
        seqs : list of list
            Each element is the adaptive DI list of a gap-free region.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.oriHMMParams : Set initial parameters
        tadlib.hitad.chromLev.Chrom.splitChrom : where gap-free regions are
                                                 returned.
        """
        import ghmm

        HMMParams = getattr(self, 'HMMParams')

        F = ghmm.Float()
        trainset = ghmm.SequenceSet(F, seqs)

        # Training ...
        model = ghmm.HMMFromMatrices(F, ghmm.GaussianMixtureDistribution(F),
                                     **HMMParams)
        model.baumWelch(trainset)

        state_num = len(HMMParams['B'])
        numdists = len(HMMParams['B'][0][0])

        A = np.zeros((state_num, state_num))
        B = np.zeros((state_num, 3, numdists))
        pi = np.zeros(state_num)

        for i in range(state_num):
            for j in range(state_num):
                A[i, j] = model.getTransition(i, j)
            for j in range(numdists): # num of distributions
                temp = model.getEmission(i, j)
                for k in range(3): # Mean, var, weight
                    B[i, k, j] = temp[k]

            pi[i] = model.getInitial(i)

        HMMParams['A'] = A
        HMMParams['B'] = B
        HMMParams['pi'] = pi

    def learning(self, regionDIs):
        """
        Prepare training data and learn HMM model parameters.
        
        Parameters
        ----------
        regionDIs : dict
            Gap-free regions and corresponding adaptive DI arrays. Returned
            by :py:meth:`tadlib.hitad.chromLev.Chrom.splitChrom`.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.paramFromModel : core part for HMM learning
                                                     and parameter updating.
        """
        import random

        self.oriHMMParams()

        seqs = []
        for region in regionDIs:
            seqs.append(list(regionDIs[region]))

        random.shuffle(seqs)

        self.paramFromModel(seqs)

    def viterbi(self, seq):
        """
        Find the most likely hidden state series given the observed *seq*
        using the viterbi algorithm.
        
        Parameters
        ----------
        seq : 1-D numbpy ndarray, float
            Adaptive DI array for any region.
        
        Returns
        -------
        path : list
            List of hidden state labels. Has the same length as the input
            *seq*.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.paramFromModel : parameter learning
                                                     using Baum-Welch algorithm
        """
        from scipy.stats import norm

        HMMParams = getattr(self, 'HMMParams')

        num_states = HMMParams['pi'].shape[0]
        numdists = HMMParams['B'].shape[2]
        seq_len = seq.shape[0]
        costs = np.zeros((num_states, seq_len))
        paths = np.zeros((num_states, seq_len - 1), np.int8)
        transition_costs = -np.log(HMMParams['A'])

        for i in range(num_states):
            for j in range(numdists):
                costs[i,:] += norm.pdf(seq, loc = HMMParams['B'][i, 0, j],
                                       scale = HMMParams['B'][i, 1, j] ** 0.5) * \
                                       HMMParams['B'][i, 2, j]
        costs = -np.log(costs)
        costs[:,0] -= np.log(HMMParams['pi'])
        for i in range(seq_len - 1):
            for j in range(num_states):
                min_cost = costs[0, i] + transition_costs[0, j]
                min_state = 0
                for k in range(1, num_states):
                    next_cost = costs[k, i] + transition_costs[k, j]
                    if next_cost < min_cost:
                        min_cost = next_cost
                        min_state = k
                costs[j, i + 1] += min_cost
                paths[j, i] = min_state
        for j in range(num_states):
            min_cost = costs[0, -1] + transition_costs[0, j]
            min_state = 0
            for k in range(1, num_states):
                next_cost = costs[k, -1] + transition_costs[k, j]
                if next_cost < min_cost:
                    min_cost = next_cost
                    min_state = k
        
        min_state = 0
        for j in range(1, num_states):
            if costs[j, -1] < costs[min_state, -1]:
                min_state = j

        path = [min_state]
        for i in range(paths.shape[1])[::-1]:
            path = [paths[path[0], i]] + path

        return path

    def _getBounds(self, path, junctions=['20','40','42']):
        """
        Call boundary sites from hidden state series. By default, these
        transition modes will be detected as boundaries: "no bias(2) 🡒
        domain start(0)", "domain end(4) 🡒 domain start(0)", and "domain end(4)
        🡒 no bias(2)".
        
        Parameters
        ----------
        path : list
            Hidden state series returned by :py:meth:`tadlib.hitad.chromLev.Chrom.viterbi`.
        junctions : list of strings
            Boundary definitions using state transition modes.
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
            
        bounds = np.r_[0, np.r_[map(len, pieces)]].cumsum()

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
            The elements are also lists in format ``[start bin, end bin,
            noise level, hierarchical level]``.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.refNoise : where domain noise level is
                                               calculated
        tadlib.hitad.aligner.BoundSet : where the meanings of the hierarchical
                                        level labels are explained in detail.
        """
        # bin-level domain (not base-pair-level domain!)
        bounds = self._getBounds(self.viterbi(seq))
        pairs = [[bounds[i], bounds[i+1]] for i in xrange(len(bounds)-1)]
        domains = []
        for b in pairs:
            # start, end, noise level, hierarchical level
            tmp = [b[0]+start, b[1]+start, 0, 0]
            domains.append(tmp)

        return domains

    def minCore(self, regionDIs):
        """
        Learner and Caller. Take adaptive DI arrays of all gap-free regions,
        learn HMM parameters, retrieve hidden state series, and finally
        output identified domain list for each region.
        
        In our implementation, this method is mainly used as a core for
        :py:meth:`tadlib.hitad.chromLev.Chrom.oriIter`.
        
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
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.oriIter : iteratively approximate adaptive
                                              window sizes and return final
                                              bottom domain lists.
        tadlib.hitad.chromLev.Chrom.splitChrom : by which the input learning
                                                 data can be calculated
        tadlib.hitad.chromLev.Chrom.pipe : take an adaptive DI array and return
                                           a list of domains after learning
                                           procedure has been completed
        """
        self.learning(regionDIs)
        minDomains = {}
        for region in sorted(regionDIs):
            seq = regionDIs[region]
            domains = self.pipe(seq, region[0])
            cr = (region[0]*self.res, region[1]*self.res)
            minDomains[cr] = []
            for domain in domains:
                domain[0] = domain[0] * self.res
                domain[1] = domain[1] * self.res
                domain[2] = self.refNoise(domain)
                minDomains[cr].append(domain)

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
        Perform size filtering on input domain lists.
        
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
        Calculate the mismatch ratio for input two domain lists. Return
        1 if *minDomains* is empty.
        
        :py:meth:`tadlib.hitad.chromLev.Chrom.oriIter` uses this method to
        determine whether to break the iteration.
        
        Parameters
        ----------
        minDomains : dict
            Target domain list calculated by the last loop.
        tmpDomains : dict
            Query domain list of current loop.
        
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
        each loop, window sizes are updated according to identified bottom
        domains and next loop will re-run the identification pipeline using
        new window sizes. The iteration terminates if domains between two
        loops are very similar (estimated by :py:meth:`tadlib.hitad.chromLev.Chrom.iterCore`).
        
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
        tadlib.hitad.chromLev.Chrom.minCore : learning HMM parameters and
                                              identify bottom domains by using
                                              updated adaptive DIs
        tadlib.hitad.chromLev.Chrom.iterCore : estimate the degree of divergence
                                               between two domain lists
        """
        for n_i in range(5):
            log.debug('Chrom %s (res %dK, %s): calling bottom domains, %d iteration, start ...',
                      self.chrom, self.res//1000, self.replabel, n_i+1)
            oriDomains = self.minCore(self.regionDIs)
            tmpDomains = self._orifilter(oriDomains)
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
            log.debug('Chrom %s (res %dK, %s): calling bottom domains, %d iteration, unaligned ratio %.3g%%',
                      self.chrom, self.res//1000, self.replabel, n_i+1, tol*100)
            if tol < 0.05:
                break

        self.minDomains = minDomains

    def fineDomain(self):
        """
        Identify hierarchical domains within each TAD.
        
        Attributes
        ----------
        hierDomains : dict
            The keys are tuples representing gap-free regions of the chromosome,
            and the values are corresponding identified hierarchical domain lists.
            Start and end of the domain are in base-pair unit.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.maxCore : identify TADs using bottom
                                              domain lists
        tadlib.hitad.chromLev.Chrom.subDomains : resolve domain hierarchy
                                                 within a given TAD
        """
        self.hierDomains = {}
        for region in sorted(self.maxDomains):
            rdomains = self.maxDomains[region]
            if not len(rdomains):
                continue
            hdomains = []
            for i in range(len(rdomains)):
                top = rdomains[i][:2]
                subdomains = {}
                self.subDomains(top, self.minDomains[region], subdomains=subdomains)
                if len(subdomains)==1:
                    noise = self.refNoise(top)
                    hdomains.append(top+[noise,0])
                else:
                    for d in subdomains:
                        noise = self.refNoise(list(d))
                        hdomains.append(list(d)+[noise,subdomains[d]])
            hdomains.sort()
            self.hierDomains[region] = hdomains

    def callDomain(self):
        """
        Direct API for our hierarchical domain identification pipeline:
        
        - Adaptively estimate window size for each bin.
          (:py:meth:`tadlib.hitad.chromLev.Chrom.minWindows`)
        - Calculate adaptive DIs. (:py:meth:`tadlib.hitad.chromLev.Chrom.calDI`)
        - Perform HMM learning, identify initial bottom domains, and iteratively
          correct adaptive window size and bottom boundary positions.
          (:py:meth:`tadlib.hitad.chromLev.Chrom.oriIter`)
        - Identify TADs based on bottom domains.
          (:py:meth:`tadlib.hitad.chromLev.Chrom.maxCore`)
        - Resolve domain hierarchy within each TAD.
          (:py:meth:`tadlib.hitad.chromLev.Chrom.fineDomain`)
        """
        log.debug('Chrom %s (res %dK, %s): domain calling start ...',
                  self.chrom, self.res//1000, self.replabel)

        self.minWindows(0, self.chromLen, self._mw)
        self.calDI(self.windows, 0)
        self.splitChrom(self.DIs)
        self.oriIter({})
        log.debug('Chrom %s (res %dK, %s): identify TADs ...',
                  self.chrom, self.res//1000, self.replabel)
        self.maxCore()
        log.debug('Chrom %s (res %dK, %s): resolve domain hierarchy within each TAD ...',
                  self.chrom, self.res//1000, self.replabel)
        self.fineDomain()
        log.debug('Chrom %s (res %dK, %s): domain calling completed',
                  self.chrom, self.res//1000, self.replabel)

        self._state = 'Completed'

    def idxmatch(self, domain):
        """
        Pair interactions of the given domain with the upstream and downstream
        interactions under the same genomic distance.
        
        Parameters
        ----------
        domain : [start, end]
            Domain interval in base-pair unit.
        
        Returns
        -------
        cur_store : tuple, (x-coordinates, y-coordinates, interaction frequencies)
            Interactions within the given domain.
        up_store : tuple, (x-coordinates, y-coordinates, interaction frequencies)
            Corresponding upstream interctions.
        down_store : tuple, (x-coordinates, y-coordinates, interaction frequencies)
            Corresponding downstream interactions.
        
        """
        ds = domain[0] // self.res + self.chromLen
        de = domain[1] // self.res + self.chromLen
        if (de - ds) < 2:
            return

        Rawx = np.tile(np.r_['0,2,0', ds:de], (1, de - ds))
        Rawy = np.tile(np.r_['0,2,1', ds:de], (de - ds, 1))
        d = Rawy - Rawx

        Ux = Rawx - d
        Uy = Rawx
        Dx = Rawy
        Dy = Dx + d

        curInters = self.rawMatrix[Rawx, Rawy].toarray()
        upInters = self.rawMatrix[Ux, Uy].toarray()
        downInters = self.rawMatrix[Dx, Dy].toarray()

        cur_store = (Rawx, Rawy, curInters)
        up_store = (Ux, Uy, upInters)
        down_store = (Dx, Dy, downInters)

        return up_store, cur_store, down_store

    def stablescore(self, domain):
        """
        Calculate TAD score for the given domain. 
        
        Parameters
        ----------
        domain : [start, end]
            Domain interval in base-pair unit.
        
        Returns
        -------
        inout : float
            TAD score defined as the enrichment between intra-domain interaction
            frequencies and inter-domain interaction frequencies controlling
            for the impact of genomic distance.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.idxmatch : Pair intra-domain and inter-domain
                                               interactions with the same
                                               genomic distance
        """
        ori = domain[0] // self.res + self.chromLen
        end = domain[1] // self.res + self.chromLen
        P = self.idxmatch(domain)
        if P is None:
            return 0, None, None
        ## Compare inside-domain interactions with outside ones, larger, better
        inside = np.r_[[]]; outside = np.r_[[]]
        # Pair with upstream interactions
        Umask = ((2*P[1][0] - ori) <= P[1][1]) & (P[1][1] - P[1][0] > self._rm) &\
                (P[1][2] != 0) & (P[0][2] != 0)
        inside = np.r_[inside, P[1][2][Umask]]
        outside = np.r_[outside, P[0][2][Umask]]
        # Pair with downstream interactions
        Dmask = ((2*P[2][0] - ori) >= P[2][1]) & (P[2][1] - P[2][0] > self._rm) &\
                (P[2][1] >= end) & (P[1][2] != 0) & (P[2][2] != 0)
        inside = np.r_[inside, P[1][2][Dmask]]
        outside = np.r_[outside, P[2][2][Dmask]]
        inout = 0
        if inside.size > 5:
            diff = (inside - outside) / (inside + outside)
            inout = diff.mean()

        return inout

    def _updateCache(self, cache, key):
        
        if not key in cache:
            cache[key] = self.stablescore(key)

    def scoreCache(self, domainlist, cache = {}):
        """
        Calculate and cache the TAD scores for any combinations of consecutive
        bottom domains.
        
        Parameters
        ----------
        domainlist : list of [start,end]
            List of bottom domain intervals in base-pair unit.
        cache : dict
            Cache of the TAD scores. The keys are in (start,end) format,
            intervals of continuous regions, and the values are corresponding
            calculated scores.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.stablescore : Calculate TAD score of any
                                                  genomic interval
            
        """
        for i in xrange(len(domainlist)):
            self._updateCache(cache, tuple(domainlist[i][:2]))
            if domainlist[i][2] > 0.2:# Noise domains will stay alone
                continue
            for j in xrange(i+1, len(domainlist)):
                if (domainlist[j][1] - domainlist[i][0] > self.maxapart):
                    break
                if domainlist[j][2] > 0.2:
                    continue
                key = (domainlist[i][0], domainlist[j][1])
                self._updateCache(cache, key)

    def maxscorepath(self, domainlist, cache):
        """
        An implementation for our proposed algorithm to find the best separation
        solution at chromosomal/domain level, given bottom domain list and
        pre-computed scores.
        
        Parameters
        ----------
        domainlist : list of [start,end]
            List of bottom domain intervals in base-pair unit.
        cache : dict
            Pre-computed scores for any combinations of continuous bottom
            domains. The keys are intervals of continuous regions in (start,end)
            format, and the values are corresponding scores.
        
        Returns
        -------
        bests : set of (sidx, eidx)
            Each element indicates one merged domain of the solution, represented
            as (start index, end index) of the input *domainlist*.
        """
        # score cache uses domain coordinate intervals
        scores = cache
        # Dynamic programming
        # paths uses domain index intervals
        paths = {}
        for i in xrange(len(domainlist)):
            Len = domainlist[i][1]//self.res - domainlist[i][0]//self.res
            paths[i] = {}
            pre = paths.get(i-1, {(-1,-1): [0,None]})
            maxp = None; maxs = float('-inf')
            for pp in pre:
                if i <= pp[1]:
                    unit = scores[(domainlist[pp[0]][0],domainlist[pp[1]][1])]
                    paths[i][pp] = [pre[pp][0]+Len*unit, pp]
                else:
                    if pre[pp][0] > maxs:
                        maxs = pre[pp][0]; maxp = pp
            for j in xrange(i, len(domainlist)):
                key = (domainlist[i][0], domainlist[j][1])
                if (key[1] - key[0] > self.maxapart) and (j > i):
                    break
                tryget = scores.get(key, None)
                if tryget is None:
                    continue
                unit = scores[(domainlist[i][0],domainlist[j][1])]
                paths[i][(i,j)] = [maxs+Len*unit, maxp]
        # Backtacking
        lastp = None; lasts = float('-inf')
        for pp in paths[len(domainlist)-1]:
            if paths[len(domainlist)-1][pp][0] > lasts:
                lasts = paths[len(domainlist)-1][pp][0]
                lastp = pp
        if lastp is None:
            return
        bestpath = [lastp]
        for i in range(len(domainlist))[::-1]:
            lastp = paths[i][lastp][1]
            if lastp == (-1,-1):
                break
            bestpath = [lastp] + bestpath

        bests = set(bestpath)

        return bests

    def maxCore(self, cache={}):
        """
        Perform TAD identification. We define TADs as domains to optimize
        chromosome separation(based on some objective function), which is
        solved by using an algorithm like dynamic programming implemented in
        :py:meth:`tadlib.hitad.chromLev.Chrom.maxscorepath`.
        
        Parameters
        ----------
        cache : dict
            TAD scores for all combinations of bottom domains will be pre-computed
            and stored in this dict. The keys are tuples (start, end)
            representing merged domain intervals. (Default: {})
        
        Attrubutes
        ----------
        maxDomains : dict
            Optimized TADs.
        
        cache : dict
            Cached TAD scores.
            
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.scoreCache : Pre-compute TAD scores for
                                                 all combinations of consecutive
                                                 bottom domains.
        tadlib.hitad.chromLev.Chrom.maxscorepath : find the best TAD list
        """
        maxDomains = {}
        for region in sorted(self.minDomains):
            rdomains = self.minDomains[region]
            if not len(rdomains):
                continue
            if len(rdomains) < 2:
                maxDomains[region] = rdomains
                continue
            self.scoreCache(rdomains, cache)
            premerge = self.maxscorepath(rdomains, cache)
            newdomain = []
            for b in sorted(premerge):
                tl = rdomains[b[0]:b[1]+1]
                newdomain.append([tl[0][0], tl[-1][1], 0, 0])
            newdomain.sort()
            if len(newdomain):
                maxDomains[region] = newdomain
                          
        self.maxDomains = self._orifilter(maxDomains)
        self.cache = cache

    def subDomains(self, domain, reflist, clv = 0, aM = None, subdomains = {}):
        """
        A recusive method (function) to identify domain hierarchy
        of a TAD (or a domain). Sub-TADs of each level are defined as the
        best separation of the outer domain.
        
        Parameters
        ----------
        domain : [start, end]
            Outer layer domain interval in base-pair unit.
        reflist : list
            List of bottom domains within the region of the outer domain.
        clv : int
            Global domain level of the outer domain. The TADs have the level
            0, the sub-TADs within a TAD have the level 1, and sub-sub-TADs
            within a sub-TAD have the level 2, and so forth. (Default: 0)
        aM : 2-D numpy ndarray or None
            Arrowhead matrix of the outer domain. If None,
            :py:meth:`tadlib.hitad.chromLev.Chrom.toArrowhead` is used to
            calculate one using *domain* interval on the fly. (Default: None)
        subdomains : dict
            A container for domains of all hierarchy. The keys are domain
            intervals in base-pair unit, and values are corresponding
            global levels.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.maxscorepath : find the best domain list
                                                   optimally separating
                                                   domain-level interactions.
        """
        def localUpdate(key):
            ori = key[0]//self.res - domain[0]//self.res
            end = key[1]//self.res - domain[0]//self.res
            # Upstream bias
            mask1 = common & (yidx >= end) & (xidx < end) & ((2*xidx - ori) >= yidx)
            # Downstream bias
            mask2 = common & (yidx < end) & (xidx >= ori) & (2*xidx >= yidx) &\
                    ((2*xidx - ori) <= yidx)
            # Update biases
            part1 = aM[mask1]
            part2 = aM[mask2]
            merge = np.r_[part1, -part2]
            if merge.size > 5:
                biases[key] = merge.mean()
            else:
                biases[key] = 0

        def getsublist(top, reflist):
            
            sublist = [d for d in reflist if ((d[0]>=top[0])and(d[1]<=top[1]))]
            sublist.sort()
            return sublist

        subdomains[(domain[0], domain[1])] = clv
        if aM is None:
            aM = self.toArrowhead(domain[0], domain[1])
        sublist = getsublist((domain[0], domain[1]), reflist)
        xidx, yidx = np.indices(aM.shape)
        common = ((yidx - xidx) > self._rm) & (aM != 1) & (aM != -1)
        biases = {}
        if len(sublist) <= 1:
            return

        for i in range(len(sublist)):
            key = (sublist[i][0], sublist[i][1])
            localUpdate(key)
            for j in range(i+1, len(sublist)):
                key = (sublist[i][0], sublist[j][1])
                localUpdate(key)
        biases[(domain[0], domain[1])] = 0

        prebests = self.maxscorepath(sublist, biases)
        #----------------------------------------------------------------------
        prebests = sorted(prebests)
        dlist = [sublist[0][0]]
        for i in range(len(prebests)-1):
            dlist.append(sublist[prebests[i][1]][1])
        dlist.append(sublist[-1][1])
        if len(dlist) == 2:
            return
        #----------------------------------------------------------------------
        for i in range(len(dlist)-1):
            nlv = clv + 1
            tmpdomain = (dlist[i], dlist[i+1])
            tmpori = tmpdomain[0]//self.res - domain[0]//self.res
            tmpend = tmpdomain[1]//self.res - domain[0]//self.res
            taM = aM[tmpori:tmpend, tmpori:tmpend]
            self.subDomains(tmpdomain, sublist, clv=nlv, aM=taM, subdomains=subdomains)

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
            return 0

        matrix = self.getSelfMatrix(domain[0], domain[1])

        total = matrix.size - np.arange(len(matrix), len(matrix)-self._rm-1, -1).sum()*2 +\
                len(matrix)

        if total < 5:
            return 0
        else:
            nx, ny = np.nonzero(matrix)
            mask = np.abs(ny-nx) > self._rm
            sigNum = mask.sum() # Number of nonzero entries apart from diag
            return 1-sigNum/total

    def toArrowhead(self, start, end):
        """
        Perform Arrowhead transformation on contact matrix of a given
        genomic region.
        
        Parameters
        ----------
        start, end : int
            The region interval.
        
        Returns
        -------
        A : 2-D numpy ndarray, float
            Transformed matrix. (A symmetric matrix)
        """
        startidx = start // self.res + self.chromLen
        endidx = end // self.res + self.chromLen

        maxd = endidx - startidx - 1

        forward = startidx - maxd
        raw = self.rawMatrix[forward:endidx, forward:endidx].toarray()

        tx, ty = np.nonzero(raw)
        raw[ty, tx] = raw[tx, ty]

        N = maxd + 1
        x, y = np.triu_indices(N)
        rawx = x + maxd
        rawy = y + maxd
        d = y - x
        minusY = rawx - d

        denominator = raw[minusY, rawx] + raw[rawx, rawy]
        denominator[denominator==0] = 1
        elemA = (raw[minusY, rawx] - raw[rawx, rawy]) / denominator

        A = np.zeros((N, N))
        A[x, y] = elemA
        A[y, x] = elemA

        return A

    def plot(self, start, end, Domains, figname = None, arrowhead = False):
        """
        Given a genomic region and a domain list, plot corresponding contact
        heatmap and all domains (represented as squares) within the region.
        
        Parameters
        ----------
        start, end : int
            The region interval.
        Domains : dict
            The keys are tuples representing gap-free regions, and values
            are corresponding identified domains.
        figname : str or None
            If not None, the figure will be saved, otherwise it will only be
            shown in an interactive window. (Default: None)
        arrowhead : bool
            If True, the Arrowhead transformed matrix will be plotted instead
            of raw contact matrix. (Default: False)
        """
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        from matplotlib.colors import Normalize

        class MidpointNormalize(Normalize):
            def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
                self.midpoint = midpoint
                Normalize.__init__(self, vmin, vmax, clip)

            def __call__(self, value, clip=None):
                x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
                return np.ma.masked_array(np.interp(value, x, y))

        def caxis_H(ax):
            """
            Axis Control for HeatMaps.
            """
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')
            ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

        def caxis_S(ax, color):
            """
            Axis Control for DI track.
            """
            for spine in ['right', 'top']:
                ax.spines[spine].set_visible(False)
            ax.tick_params(axis = 'both', bottom = False, top = False, left = False,
                           right = False, labelbottom = False, labeltop = False,
                           labelleft = False, labelright = False)
            ax.spines['left'].set_lw(1.5)
            ax.spines['left'].set_color('#1B9E77')
            ax.spines['left'].set_alpha(0.9)
            ax.spines['left'].set_linestyle('dotted')

            ax.spines['bottom'].set_lw(0.5)
            ax.spines['bottom'].set_color(color)
            ax.spines['bottom'].set_alpha(0.9)

        def properU(pos):
            """
            Express a genomic position in a proper unit (KB, MB, or both).

            """
            i_part = int(pos) // 1000000 # Integer Part
            d_part = (int(pos) % 1000000) // 1000 # Decimal Part

            if (i_part > 0) and (d_part > 0):
                return ''.join([str(i_part), 'M', str(d_part), 'K'])
            elif (i_part == 0):
                return ''.join([str(d_part), 'K'])
            else:
                return ''.join([str(i_part), 'M'])

        def boundDrawer(ax, regions, startidx, color = '#A5ACAF'):
            pairs = [(bi[0] - startidx, bi[1] - startidx) for bi in regions]
            for corner in pairs:
                if (corner[0] <= 0):
                    ax.plot([0, corner[1]], [corner[1], corner[1]], color = color,
                            linewidth = 1.5)
                    ax.plot([corner[1], corner[1]], [0, corner[1]], color = color,
                            linewidth = 1.5)
                elif (corner[1] >= interval):
                    ax.plot([corner[0], corner[0]], [corner[0], interval],
                            color = color, linewidth = 1.5)
                    ax.plot([corner[0], interval], [corner[0], corner[0]],
                            color = color, linewidth = 1.5)
                else:
                    ax.plot([corner[0], corner[0]], [corner[0], corner[1]],
                            color = color, linewidth = 1.5)
                    ax.plot([corner[0], corner[1]], [corner[0], corner[0]],
                            color = color, linewidth = 1.5)
                    ax.plot([corner[0], corner[1]], [corner[1], corner[1]],
                            color = color, linewidth = 1.5)
                    ax.plot([corner[1], corner[1]], [corner[0], corner[1]],
                            color = color, linewidth = 1.5)

        # Basic Settings
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        size = (12, 11)
        width = 0.618; Left = (1 - width) / 2
        HB = 0.1; HH = width * size[0] / size[1] # HeatMap Bottom / HeatMap Height
        SB = HB + HH # Bottom of tracks
        ST = 0.9 # The Top of tracks
        SH = (ST - SB) / 2
        raw_cmap = LinearSegmentedColormap.from_list('interaction',
                                                     ['#FFFFFF','#CD0000'])
        arrow_cmap = plt.cm.RdBu_r
        DI_color = '#1B9E77'

        fig = plt.figure(figsize = size)
        if arrowhead:
            Matrix = self.toArrowhead(start, end)
            norm = MidpointNormalize(vmin = Matrix.min(), midpoint = 0, vmax = Matrix.max())
            Params = {'norm': norm, 'cmap': arrow_cmap}
        else:
            Matrix = self.getSelfMatrix(start, end)
            nonzero = Matrix[np.nonzero(Matrix)]
            p = np.percentile(nonzero, 95)
            Params = {'vmax': p, 'cmap': raw_cmap}

        startidx = start // self.res
        endidx = end // self.res

        interval = Matrix.shape[0]
        ax = fig.add_axes([Left, HB, width, HH])
        sc = ax.imshow(Matrix, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), origin = 'lower',
                       **Params)
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ticks = list(np.linspace(0, interval, 6).astype(int))

        pos = [(start + t * self.res) for t in ticks]
        labels = [properU(i) for i in pos]

        Domains = np.r_[self.getDomainList(Domains)]
        Domains = np.r_['1,2,0', Domains[:,0], Domains[:,1], Domains[:,-1]]
        Domains[:,0] = Domains[:,0] // self.res
        Domains[:,1] = Domains[:,1] // self.res
        # Top-level
        mask = (Domains[:,2] == 0) & (Domains[:,1] > startidx) & (Domains[:,0] < endidx)
        extract = Domains[mask]
        boundDrawer(ax, extract, startidx, color = '#60636A')
        # Lower-level
        mask = (Domains[:,1] > startidx) & (Domains[:,0] < endidx) & (Domains[:,2] > 0)
        extract = Domains[mask]
        boundDrawer(ax, extract, startidx, color = '#A5ACAF')

        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        caxis_H(ax)

        axc = fig.add_axes([Left + width + 0.08, HB, 0.03, HH])
        fig.colorbar(sc, cax = axc)

        track = self.DIs[startidx:endidx]
        ax = fig.add_axes([Left, SB, width, SH])
        ax.fill_between(np.arange(track.size), track, color = DI_color,
                        edgecolor = 'none')
        ax.set_ylabel('Adaptive DI', labelpad = 43, rotation = 'horizontal',
                      style = 'italic', size = 12)
        ax.set_xlim(xlim)
        caxis_S(ax, DI_color)

        if figname != None:
            fig.savefig(figname)
            plt.close(fig)
        else:
            plt.show()

class MultiReps(DomainAligner):
    """
    We define *MultiReps* to:
    
    - Hold Hi-C data of the same chromosome from different biological replicates
      at the same time
    - Provide an interface to identify hierarchical domains using different
      replicate data independently and maintain reproducible domains of all
      levels to form final domain list.
    
    Parameters
    ----------
    chrom : str
        Chromosome label.
    res : int
        Hi-C data resolution in base-pair unit.
    datasets : dict
        The keys are unique replicate labels, and the values are constructed
        *Chrom* objects using Hi-C data of corresponding biological replicates.
    """
    def __init__(self, chrom, res, datasets):

        self.chrom = chrom
        self.res = res
        self.Queue = datasets # key: rep label, value: Chrom object
        self._state = 'Unfinished'

    def getDomainList(self, byregion):
        """
        Combine by-region domains into a single list.
        
        Each domain in the returned list is represented in the format
        ``[chromosome label, start bp, end bp, hierarchical level]``.
        """
        domainlist = []
        for region in sorted(byregion):
            tmp = byregion[region]
            for d in tmp:
                if d[1] - d[0] >= 3*self.res:
                    domainlist.append([self.chrom]+d[:2]+[d[-1]])
        return domainlist

    def callDomain(self):
        """
        Invoke *callDomain* method for each loaded :py:class:`tadlib.hitad.chromLev.Chrom`
        object, and find reproducible domains between replicates by using our
        domain-level alignment strategy.
        """
        for rep in self.Queue:
            if self.Queue[rep]._state != 'Completed':
                self.Queue[rep].callDomain()
        
        reps = sorted(self.Queue)
        tl = self.getDomainList(self.Queue[reps[0]].hierDomains)
        tg = DomainSet(reps[0], tl, self.res)
        for rep in reps[1:]:
            ql = self.getDomainList(self.Queue[rep].hierDomains)
            qy = DomainSet(rep, ql, self.res)
            tl = self._interp(self.reproducible(tg, qy))
            tg = DomainSet('mtg', tl, self.res)
        
        self.mergedDomains = tl
        
        self._state = 'Completed'

    def reproducible(self, tg, qy):
        """
        Subprocess of *callDomain* for replicate alignment and reproducibility
        determination.
        """
        work = repAligner(tg, qy)
        tcache = work._align(tg, qy)
        qcache = work._align(qy, tg)
        source = {}
        tbase = self._customize(tcache, qcache, tg, qy, source)
        qbase = self._customize(qcache, tcache, qy, tg, source)
        pool = set()
        self._extract(tbase, pool, tg)
        self._extract(qbase, pool, qy)
        dlist = hierFormat(pool)
        
        return dlist
    
    def _iscompatible(self, dlist, ref):
        
        dlist.sort()
        lvs = []
        for d in dlist:
            lvs.append(ref.boundclass[(d[0], d[1])])
            lvs.append(ref.boundclass[(d[0], d[2])])
        if len(lvs) <= 2:
            return True
        else:
            inlv = min(lvs[1:-1])
            blv = max(lvs[0], lvs[-1])
            return inlv >= blv
    
    def _customize(self, cache, ref, tg, qy, source):
        
        tcore = {}
        for k in sorted(cache):
            if k[1] is None:
                continue
            sk = k[::-1]
            if not sk in ref:
                continue
            if sk in source:
                continue
            hit = ref[sk]
            target = cache[k]
            tl, ql = target.info
            if (len(tl)==1) and (len(ql)>1):
                # wait for reverse alignment
                continue
            tcore[k] = Container(target.info)
            for lv in target:
                for p in target[lv]:
                    tl, ql = target[lv][p].info
                    tlvs = set([d[-1] for d in tl])
                    tcheck = self._iscompatible(tl, tg)
                    qcheck = self._iscompatible(ql, qy)
                    found = self._search(hit, p)
                    if found and ((len(tl)==1) or (len(ql)==1)) and tcheck and qcheck:
                        mlv = min(tlvs)
                        if (mlv==0) and (len(tl)>1):
                            for d in tl:
                                nk = (tuple(d[:-1]),None)
                                if not d[-1] in tcore[k]:
                                    tcore[k][d[-1]] = {nk:Container([[d],[]])}
                                else:
                                    tcore[k][d[-1]][nk] = Container([[d],[]])
                        else:
                            if not mlv in tcore[k]:
                                tcore[k][mlv] = {p:target[lv][p]}
                            else:
                                tcore[k][mlv][p] = target[lv][p]
                    else:
                        for d in tl:
                            nk = (tuple(d[:-1]),None)
                            if not d[-1] in tcore[k]:
                                tcore[k][d[-1]] = {nk:Container([[d],[]])}
                            else:
                                tcore[k][d[-1]][nk] = Container([[d],[]])
        
        source.update(tcore)
        
        return tcore

    def _extract(self, base, pool, tg):
        
        for k in sorted(base):
            tl, ql = base[k].info
            topobounds = {}  # bounds of higher levels should be passed down
            if len(tl)==len(ql)==1:
                d = tuple(tl[0])
                topobounds[d] = tg.Domains[d]
                new = self.intersect(k[0],k[1])
                pool.add(tuple(new))
                self._updateTopo(topobounds, d, new)
            else:
                for d in map(tuple, tl):
                    topobounds[d] = tg.Domains[d]
                    self._updateTopo(topobounds, d, [])

            for lv in base[k]:
                for t in base[k][lv]:
                    tl, ql = base[k][lv][t].info
                    if (len(tl)==1) and (len(ql)==0):
                        new = tl[0][:-1]
                    if (len(tl)==1) and (len(ql)>=1):
                        new = self.intersect(t[0],t[1])
                    if (len(tl)>1) and (len(ql)==1):
                        new = self.intersect(t[0],t[1])
                    bounds = self._getTopo(topobounds, tuple(tl[0]))
                    assert not bounds is None
                    if len(bounds):
                        if len(ql) == 0:
                            continue
                        new = self.intersect(new, bounds)
                    if new[2] - new[1] < 3*self.res:
                        continue
                    if not len(bounds):
                        if (len(tl)==1) and (len(ql)<=1):
                            pool.add(tuple(new))
                            self._updateTopo(topobounds, tuple(tl[0]), new)
                    else:
                        pool.add(tuple(new))
                        for d in tl:
                            self._updateTopo(topobounds, tuple(d), new)
    
    def _updateTopo(self, topo, d, bounds):
        
        for k in topo:
            if not d is None:
                if k == d:
                    topo[k].bounds = bounds
                    d = self._updateTopo(topo[k], None, bounds)
                else:
                    d = self._updateTopo(topo[k], d, bounds)
                if d is None:
                    break
            else:
                topo[k].bounds = bounds
                self._updateTopo(topo[k], None, bounds)
        
        return d
    
    def _getTopo(self, topo, d):
        
        bounds = None
        for k in topo:
            if k == d:
                bounds = topo[k].bounds
            else:
                bounds = self._getTopo(topo[k], d)
            if not bounds is None:
                break
        
        return bounds
    
    def intersect(self, d1, d2):
        
        return [d1[0], max(d1[1],d2[1]), min(d1[2],d2[2])]

    def _interp(self, rawlist):
        """
        Fill holes for domains of all levels. We skip gap regions detected
        on any biological replicate dataset.
        
        Consider [['1',0,200000],['1',350000,600000]], then the region
        [200000,350000] is defined as a hole in original domain list.
        
        Parameters
        ----------
        rawlist : list
            Original domain list. Each domain is represented by
            [chrom,start,end,label], in which *start* and *end* should be
            in base-pair unit, and *label* indicates the hierarchical level
            of the domain.
            
        Returns
        -------
        domainlist : list
            Domain list after hole filling.
        
        See Also
        --------
        tadlib.hitad.chromLev.Chrom.splitChrom : where gap regions are detected
        """
        self.gapbins = set()
        for rep in self.Queue:
            self.gapbins.update(self.Queue[rep].gapbins)
        
        bo = BoundSet('pse', rawlist, self.res)
        D = set()
        newlist = []
        for d in rawlist:
            D.add(tuple(d[:3])) # Hash the domain list
            newlist.append(d[:3])

        for i in xrange(len(bo.Bounds)-1):
            if bo.Bounds[i+1][0] != bo.Bounds[i][0]:
                continue
            td = (bo.Bounds[i][0], bo.Bounds[i][1], bo.Bounds[i+1][1])
            if td[2] - td[1] < 3*self.res:
                continue
            dbset = set(range(td[1],td[2],self.res))
            gapratio = len(dbset.intersection(self.gapbins)) / len(dbset)
            if gapratio >= 0.2:
                continue
            if not td in D:
                newlist.append(list(td))

        domainlist = hierFormat(newlist)

        return domainlist

class repAligner(DomainAligner):
    """
    A customized *tadlib.hitad.aligner.DomainAligner* for hierarchical
    domain alignment between two replicates.
    
    The API stays the same.
    """
    def __init__(self, *args):
        DomainAligner.__init__(self, *args)
    
    def _align(self, tg, qy):
        
        # First, align the top-level domains
        ttl = [d for d in tg.Domains if d[-1]==0]
        tql = [d for d in qy.Domains if d[-1]==0]
        ttg = DomainSet(tg.Label, ttl, tg.res, False)
        tqy = DomainSet(qy.Label, tql, qy.res, False)
        pairs, _ = self._aligncore(ttg, tqy)
        
        cache = {}
        for k in sorted(pairs):
            tl = tg.getregion(*k[0])
            ftg = DomainSet(tg.Label, tl, tg.res)
            ql = qy.getregion(*k[1])
            fqy = DomainSet(qy.Label, ql, qy.res)
            cache[k] = Container(pairs[k])
            # Parse the hierarchy step by step
            ori = 1 if len(pairs[k][0])==1 else 0
            for tv in range(ori, max(ftg.levs)+1):
                cache[k][tv] = {}
                tl = ftg.getregion(*(k[0]+(tv,)))
                ntg = DomainSet(ftg.Label, tl, ftg.res, False)
                ql = self._localhits(ntg, fqy)
                nqy = DomainSet(fqy.Label, ql, fqy.res, False)
                npairs, _ = self._aligncore(ntg, nqy)
                for p in npairs:
                    cache[k][tv][p] = Container(npairs[p])
        
        return cache
    
    def _toberobust(self, tl, ql, tg, qy, t_ref, q_ref, tol=0.9):
        
        if t_ref is None:
            tl = tg.getregion(tl[0][0], tl[0][1], tl[-1][2])
            ql = qy.getregion(ql[0][0], ql[0][1], ql[-1][2])
        else:
            tl = t_ref.getregion(tl[0][0], tl[0][1], tl[-1][2], 0)
            ql = q_ref.getregion(ql[0][0], ql[0][1], ql[-1][2], 0)
        
        tl.sort()
        ql.sort()
        bylev = [[]]
        olev = ql[0][-1]
        for d in ql:
            if d[-1] == olev:
                bylev[-1].append(d)
            else:
                olev = d[-1]
                bylev.append([d])
        if len(bylev) > 1:
            # query domain levels are heterogeneous
            ref = self.overlap([tl[0][1],tl[-1][2]], [ql[0][1],ql[-1][2]])
            score = 0
            hit = []
            for c in bylev:
                tmp = self.overlap([tl[0][1],tl[-1][2]], [c[0][1],c[-1][2]])
                if tmp / ref > score:
                    score = tmp / ref
                    hit = c
            if score >= tol:
                ql = hit
                
        tk = (tl[0][0], tl[0][1], tl[-1][2])
        qk = (ql[0][0], ql[0][1], ql[-1][2])
        
        return tk, qk, [tl, ql]