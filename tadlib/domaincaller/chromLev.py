# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:47:34 2016

@author: wxt
"""

from __future__ import division
import numpy as np
from tadlib.hitad import chromLev
import copy

chromCore = chromLev.Chrom

np.seterr(divide = "ignore")

class Chrom(chromCore):

    def __init__(self, chrom, res, hicdata, window=2000000, minsize=5):

        self.chrom = chrom
        self.res = res
        self._rm = 1
        self.chromLen = hicdata.shape[0]
        self.hmm = None
        self._dw = window // res
        self.minsize = minsize

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

        # fixed windows for DI
        self.windows = np.ones(self.chromLen, dtype=np.int32) * self._dw
    
    def callDomains(self):

        self.calDI(self.windows, 0)
        self.splitChrom(self.DIs)
        self.domains = self.getDomainList(self.minCore(self.regionDIs))
    
    def DIsplit(self, arr):

        intervals = []
        # runs of negative values
        tmp = set(np.where(arr>=0)[0])
        breaks = copy.deepcopy(tmp)
        for b in tmp:
            breaks.add(b+1)
        ori = np.split(arr, sorted(breaks))
        start = 0
        for a in ori:
            if (len(a) > 1) and (a[0] < 0):
                rb = start+len(a)
                si = np.where(a[:-1]*5<a[1:])[0] # default factor: 5
                if si.size > 0:
                    lb = start + si[0]
                    if rb - lb > 1:
                        intervals.append((lb, rb))
            start += len(a)
        
        # runs of positive values
        tmp = set(np.where(arr<=0)[0])
        breaks = copy.deepcopy(tmp)
        for b in tmp:
            breaks.add(b+1)
        ori = np.split(arr, sorted(breaks))
        start = 0
        for a in ori:
            if (len(a) > 1) and (a[0] > 0):
                lb = start
                ei = np.where(a[1:]*5>a[:-1])[0]
                if ei.size > 0:
                    rb = start + ei[-1] + 2
                    if rb - lb > 1:
                        intervals.append((lb, rb))
            start += len(a)
        
        intervals.sort()
        
        return intervals
    
    def refine_bounds(self, seq, b):

        subseq = seq[b[0]:b[1]]
        parse = self.DIsplit(subseq)
        starts = ()
        min_DI = 0
        ends = ()
        max_DI = 0
        for i in parse:
            tmp = subseq[i[0]:i[1]].min()
            if tmp < min_DI:
                min_DI = tmp
                starts = i
            tmp = subseq[i[0]:i[1]].max()
            if tmp > max_DI:
                max_DI = tmp
                ends = i
        
        if len(starts) and len(ends):
            if starts[1] > ends[0]:
                return []
        
        left = b[0]
        if len(starts):
            if starts == parse[0]:
                left = b[0] + starts[0]
            else:
                tmp = subseq[:starts[0]].min()
                if tmp >= 0:
                    left = b[0] + starts[0]
                else:
                    fold = min_DI / tmp
                    if fold > 2:
                        left = b[0] + starts[0]
        
        right = b[1]
        if len(ends):
            if ends == parse[-1]:
                right = b[0] + ends[1]
            else:
                tmp = subseq[ends[1]:].max()
                if tmp <= 0:
                    right = b[0] + ends[1]
                else:
                    fold = max_DI / tmp
                    if fold > 2:
                        right = b[0] + ends[1]
        
        return [left, right]
        
    def pipe(self, seq, start):

        # bin-level domain (not base-pair-level domain!)
        bounds = self._getBounds(self.viterbi(seq), junctions=['30'])
        pairs = [[bounds[i], bounds[i+1]] for i in range(len(bounds)-1)]
        domains = []
        for b in pairs:
            rb = self.refine_bounds(seq, b)
            if len(rb):
                # start, end, noise level, hierarchical level
                tmp = [rb[0]+start, rb[1]+start, 0, 0]
                domains.append(tmp)

        return domains
    
    def viterbi(self, seq):
        
        path = [int(s.name) for i, s in self.hmm.viterbi(seq)[1][1:-1]]

        return path