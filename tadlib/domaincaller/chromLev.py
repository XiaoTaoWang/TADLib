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