# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:05:35 2016

@author: wxt

"""
from __future__ import division
import bisect
from tadlib.hitad.graph import Graph

#------------------------------------------------------------------------------
# Boundary-Based Alignment
class SingleBound(object):
    """
    *SingleBound* is defined to:
    
    - Represent a single bound (chrom, pos)
    - Map the bound to a pool of bounds
    
    Parameters
    ----------
    chrom : str
        Chromosome label.
    
    pos : int
        Bound position on the chromosome.
    
    Attributes
    ----------
    chrom : str
        Chromosome label.
    
    pos : int
        Position on the chromosome.
    
    cache : dict
        Container for matched details.
    
    """
    def __init__(self, chrom, pos):
        self.chrom = chrom
        self.pos = pos
        self.cache = {}
    
    def align(self, qn, qb, tol):
        """
        Map the bound to *qb* using the binary search method.
        
        Parameters
        ----------
        qn : str
            Unique identifier for *qb*.
        
        qb : list of tuples
            Reference bound (remember reference genome?) list. Each element
            is a tuple (chrom, pos) representing a single bound. And the list
            must be sorted in advance for binary search.

        tol : int
            Mismatch tolerance. If the genomic distance between the bound
            and the best hit is less than this value, we say we have found
            a match, otherwise the bound is missed in *qb*.

        Notes
        -----
        Internally, the matched details will be stored in *cache* under the
        key *qn*, the value is also a dict recording the matched bound
        index (midx) and the indices of the matched neighbors (nindices).
        
        """
        self.cache[qn] = {}
        tidx = bisect.bisect(qb, (self.chrom, self.pos))
        lidx = ridx = -1
        assert len(qb)>0, 'Empty query domain list'
        if (self.chrom, self.pos) == qb[tidx-1]:
            lidx = ridx = tidx - 1
        else:
            if tidx <= len(qb) - 1:
                if self.chrom != qb[tidx-1][0]:
                    if self.chrom == qb[tidx][0]:
                        lidx = ridx = tidx
                else:
                    lidx = tidx - 1
                    if self.chrom == qb[tidx][0]:
                        ridx = tidx
                    else:
                        ridx = lidx
            else:
                if self.chrom == qb[tidx-1][0]:
                    lidx = ridx = tidx - 1
        if tidx == 0:
            if self.chrom == qb[tidx][0]:
                lidx = ridx = tidx
        
        self.cache[qn]['midx'] = None
        self.cache[qn]['nindices'] = None
        if lidx != -1:
            midx = lidx
            shift = abs(self.pos - qb[midx][1])
            rshift = abs(self.pos - qb[ridx][1])
            if rshift < shift:
                midx = ridx
            # Matched bound infomation
            self.cache[qn]['midx'] = midx
            nindices = set([lidx, ridx])
            tidx = lidx - 1
            while tidx >= 0:
                if qb[tidx][0]!=self.chrom:
                    break
                nindices.add(tidx)
                tshift = abs(self.pos - qb[tidx][1])
                if tshift>=tol:
                    break
                tidx -= 1
            tidx = ridx + 1
            while tidx < len(qb):
                if qb[tidx][0]!=self.chrom:
                    break
                nindices.add(tidx)
                tshift = abs(self.pos - qb[tidx][1])
                if tshift>=tol:
                    break
                tidx += 1
            self.cache[qn]['nindices'] = sorted(nindices)

class BoundSet(object):
    """
    As the name suggests, we use *BoundSet* to hold all bounds of a domain
    list.
    
    Parameters
    ----------
    en : str
        Unique identifier for current set of bounds.
    
    domainlist : list
        List of the domains. Each domain is represented by
        ``[chrom,start,end,level]``. I think *chrom*, *start* and *end* are
        self-explanatory, all you need to keep in mind is that *start and
        *end* should be in base-pair unit. *level* indicates the hierarchical
        level of the domain. In our work, TAD is denoted as 0, sub-TAD is
        denoted as 1, and subsequent domain level is denoted as 2, etc.
    
    res : int
        Resolution of the Hi-C data in base-pair unit.
    
    Attributes
    ----------
    Label : str
        Unique identifier.
    
    boundclass : dict
        The keys are bound representations (chrom,pos), and the values indicate
        corresponding hierarchical level notations. The Level of a bound is
        determined by the domain with the lowest level notation. For example,
        if we have two domains, ['1',100000,500000,0] and ['1',100000,200000,1],
        according to our definition, the level of ('1',200000) is 1, but the
        level of ('1',100000) is 0.
    
    Bounds : list
        Sorted bound list. This attribute can be used as the reference
        bound list in :py:meth:`tadlib.hitad.aligner.SingleBound.align`
        directly.
    
    """
    
    def __init__(self, en, domainlist, res):
        self.Label = en
        self.boundclass = {}
        for domain in domainlist:
            chrom, start, end, label = domain
            if end - start < 5*res:
                continue
            if (chrom, start) in self.boundclass:
                if label < self.boundclass[(chrom, start)]:
                    self.boundclass[(chrom, start)] = label
            else:
                self.boundclass[(chrom, start)] = label
            if (chrom, end) in self.boundclass:
                if label < self.boundclass[(chrom, end)]:
                    self.boundclass[(chrom, end)] = label
            else:
                self.boundclass[(chrom, end)] = label 
        
        self.Bounds = sorted(self.boundclass) # bounds of all hierarchy
        self._boundIdx = dict(zip(self.Bounds, range(len(self.Bounds))))

#------------------------------------------------------------------------------
#  A kind of domain-based hierarchical alignment
class DomainSet(BoundSet):
    """
    Parse and hold a hierarchical domain set.
    
    Parameters
    ----------
    en : str
        Unique identifier for input domain set.
    domainlist : list
        List of domains. See :py:class:`tadlib.hitad.aligner.BoundSet` for
        details.
    res : int
        Resolution of the Hi-C data in base-pair unit.
    hier : bool
        Whether *domainlist* contains multiple-level domains or not.
        (Default: True)
    
    Attributes
    ----------
    res : int
        Resolution of the Hi-C data in base-pair unit.
    levs : set of int
        All possible domain levels contained in *domainlist*.
    bychroms : dict
        Bychromosomal rearrangement of *domainlist*. The keys are chromosome
        labels(1,2,...,22,X,Y), and the values are list of [start,end,level].
    pretree : dict
        Nested domain list within any domain interval. Returned by
        :py:meth:`tadlib.hitad.aligner.DomainSet.NestedDomains`.
    subpool : dict
        Domain list within any domain interval.
        (:py:meth:`tadlib.hitad.aligner.DomainSet.NestedDomains`)
    lidx : dict
        The smallest indices of left domain boundaries in a by-chromosomal
        domain list. (:py:meth:`tadlib.hitad.aligner.DomainSet.NestedDomains`)
    ridx : dict
        The largest indices of right domain boundaries in a by-chromosomal
        domain list. (:py:meth:`tadlib.hitad.aligner.DomainSet.NestedDomains`)
    Domains : dict
        Store each TAD and its nested domains as a tree. Each node in the
        tree indicates one domain, in particular, the root corresponds to
        the TAD, and the leaves correspond to bottom domains.
    
    """
    def __init__(self, en, domainlist, res, hier = True):
        
        self.res = res
        BoundSet.__init__(self, en, domainlist, res)
        bychroms = {}
        self.levs = set()
        for domain in domainlist:
            chrom, start, end, label = domain
            if end - start < 5*res:
                continue
            if not chrom in bychroms:
                bychroms[chrom] = [[start, end, label]]
            else:
                bychroms[chrom].append([start, end, label])
            self.levs.add(label)
        self.bychroms = bychroms
        self.pretree, self.subpool, self.lidx, self.ridx = self.NestedDomains(bychroms)
        self.Domains = {}
        for d in self.pretree:
            if hier:
                if d[-1] > 0:
                    continue
            hit = self.pretree[d]
            self.Domains[d] = Node()
            self.genDomainTree(self.Domains[d], self.pretree, hit)
    
    def getBottoms(self):
        """
        Link bottom domains to corresponding outer TADs.
        
        Attributes
        ----------
        bottoms : dict
            Used to quickly retrieve (bottom domain, TAD) pairs.
            
        """
        self.bottoms = {}
        for d in self.subpool:
            if d[-1] == 0:
                for sub in self.subpool[d]:
                    if not len(self.pretree[tuple(sub)]):
                        self.bottoms[tuple(sub[:-1])] = list(d)
    
    def NestedDomains(self, bychroms):
        """
        Pre-parse domain lists for accelerating subsequent calculations.
        
        Parameters
        ----------
        bychroms : dict
            By-chromosomal domain lists.
        
        Returns
        -------
        tmpdict : dict
            Nested domain list within any domain interval. If a domain have
            no nested domains, then its value is an empty list.
        subpool : dict
            Domain list within any domain interval. Different from *tmpdict*,
            if a domain have no nested domains, the value is a list only
            containing itself.
        lidx : dict
            The smallest indices of left domain boundaries in a by-chromosomal
            domain list.
        ridx : dict
            The largest indices of right domain boundaries in a by-chromosomal
            domain list.
        """
        tmpdict = {} # Pre-tree
        subpool = {}
        lidx = {}; ridx = {}
        for c in bychroms:
            sidx = 0
            bychroms[c].sort()
            for q in bychroms[c]:
                pres = []; pool = []
                label = 1
                key = (c,)+tuple(q)
                lk = (c, q[0])
                rk = (c, q[1])
                lidx[lk] = min(lidx.get(lk, len(bychroms[c])-1), sidx)
                lidx[rk] = min(lidx.get(rk, len(bychroms[c])-1), sidx)
                for i in range(sidx, len(bychroms[c])):
                    if (bychroms[c][i][0]>=q[0]) and (bychroms[c][i][1]<=q[1]):
                        if bychroms[c][i][-1]==q[-1]+1:
                            pres.append([c]+bychroms[c][i])
                        pool.append([c]+bychroms[c][i])
                    if (bychroms[c][i][0]>=q[0]) and label:
                        sidx = i
                        label = 0
                    if bychroms[c][i][0]>q[1]:
                        break
                ridx[lk] = max(ridx.get(lk, 0), i)
                ridx[rk] = max(ridx.get(rk, 0), i)
                tmpdict[key] = pres
                subpool[key] = pool
        
        return tmpdict, subpool, lidx, ridx
    
    def genDomainTree(self, node, pretree, cur):
        """
        Recursively generate a tree/sub-tree taking *node* as starting point.
        
        Parameters
        ----------
        node : Node
            A dict-like container for current domain.
        pretree : dict
            Nested domain list within any domain interval.
            (:py:meth:`tadlib.hitad.aligner.DomainSet.NestedDomains`)
        cur : list
            Nested domain list within current domain interval.
            
        """
        for d in cur:
            node[tuple(d)] = Node(list(d[:-1]))
            hit = pretree[tuple(d)]
            self.genDomainTree(node[tuple(d)], pretree, hit)
    
    def getregion(self, chrom, start, end, lev=None):
        """
        Extract all domains (or domains at specific level) within a given
        region.
        
        Parameters
        ----------
        chrom : str
            Chromosome label.
        start, end : int
            Domain interval in base-pair unit.
        lev : int or None
            Specify the desired domain level. (Default: None, domains of all
            levels will be returned)
        
        Returns
        -------
        rdomains : list
            Sorted domain list. Each element corresponds to one domain in the
            format ``[chrom,start,end,level]``.
        """
        
        rdomains = []
        sidx = self.lidx[(chrom, start)]
        eidx = self.ridx[(chrom, end)]+1
        candis = self.bychroms[chrom][sidx:eidx]
        cache = {}
        for d in candis:
            tmp = [chrom] + d
            if d[0] < start:
                continue
            if d[1] > end:
                continue
            for sub in self.subpool[tuple(tmp)]:
                if not len(self.pretree[tuple(sub)]):
                    if tuple(sub) in cache:
                        cache[tuple(sub)][tmp[-1]] = tmp
                    else:
                        cache[tuple(sub)] = {tmp[-1]:tmp}
        pool = set()
        for key in cache:
            for l in cache[key]:
                if not lev is None:
                    if (l != lev) and ((max(cache[key]) >= lev) or (l < max(cache[key]))):
                        continue
                if not tuple(cache[key][l]) in pool:
                    rdomains.append(cache[key][l])
                    pool.add(tuple(cache[key][l]))
        rdomains.sort()
            
        return rdomains
        
class SingleDomain(object):
    """
    We use *SingleDomain* to:
    
    1. Represent a single domain (chrom, start, end).
    2. Map the domain to another domain set
    
    Parameters
    ----------
    chrom : str
        Chromosome label.
    start, end : int
        Interval of the domain in base-pair unit.
    
    Attributes
    ----------
    chrom : str
        Chromosome label.
    interval : list
        [start, end]
    cache : dict
        Container for matched details.
        
    """
    def __init__(self, chrom, start, end):
        
        self.chrom = chrom
        self.interval = [start, end]
        self.cache = {}
    
    def overlap(self, ti, qi):
        """
        Calculate overlap ratio of any two regions.
        
        Parameters
        ----------
        ti, qi : list
            Interval ([start,end]) of the region.
        
        Returns
        -------
        OR : float, 0-1
            Overlap ratio.
        """
        if (ti[1]<=qi[0]) or (qi[1]<=ti[0]):
            return 0
        mi = ti + qi
        mi.sort()
        OR = (mi[2]-mi[1])/(mi[3]-mi[0]) # intersect / union
        return OR
    
    def align(self, qy):
        """
        Find the domain *D* in *qy* maximizing the overlap ratio. Binary
        search method is used internally for accelerating the search process.
        
        Parameters
        ----------
        qy : a :py:class:`tadlib.hitad.aligner.DomainSet` instance
            Reference domain set. (Recall sequence mapping and reference
            genome)
        
        Notes
        -----
        The matched details are stored in *cache* using the unique identifier
        of *qy* (``qy.Label``) as the key, the value is also a dict with 2
        keys: *hitdomain* records the matched domain interval in
        (chrom,start,end) format, and *hitoverlap* records the overlap ratio
        between the hitdomain and current query domain.
        """
        qn = qy.Label
        self.cache[qn] = {}
        qb = qy.Bounds
        res = qy.res
        self.cache[qn]['hitdomain'] = None
        self.cache[qn]['hitoverlap'] = 0

        left = SingleBound(self.chrom, self.interval[0])
        right = SingleBound(self.chrom, self.interval[1])
        left.align(qn, qb, 5*res)
        if left.cache[qn]['nindices'] is None:
            return
        right.align(qn, qb, 5*res)
        if right.cache[qn]['nindices'] is None:
            return
            
        lidx = left.cache[qn]['nindices'][0]
        ridx = right.cache[qn]['nindices'][-1]
        candis = []
        for i in range(lidx, ridx):
            for j in range(i+1, ridx+1):
                tmp = [qb[t][1] for t in range(i,j+1)]
                score = self.overlap(self.interval, [tmp[0],tmp[-1]])
                candis.append((score, tmp))
        
        if not len(candis):
            return
        besthit = max(candis)
        if besthit[0] == 0:
            return
        hitdomain = (self.chrom, besthit[1][0], besthit[1][-1])
        self.cache[qn]['hitdomain'] = list(hitdomain)
        self.cache[qn]['hitoverlap'] = besthit[0]

class Container(dict):
    """
    Dict-like. Used in the organizing of domain alignment results.
    
    Parameters
    ----------
    info : list
        Pair of domain lists from two domain sets.
    
    """
    def __init__(self, info):
        self.info = info

class Node(dict):
    """
    Dick-like. We use it to represent nodes of a hierarchical domain tree
    in *DomainSet*.
    
    Parameters
    ----------
    bounds : list or None
        Domain interval represented by [chrom,start,end].
    
    """
    def __init__(self, bounds=None):
        self.bounds = bounds

class DomainAligner(SingleDomain):
    """
    This class is the work horse we define to:
        
    1. Hold multiple :py:class:`tadlib.hitad.aligner.DomainSet` instances
       at the same time.
    2. Perform domain-based hierarchical alignment between any two
       :py:class:`tadlib.hitad.aligner.DomainSet`.
    3. Define and extract domain-level change types from alignment results
       between two :py:class:`tadlib.hitad.aligner.DomainSet`.
    
    Parameters
    ----------
    args : two or more :py:class:`tadlib.hitad.aligner.DomainSet` instances
    
    Attributes
    ----------
    DomainSets : dict
        Pool of :py:class:`tadlib.hitad.aligner.DomainSet`. The keys are
        unique identifiers extracted from :py:class:`tadlib.hitad.aligner.DomainSet`,
        and the values are corresponding :py:class:`tadlib.hitad.aligner.DomainSet`
        instances.
    Results : dict
        Container for alignment results between any pair of
        :py:class:`tadlib.hitad.aligner.DomainSet` instances.
    """
    def __init__(self, *args):
        self.DomainSets = {}
        for domains in args:
            self.DomainSets[domains.Label] = domains
        self.Results = {}
    
    def _getTree(self, start, end, ref, pool):
        
        # all domains contained in *ref* must lie on the same chromosome
        for d in ref: # four columns
            if d[1] >= end:
                continue
            if d[2] <= start:
                continue
            if (d[1]>=start) and (d[2]<=end):
                pool[tuple(d[:3])] = d[-1]
                continue
            self._getTree(start, end, ref[d], pool)
    
    def _oneway(self, tg, qy, t_ref, q_ref, vs, es):
        for td in sorted(tg.Domains):
            if t_ref is None:
                tk = (tg.Label,) + td
            else:
                tk = (tg.Label,) + tuple(t_ref.bottoms[td[:-1]])
            sd = SingleDomain(*td[:3])
            sd.align(qy)
            pse = sd.cache[qy.Label]['hitdomain']
            vs.add(tk)
            if pse is None:
                continue
            args = tuple(pse)
            qds = qy.getregion(*args)
            if not len(qds):
                continue
            if q_ref is None:
                qks = [(qy.Label,)+tuple(d) for d in qds]
            else:
                qks = [(qy.Label,)+tuple(q_ref.bottoms[tuple(d[:-1])]) for d in qds]
            for k in qks:
                vs.add(k)
                op = self.overlap([td[1],td[2]], [k[2],k[3]])
                es[(tk,k)] = op
    
    def _localhits(self, tg, qy):
        
        table = {}
        for td in tg.Domains:
            sd = SingleDomain(*td[:3])
            sd.align(qy)
            pse = sd.cache[qy.Label]['hitdomain']
            if pse is None:
                continue
            self._getTree(pse[1], pse[2], qy.Domains, table)
        
        if not len(table):
            return []
            
        rmin = min([d[1] for d in table])
        rmax = max([d[2] for d in table])
        reorg = hierFormat(table.keys())
        nqy = DomainSet(qy.Label, reorg, qy.res)
        pool = {}
        self._getTree(rmin, rmax, nqy.Domains, pool)
        ql = []
        for d in pool:
            ql.append(list(d)+[table[d]])
            
        return ql
    
    def _aligncore(self, tg, qy, t_ref=None, q_ref=None):
        
        # *tg* and **qy only contain single-level domains
        vs = set()
        if len(qy.Domains):
            pes = {} # The initial edges are directed
            self._oneway(tg, qy, t_ref, q_ref, vs, pes)
            self._oneway(qy, tg, q_ref, t_ref, vs, pes)
            es = {} # Only symmetric edges are retained
            for e in pes:
                if e[::-1] in pes:
                    es[e] = pes[e]
        else:
            es = {}
            for d in tg.Domains:
                if t_ref is None:
                    vs.add((tg.Label,)+d)
                else:
                    vs.add((tg.Label,)+tuple(t_ref.bottoms[d[:-1]]))
        g = Graph(vs, es)
        pairs = {}
        self._pairfromgraph(g, tg, qy, t_ref, q_ref, pairs)
            
        return pairs, g
    
    def _pairfromgraph(self, g, tg, qy, t_ref, q_ref, pairs):
        
        tmp = g.find_connected_components()
        for c in tmp:
            med = [[],[]]
            for d in c:
                if d[0]==tg.Label:
                    med[0].append(list(d[1:]))
                else:
                    med[1].append(list(d[1:]))
            med[0].sort()
            med[1].sort()
            if (len(med[0])==0) or (len(med[1])==0):
                continue
            tk, qk, med = self._toberobust(med[0], med[1], tg, qy, t_ref, q_ref)
            pairs[(tk,qk)] = med
        
        return pairs
    
    def _toberobust(self, tl, ql, tg, qy, t_ref, q_ref):
        
        if t_ref is None:
            tl = tg.getregion(tl[0][0], tl[0][1], tl[-1][2])
            ql = qy.getregion(ql[0][0], ql[0][1], ql[-1][2])
        else:
            tl = t_ref.getregion(tl[0][0], tl[0][1], tl[-1][2], 0)
            ql = q_ref.getregion(ql[0][0], ql[0][1], ql[-1][2], 0)
        
        tl.sort()
        ql.sort()
                
        tk = (tl[0][0], tl[0][1], tl[-1][2])
        qk = (ql[0][0], ql[0][1], ql[-1][2])
        
        return tk, qk, [tl, ql]
    
    def align(self, tn, qn):
        """
        Construct hierarchical alignment between *tn* and *qn*.
        
        Parameters
        ----------
        tn, qn : str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances which are collected by *arg* during initialization.
        
        Notes
        -----
        The alignment results are organized in a hierarchical way in a
        dictionary. The keys are matched chromosome region (corresponds to
        either one TAD or several continuous TADs) pairs at the TAD level,
        and the values are :py:class:`tadlib.hitad.aligner.Container` instances
        with *info* attribute set to be pair of detailed TAD lists; the keys
        of these :py:class:`tadlib.hitad.aligner.Container` indicate domain
        levels, the values again are dictionaries containing sub-alignment
        results within the upper-layer TAD region.
        
        You can access the results from *Results* attribute:
        ``self.Results[tn][qn]`` or ``self.Results[qn][tn]``.
        """
        tg = self.DomainSets[tn]
        qy = self.DomainSets[qn]
        # Originally, *tg* and *qy* are always hierarchical
        tcache = self._align(tg, qy)
        qcache = self._align(qy, tg)
        tcore = self._crosscorrect(tcache, qcache)
        qcore = self._crosscorrect(qcache, tcache)
        
        if not tn in self.Results:
            self.Results[tn] = {qn:tcore}
        else:
            self.Results[tn][qn] = tcore

        if not qn in self.Results:
            self.Results[qn] = {tn:qcore}
        else:
            self.Results[qn][tn] = qcore   
    
    def _align(self, tg, qy):
        """
        Match domains of *tg* with domains in *qy*.
        
        Parameters
        ----------
        tg, qy : :py:class:`tadlib.hitad.aligner.DomainSet` instances
        
        Returns
        -------
        cache : dict
            Hierarchically organized matched domain pairs.
        
        """
        # First, align the top-level domains
        tg.getBottoms()
        qy.getBottoms()
        ttl = hierFormat(sorted(tg.bottoms))
        tql = hierFormat(sorted(qy.bottoms))
        ttg = DomainSet(tg.Label, ttl, tg.res, False)
        tqy = DomainSet(qy.Label, tql, qy.res, False)
        pairs, _ = self._aligncore(ttg, tqy, tg, qy)
        
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
    
    def _crosscorrect(self, cache, ref):
        
        tcore = {}
        for k in sorted(cache):
            sk = k[::-1]
            if not sk in ref:
                continue
            hit = ref[sk]
            target = cache[k]
            tcore[k] = Container(target.info)
            for lv in target:
                for p in target[lv]:
                    found = self._search(hit, p)
                    if found:
                        tl, ql = target[lv][p].info
                        mlv = min([d[-1] for d in tl])
                        if not mlv in tcore[k]:
                            tcore[k][mlv] = {p:target[lv][p]}
                        else:
                            tcore[k][mlv][p] = target[lv][p]
        return tcore
    
    def _search(self, hit, k):
        
        found = False
        for lv in sorted(hit):
            for h in hit[lv]:  
                if h[::-1] == k:
                    found = True
                    break
            if found:
                break
            
        return found
    
    def conserved(self, tn, qn):
        """
        Return conserved TAD pairs.
        
        Parameters
        ----------
        tn, qn : str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`.
        
        Returns
        -------
        pairs : set of tuples
            Each tuple has two elements (domain intervals), corresponding to
            *tn* and *qn* respectively.
        """
        pool = self.Results[tn][qn]
        pairs = set()
        for k in pool:
            tl, ql = pool[k].info
            if (len(tl)==1) and (len(ql)==1):
                pairs.add(k)
        pairs = pairs - self.inner_changed(tn,qn)
        
        return pairs
    
    def _lowlevel_changed(self, tn, qn):
        
        pool = self.Results[tn][qn]
        dset = self.DomainSets[tn]
        pairs = set()
        for k in sorted(pool):
            tl, ql = pool[k].info
            if (len(tl)>1) or (len(ql)>1):
                continue
            alls = dset.getregion(*tl[0][:3])
            labels = {tuple(d):0 for d in alls if d[-1]>0}
            if not len(labels):
                continue
            for lv in pool[k]:
                for t in pool[k][lv]:
                    tl, ql = pool[k][lv][t].info
                    if (len(tl)==1) and (len(ql)==1) and (tl[0][-1]==ql[0][-1]):
                        labels[tuple(tl[0])] = 1
            if not all(labels.values()):
                pairs.add(k)
            
        return pairs
    
    def inner_changed(self, tn ,qn):
        """
        Return semi-conserved TAD pairs.
        
        Parameters
        ----------
        tn, qn : str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`.
        
        Returns
        -------
        pairs : set of tuples
            Each tuple has two elements (domain intervals), corresponding to
            *tn* and *qn* respectively.
        """
        pairs = self._lowlevel_changed(tn, qn)
        r_pairs = self._lowlevel_changed(qn, tn)
        pairs.update(set([k[::-1] for k in r_pairs]))
        
        return pairs
    
    def merged(self, tn, qn):
        """
        Return merged region pairs and merged TAD details.
        
        Parameters
        ----------
        tn, qn : str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`.
        
        Returns
        -------
        pairs : dict
            The keys are merged region pairs in tuple, and the values are
            corresponding TAD list pairs within the region.
            
        """
        pool = self.Results[tn][qn]
        pairs = {}
        for k in pool:
            tl, ql = pool[k].info
            if (len(tl)>1) and (len(ql)==1):
                pairs[k] = [tl, ql]
                
        return pairs
    
    def split(self, tn, qn):
        """
        Return split region pairs and split TAD details.
        
        Parameters
        ----------
        tn, qn : str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`.
        
        Returns
        -------
        pairs : dict
            The keys are split region pairs in tuple, and the values are
            corresponding TAD list pairs within the region.
        """
        r_pairs = self.merged(qn, tn)
        pairs = {}
        for k in r_pairs:
            pairs[k[::-1]] = r_pairs[k][::-1]
        return pairs

class BoundAligner(DomainAligner):
    """
    Based on our hierarchical domain alignment scheme, we also define several
    change types on boundary level between two datasets, including conserved
    TAD boundary, conserved sub-TAD boundary, disappeared TAD boundary,
    disappeared sub-TAD boundary, and TAD-to-sub-TAD boundary switch.
    
    Boundaries are expressed in (chrom,pos) format in this class.
    
    Parameters
    ----------
    args : two or more :py:class:`tadlib.hitad.aligner.DomainSet` instances
    
    Attributes
    ----------
    byclass : dict
        Cache boundary pairs of each change type between datasets.
    
    """
    def __init__(self, *args):
        DomainAligner.__init__(self, *args)
        self.byclass = {}
        self.pairwise_alignment()
    
    def pairwise_alignment(self):
        
        completed = set()
        for tn in self.DomainSets:
            self.byclass[tn] = {}
            for qn in self.DomainSets:
                if tn==qn:
                    continue
                key = tuple(sorted((tn,qn)))
                if not key in completed:
                    self.align(tn, qn)
                    completed.add(key)
                self.byclass[tn][qn] = {}
                self.all_in_one(tn, qn, self.byclass[tn][qn])
    
    def all_in_one(self, tn, qn, cache):
        """
        Parse domain alignment results between *tn* and *qn* and cache all
        detected cases of 6 defined change types.
        
        Parameters
        ----------
        tn, qn : str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances.
        cache : dict
            An empty dictionary.
        """
        tad_c = {}
        sub_c = {}
        sub2tad = {}
        tad2sub = {}
        tad_dp = set()
        sub_dp = set()
        tset = self.DomainSets[tn]
        qset = self.DomainSets[qn]
        topbounds = [b for b in tset.Bounds if tset.boundclass[b]==0]
        lowbounds = [b for b in tset.Bounds if tset.boundclass[b]>0]
        pool = self.Results[tn][qn]
        for k in sorted(pool):
            tl, ql = pool[k].info
            tt_1 = (tl[0][0], tl[0][1]) # left target top bounds 
            tt_2 = (tl[-1][0], tl[-1][2])
            qt_1 = (ql[0][0], ql[0][1]) # left query top bounds
            qt_2 = (ql[-1][0], ql[-1][2])
            if (tset.boundclass[tt_1]==0) and (qset.boundclass[qt_1]==0):
                tad_c[tt_1] = qt_1
            elif (tset.boundclass[tt_1]>0) and (qset.boundclass[qt_1]>0):
                sub_c[tt_1] = qt_1
            elif (tset.boundclass[tt_1]==0) and (qset.boundclass[qt_1]>0):
                tad2sub[tt_1] = qt_1
            else:
                sub2tad[tt_1] = qt_1
            if (tset.boundclass[tt_2]==0) and (qset.boundclass[qt_2]==0):
                tad_c[tt_2] = qt_2
            elif (tset.boundclass[tt_2]>0) and (qset.boundclass[qt_2]>0):
                sub_c[tt_2] = qt_2
            elif (tset.boundclass[tt_2]==0) and (qset.boundclass[qt_2]>0):
                tad2sub[tt_2] = qt_2
            else:
                sub2tad[tt_2] = qt_2
            for lv in pool[k]:
                for t in pool[k][lv]:
                    tl, ql = pool[k][lv][t].info
                    tl_1 = (tl[0][0], tl[0][1]) # left target low bounds
                    tl_2 = (tl[-1][0], tl[-1][2])
                    ql_1 = (ql[0][0], ql[0][1]) # left query low bounds
                    ql_2 = (ql[-1][0], ql[-1][2])
                    if (tset.boundclass[tl_1]==0) and (qset.boundclass[ql_1]==0):
                        tad_c[tl_1] = ql_1
                    elif (tset.boundclass[tl_1]>0) and (qset.boundclass[ql_1]>0):
                        sub_c[tl_1] = ql_1
                    elif (tset.boundclass[tl_1]==0) and (qset.boundclass[ql_1]>0):
                        tad2sub[tl_1] = ql_1
                    else:
                        sub2tad[tl_1] = ql_1
                    if (tset.boundclass[tl_2]==0) and (qset.boundclass[ql_2]==0):
                        tad_c[tl_2] = ql_2
                    elif (tset.boundclass[tl_2]>0) and (qset.boundclass[ql_2]>0):
                        sub_c[tl_2] = ql_2
                    elif (tset.boundclass[tl_2]==0) and (qset.boundclass[ql_2]>0):
                        tad2sub[tl_2] = ql_2
                    else:
                        sub2tad[tl_2] = ql_2
        tad_dp = set(topbounds) - set(tad_c) - set(tad2sub)
        sub_dp = set(lowbounds) - set(sub_c) - set(sub2tad)
        cache['Conserved TAD boundaries'] = tad_c
        cache['Conserved sub-TAD boundaries'] = sub_c
        cache['TAD to sub-TAD'] = tad2sub
        cache['sub-TAD to TAD'] = sub2tad
        cache['Disappeared TAD'] = tad_dp
        cache['Disappeared sub-TAD'] = sub_dp
    
    def conserved_tad_bounds(self, tn, qn):
        """
        Return pairs of conserved TAD boundaries.
        
        Parameters
        ----------
        tn, qn: str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances.
        
        Returns
        -------
        pairs : dict
            Keys and values indicate TAD boundaries in *tn* and *qn*,
            respectively.
        """
        return self.byclass[tn][qn]['Conserved TAD boundaries']
    
    def conserved_sub_bounds(self, tn, qn):
        """
        Return pairs of conserved sub-TAD boundaries.
        
        Parameters
        ----------
        tn, qn: str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances.
        
        Returns
        -------
        pairs : dict
            Keys and values indicate sub-TAD boundaries in *tn* and *qn*,
            respectively.
        """
        return self.byclass[tn][qn]['Conserved sub-TAD boundaries']
    
    def tad2sub(self, tn, qn):
        """
        Return TAD to sub-TAD switch cases.
        
        Parameters
        ----------
        tn, qn: str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances.
        
        Returns
        -------
        pairs : dict
            Keys are TAD boundaries in *tn*, and values indicate corresponding
            sub-TAD boundaries in *qn*.
        """
        return self.byclass[tn][qn]['TAD to sub-TAD']
    
    def sub2tad(self, tn, qn):
        """
        Return sub-TAD to TAD switch cases.
        
        Parameters
        ----------
        tn, qn: str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances.
        
        Returns
        -------
        pairs : dict
            Keys are sub-TAD boundaries in *tn*, and values indicate corresponding
            TAD boundaries in *qn*.
        """
        return self.byclass[tn][qn]['sub-TAD to TAD']
    
    def disappeared_tad(self, tn, qn):
        """
        TAD boundaries that exist in *tn*, but disappear in *qn*.
        
        Parameters
        ----------
        tn, qn: str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances.
        
        Returns
        -------
        pairs : set of tuples
            TAD boundary positions in *tn*.
        """
        return self.byclass[tn][qn]['Disappeared TAD']
    
    def disappeared_sub(self, tn, qn):
        """
        Sub-TAD boundaries that exist in *tn*, but disappear in *qn*.
        
        Parameters
        ----------
        tn, qn: str
            Unique identifiers of :py:class:`tadlib.hitad.aligner.DomainSet`
            instances.
        
        Returns
        -------
        pairs : set of tuples
            Sub-TAD boundary positions in *tn*.
        """
        return self.byclass[tn][qn]['Disappeared sub-TAD']
        
#------------------------------------------------------------------------------
# Functions for parsing domain data        
def readHierDomain(domainfile, pre=''):
    """
    Load hierarchical domain list from a text file.
    
    The source file should contain 4 columns indicating chromosome label
    (1,2,...,X,Y), domain start (bp), domain end (bp), and hierarchical level
    (0,1,2,...), respectively.
    
    In our paper, TAD is denoted as level 0, sub-TAD is denoted as level 1,
    and subsequent domain level is denoted as level 2, etc.
    
    Parameters
    ----------
    domainfile : str
        Domain file path.
    
    Returns
    -------
    domainlist : list
        Each element of the list indicates one domain represented by
        [chrom,start,end,level].
    """
    domainlist = []
    with open(domainfile) as source:
        for line in source:
            parse = line.rstrip().split()
            chrom, start, end, label = [parse[i] for i in range(4)]
            chrom = chrom.lstrip(pre)
            if int(end) - int(start) <= 0:
                continue
            domainlist.append([chrom, int(start), int(end), int(label)])
    
    return domainlist

def readPlainDomain(domainfile, pre='chr'):
    """
    Load domain list from a text file.
    
    The source file should contain 3 columns indicating chromosome name,
    domain start (bp) and domain end (bp), respectively.
    
    Parameters
    ----------
    domainfile : str
        Domain file path.
    pre : str
        Leading string of the chromosome name. (Default: chr)
    
    Returns
    -------
    domainlist : list
        Each element indicates one domain represented by
        [chrom(leading string removed),start,end].
    
    See Also
    --------
    tadlib.hitad.aligner.hierFormat : parse hierarchical relationships between
                                      domains
    """
    domainlist = []
    with open(domainfile, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom, start, end = [parse[i] for i in range(3)]
            chrom = chrom.lstrip(pre)
            if int(end) - int(start) <= 0:
                continue
            domainlist.append([chrom, int(start), int(end)])
            
    return domainlist
    

def hierFormat(domainlist):
    """
    Resolve the nested/hierarchical relationships between domains, and
    transform the input [chrom,start,end] format domains into a format
    including hierarchical level information.
    
    Parameters
    ----------
    domainlist : list
        Domains with the format [chrom,start,end].
    
    Returns
    -------
    domainlist : list
        Domains with the format [chrom,start,end,level].
    
    """
    bychroms = {}
    hierlabel = {}
    for chrom, start, end in domainlist:
        key = (start, end)
        if not chrom in bychroms:
            bychroms[chrom] = [key]; hierlabel[chrom] = {key:0}
        else:
            bychroms[chrom].append(key)
            hierlabel[chrom][key] = 0
    for chrom in bychroms:
        bychroms[chrom].sort()
    nested = {}
    for c in bychroms:
        sidx = 0
        nested[c] = {}
        for q in bychroms[c]:
            pool = []
            label = 1
            for i in range(sidx, len(bychroms[c])):
                if (bychroms[c][i][0]>=q[0]) and (bychroms[c][i][1]<=q[1]):
                    pool.append(bychroms[c][i])
                if (bychroms[c][i][0]>=q[0]) and label:
                    sidx = i
                    label = 0
                if bychroms[c][i][0]>q[1]:
                    break
            nested[c][q] = pool
            for p in pool:
                if p != q:
                    hierlabel[c][p] += 1
    domainlist = []
    for c in hierlabel:
        for d in hierlabel[c]:
            line = [c, d[0], d[1], hierlabel[c][d]]
            domainlist.append(line)
    domainlist.sort()
    
    return domainlist