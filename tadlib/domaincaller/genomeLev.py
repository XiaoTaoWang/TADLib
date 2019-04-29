# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:36:44 2016

@author: wxt
"""
from __future__ import division
import os, pickle, tempfile, time, logging, cooler
import multiprocessing as mp
import numpy as np
from scipy.sparse import triu
from pomegranate import NormalDistribution, HiddenMarkovModel, GeneralMixtureModel, State
from tadlib.hitad import genomeLev
from tadlib.domaincaller.chromLev import Chrom

total_cpu = mp.cpu_count()

log = logging.getLogger(__name__)

genomeCore = genomeLev.Genome

class Genome(genomeCore):
    
    def __init__(self, uri, balance_type='weight', window=2000000, minsize=100000, cache=None,
                 exclude=['chrY','chrM'], DIcol='DIs'):

        if cache is None:
            self._cache = tempfile.gettempdir()
        else:
            self._cache = os.path.abspath(os.path.expanduser(cache))
            if not os.path.isdir(cache):
                os.makedirs(cache)
        
        self.exclude = exclude
        self.DIcol = DIcol

        if balance_type.lower() == 'raw':
            correct = False
        else:
            correct = balance_type
        
        lib = cooler.Cooler(uri)
        res = lib.binsize
        minsize = int(np.ceil(minsize/res))
        # Before starting calling, we make cache data under the cache folder
        self.data = {}
        self.chroms = []
        for c in lib.chromnames:
            if c in self.exclude:
                continue
            self.chroms.append(c)
            log.debug('Chrom {0}:'.format(c))
            if 'weight' in lib.bins().keys(): # ICE correction
                tdata = triu(lib.matrix(balance=correct, sparse=True).fetch(c)).tocsr()
            else:
                tdata = triu(lib.matrix(balance=correct, sparse=True).fetch(c)).tocsr()
            work = Chrom(c, res, tdata, window=window, minsize=minsize)
            tl = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
            kw = {'suffix':tl, 'dir':self._cache}
            fd, tmpfil = tempfile.mkstemp(**kw)
            os.close(fd)
            with open(tmpfil, 'wb') as output:
                log.debug('  Cache Chrom object into {0} ...'.format(tmpfil))
                pickle.dump(work, output)
            self.data[c] = tmpfil

            time.sleep(3)
        
        self.cool = lib
        
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
        """
        hmm = HiddenMarkovModel()
        # GMM emissions
        # 4 Hidden States:
        # 0--start, 1--downstream, 2--upstream, 3--end
        numdists = 3 # Three-distribution Gaussian Mixtures
        var = 7.5 / (numdists - 1)
        means = [[], [], [], []]
        for i in range(numdists):
            means[3].append(i * 7.5 / ( numdists - 1 ) + 2.5)
            means[2].append(i * 7.5 / ( numdists - 1 ))
            means[1].append(-i * 7.5 / ( numdists - 1 ))
            means[0].append(-i * 7.5 / ( numdists - 1 ) - 2.5)
        states = []
        for i, m in enumerate(means):
            tmp = []
            for j in m:
                tmp.append(NormalDistribution(j, var))
            mixture = GeneralMixtureModel(tmp)
            states.append(State(mixture, name=str(i)))
        hmm.add_states(*tuple(states))

        # Transmission matrix
        #A = [[0., 1., 0., 0.],
        #    [0., 0.5, 0.5, 0.],
        #    [0., 0., 0.5, 0.5],
        #    [1., 0., 0., 0.]]
        hmm.add_transition(states[0], states[1], 1)
        hmm.add_transition(states[1], states[1], 0.5)
        hmm.add_transition(states[1], states[2], 0.5)
        hmm.add_transition(states[2], states[2], 0.5)
        hmm.add_transition(states[2], states[3], 0.5)
        hmm.add_transition(states[3], states[0], 1)

        #pi = [0.2, 0.3, 0.3, 0.2]
        hmm.add_transition(hmm.start, states[0], 1)
        hmm.add_transition(states[3], hmm.end, 1)

        hmm.bake()

        return hmm

    def train_data(self):

        seqs = []
        for c in self.chroms:
            tmpfil = self.data[c]
            with open(tmpfil, 'rb') as source:
                tmpcache = pickle.load(source)
            tmpcache.calDI(tmpcache.windows, 0)
            tmpcache.splitChrom(tmpcache.DIs)
            for region in tmpcache.regionDIs:
                seqs.append(tmpcache.regionDIs[region])
        
        return seqs
    
    def learning(self, cpu_core = 1):
        
        model = self.oriHMMParams()
        seqs = self.train_data()
        model.fit(seqs, algorithm='baum-welch', max_iterations=10000,
                  stop_threshold=1e-5, n_jobs=cpu_core, verbose=False)
        for c in self.chroms:
            tmpfil = self.data[c]
            with open(tmpfil, 'rb') as source:
                tmpcache = pickle.load(source)
            tmpcache.hmm = model
            # same hmm for all chromosomes
            with open(tmpfil, 'wb') as output:
                pickle.dump(tmpcache, output)
        
    def callDomains(self, cpu_core = 1):
        
        class SubProcess(mp.Process):
            
            def __init__(self, task_queue):
                mp.Process.__init__(self)
                self.task_queue = task_queue
            
            def run(self):
                proc_name = self.name
                while True:
                    current = self.task_queue.get()
                    if current is None:
                        log.debug('{0}: completed'.format(proc_name))
                        self.task_queue.task_done()
                        break
                    log.debug('{0}: Chrom {1}'.format(proc_name, current))
                    worker(current)
                    self.task_queue.task_done()
                
        def worker(chrom):
            tmpfil = self.data[chrom]
            with open(tmpfil, 'rb') as source:
                curChrom = pickle.load(source)
            curChrom.callDomains()
            # Update cache data
            with open(tmpfil, 'wb') as output:
                pickle.dump(curChrom, output)
            
        cpu_core = min(cpu_core, total_cpu)
        
        log.debug('Spawn {0} subprocesses ...'.format(cpu_core))
        
        tasks = mp.JoinableQueue()
        procs = [SubProcess(tasks) for i in range(cpu_core)]
        for p in procs:
            p.start()
        for chrom in self.data:
            tasks.put(chrom)
        for p in procs:
            tasks.put(None)
        
        tasks.join()

        self.Results = []
        for chrom in self.data:
            tmpfil = self.data[chrom]
            with open(tmpfil, 'rb') as source:
                curChrom = pickle.load(source)
            for d in curChrom.domains:
                if d[2] > 0.5:
                    continue
                self.Results.append([chrom, d[0], d[1], 0])
        
        log.debug('Add the DI column into cools ...')
        DIs = np.r_[[]]
        for c in self.cool.chromnames:
            if not c in self.exclude:
                tmpfil = self.data[c]
                with open(tmpfil, 'rb') as source:
                    tmpcache = pickle.load(source)
                DIs = np.r_[DIs, tmpcache.DIs]
            else:
                DIs = np.r_[DIs, np.zeros(len(self.cool.bins().fetch(c)))]
        with self.cool.open('r+') as grp:
            if self.DIcol in grp['bins']:
                del grp['bins'][self.DIcol]
            h5opts = dict(compression='gzip', compression_opts=6)
            grp['bins'].create_dataset(self.DIcol, data=DIs, **h5opts)
    
    def outputDomain(self, filename):
        
        with open(filename, 'w') as output:
            for d in self.Results:
                line = '{0}\t{1}\t{2}\t{3}\n'.format(*tuple(d))
                output.write(line)
    
    def wipeDisk(self):
        """
        Remove catched (pickled) :py:class:`tadlib.hitad.chromLev.Chrom`
        objects before exiting.
        """
        for chrom in self.data:
            os.remove(self.data[chrom])
    