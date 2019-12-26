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
from tadlib.hitad.chromLev import Chrom, MultiReps

total_cpu = mp.cpu_count()

log = logging.getLogger(__name__)

class Genome(object):
    """
    *Genome* is built on top of :py:class:`tadlib.hitad.chromLev.Chrom`. We
    use it to:
        
    - Load bin-level Hi-C data
    - Initialize, pickle and organize *Chrom* objects
    - Call hierarchical domains of each chromosome in parallel
    
    Parameters
    ----------
    datasets : 2-level dict, {resolution(int):{biological_replicate_label(str):data_path,...}}
        *data_path* indicates the *cool* URI under corresponding resolution and biological
        replicate label.
    
    maxsize : int
        Maximum allowable domain size in base-pair unit. (Default: 4000000)
    
    cache : str or None
        Cache folder path. If None, the folder returned by :py:func:`tempfile.gettempdir` will
        be used. (Default: None)
    
    Attributes
    ----------
    data : 3-level dict. {chrom(str):{resolution(int):{biological_replicate_label(str):cachedfile,...}},...}
        Different from the input datasets, it organizes data by chromosome
        label, and each bottom-level value indicates one pickled
        :py:class:`tadlib.hitad.chromLev.Chrom` file under the *cache* folder.
    
    Results : list
        Final consistent domain list merged from all chromosome results.
    """
    
    def __init__(self, datasets, balance_type='weight', maxsize=4000000, cache=None,
                 exclude=['chrM', 'chrY'], DIcol='DIs'):
        
        data = datasets
        self.exclude = exclude
        self.DIcol = DIcol

        if balance_type.lower() == 'raw':
            correct = False
        else:
            correct = balance_type
        
        # We don't read data in memory at this point.
        # We only construct the mapping for loading convenience
        self.data = {}
        for res in data:
            for rep in data[res]:
                lib = cooler.Cooler(data[res][rep])
                for i in lib.chromnames:
                    if i in self.exclude:
                        continue
                    if not i in self.data:
                        self.data[i] = {res:{rep:lib}}
                    else:
                        if res in self.data[i]:
                            self.data[i][res][rep] = lib
                        else:
                            self.data[i][res] = {rep:lib}
        
        if cache is None:
            self._cache = tempfile.gettempdir()
        else:
            self._cache = os.path.abspath(os.path.expanduser(cache))
            if not os.path.isdir(cache):
                os.makedirs(cache)
        
        # Before starting calling, we make cache data under the cache folder
        for chrom in self.data:
            log.debug('Chrom {0}:'.format(chrom))
            ms = self.data[chrom]
            for res in ms:
                for rep in ms[res]:
                    log.debug('  resolution: {0}, {1}'.format(res, rep))
                    tdata = triu(ms[res][rep].matrix(balance=correct, sparse=True).fetch(chrom)).tocsr()
                    work = Chrom(chrom, res, tdata, rep, maxsize)
                    work.Label = rep
                    tl = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
                    kw = {'suffix':tl, 'dir':self._cache}
                    fd, tmpfil = tempfile.mkstemp(**kw)
                    os.close(fd)
                    with open(tmpfil, 'wb') as output:
                        log.debug('  Cache Chrom object into {0} ...'.format(tmpfil))
                        pickle.dump(work, output)
                    self.data[chrom][res][rep] = tmpfil

                    time.sleep(3)

        self.cools = datasets
    
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
    
    def train_data(self, res, rep):

        lib = cooler.Cooler(self.cools[res][rep])
        seqs = []
        for c in lib.chromnames:
            if c in self.exclude:
                continue
            tmpfil = self.data[c][res][rep]
            with open(tmpfil, 'rb') as source:
                tmpcache = pickle.load(source)
            tmpcache.minWindows(0, tmpcache.chromLen, tmpcache._dw)
            tmpcache.calDI(tmpcache.windows, 0)
            tmpcache.splitChrom(tmpcache.DIs)
            for region in tmpcache.regionDIs:
                withzeros = tmpcache.regionDIs[region]
                nozeros = withzeros[withzeros!=0]
                if nozeros.size > 20:
                    seqs.append(nozeros)
        
        return seqs
    
    def learning(self, cpu_core = 1):
        """
        Prepare training data and learn HMM model parameters for each dataset.

        Parameters
        ----------
        cpu_core : int
            Number of processes to launch.
        """
        for res in self.cools:
            for rep in self.cools[res]:
                log.debug('  resolution: {0}, {1}'.format(res, rep))
                model = self.oriHMMParams()
                seqs = self.train_data(res, rep)
                model.fit(seqs, algorithm='baum-welch', max_iterations=10000,
                          stop_threshold=1e-5, n_jobs=cpu_core, verbose=False)
                lib = cooler.Cooler(self.cools[res][rep])
                for c in lib.chromnames:
                    if c in self.exclude:
                        continue
                    tmpfil = self.data[c][res][rep]
                    with open(tmpfil, 'rb') as source:
                        tmpcache = pickle.load(source)
                    tmpcache.hmm = model
                    # same hmm for all chromosomes
                    with open(tmpfil, 'wb') as output:
                        pickle.dump(tmpcache, output)
        
    def callHierDomain(self, cpu_core = 1):
        """
        Identify hierarchical domains of each chromosome independently
        and concurrently, find consistent domains between biological
        replicates, and finally combine results of all chromosomes.
        
        Parameters
        ----------
        cpu_core : int
            Number of processes to launch. For now, *hitad* only supports
            parallel computing on a quite high layer, that is, it simply
            allocates an uncompleted :py:class:`tadlib.hitad.chromLev.Chrom`
            object to an idle processor and invokes its *callDomain* method.
        """
        
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
                    chrom, res, rep = current
                    log.debug('{0}: Chrom {1} (res {2}, {3})'.format(proc_name, chrom, res, rep))
                    worker(chrom, res, rep)
                    self.task_queue.task_done()
                
        def worker(chrom, res, rep):
            tmpfil = self.data[chrom][res][rep]
            with open(tmpfil, 'rb') as source:
                curChrom = pickle.load(source)
            curChrom.callDomain()
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
            for res in self.data[chrom]:
                for rep in self.data[chrom][res]:
                    tasks.put((chrom, res, rep))
        for p in procs:
            tasks.put(None)
        
        tasks.join()
        
        log.debug('Extract reproducible domains from replicates ...')
        self.Results = []
        self.DIs = {}
        for chrom in self.data:
            log.debug('Chrom {0} ...'.format(chrom))
            pool = {}
            for res in self.data[chrom]:
                reps = {}
                for rep in self.data[chrom][res]:
                    tmpfil = self.data[chrom][res][rep]
                    with open(tmpfil, 'rb') as source:
                        reps[rep] = pickle.load(source)
                mrep = MultiReps(chrom, res, reps)
                mrep.callDomain()
                pool[res] = mrep
            minres = min(pool)
            self.Results.extend(pool[minres].mergedDomains)
        
        log.debug('Add the DI column into cools ...')
        for res in self.cools:
            for rep in self.cools[res]:
                lib = cooler.Cooler(self.cools[res][rep])
                DIs = np.r_[[]]
                for c in lib.chromnames:
                    if not c in self.exclude:
                        tmpfil = self.data[c][res][rep]
                        with open(tmpfil, 'rb') as source:
                            tmpcache = pickle.load(source)
                        DIs = np.r_[DIs, tmpcache.DIs]
                    else:
                        DIs = np.r_[DIs, np.zeros(len(lib.bins().fetch(c)))]
                with lib.open('r+') as grp:
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
            for res in self.data[chrom]:
                for rep in self.data[chrom][res]:
                    os.remove(self.data[chrom][res][rep])
    