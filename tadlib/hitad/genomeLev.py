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
from chromLev import Chrom, MultiReps

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
    
    def __init__(self, datasets, maxsize=4000000, cache=None):
        
        data = datasets
        
        # We don't read data in memory at this point.
        # We only construct the mapping for loading convenience
        self.data = {}
        for res in data:
            for rep in data[res]:
                lib = cooler.Cooler(data[res][rep])
                for i in lib.chromnames:
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
                    if 'weight' in ms[res][rep].bins().keys(): # ICE correction
                        tdata = triu(ms[res][rep].matrix(balance=True, sparse=True).fetch(chrom).tocsr())
                    else:
                        tdata = triu(ms[res][rep].matrix(balance=False, sparse=True).fetch(chrom).tocsr())
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
    
    def outputDomain(self, filename):
        
        with open(filename, 'wb') as output:
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
    