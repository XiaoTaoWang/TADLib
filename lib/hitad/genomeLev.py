# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:36:44 2016

@author: wxt
"""
from __future__ import division
import glob, re, os, sys, cPickle, tempfile, time, zipfile, logging
import multiprocessing as mp
import numpy as np
from chromLev import Chrom, MultiReps
from numpy.lib.format import write_array

total_cpu = mp.cpu_count()

log = logging.getLogger(__name__)

class Genome(object):
    """
    *Genome* is built on top of :py:class:`tadlib.hitad.chromLev.Chrom`. We
    use it to:
        
    - Load bin-level HiC data
    - Initialize, pickle and organize *Chrom* objects
    - Call hierarchical domains of each chromosome in parallel
    - Optionally save HiC data into *Numpy* .npz files for accelerating
      IO in future.
    
    Parameters
    ----------
    datasets : 2-level dict, {resolution(int):{biological_replicate_label(str):data_path,...}}
        *resolution* should be in base-pair unit. *data_path* indicates the
        absolute HiC data path under corresponding resolution and biological
        replicate label.
        
        If your HiC data are stored in *NPZ* format, *data_path* should point
        to the npz file. Otherwise, you may provide with data in *TXT* format,
        in this case, HiC data of each chromosome must be stored separately,
        and data with the same resolution and replicate label should be placed
        in the same folder, and naturally *data_path* should point to the folder.
        
        You can generate *NPZ* files in two ways:1.By runHiC pipeline. runHiC
        is a user-friendly command-line software developed by our lab for HiC
        data processing. Refer to the `link <https://github.com/XiaoTaoWang/HiC_pipeline>`_
        for more details. 2.By hitad itself, provide TXT HiC data and run
        *hitad* with ``--npzpre`` specified, or, if you are familiar with Python
        environment, just open a Python interpreter and initialize a *Genome*
        object and set the parameter *npzpre* explicitly.
    
    maxsize : int
        Maximum allowable domain size in base-pair unit. (Default: 4000000)
    
    chroms : list
        List of chromosome labels. Only HiC data within the specified chromosomes
        will be included. Specially, '#' stands for chromosomes with numerical
        labels. If an empty list is provided, all chromosome data will be loaded.
        (Default: ['#', 'X'])
    
    npzpre : str or None
        If not None, loaded HiC data will be stored in *Numpy* .npz format
        to accelerate IO extremely for further use, and this parameter indicates
        the prefix of these *NPZ* filenames. Path may be included, if no path
        is contained, *NPZ* files will be placed under current working directory.
        (Default: None)
    
    cache : str or None
        :py:class:`tadlib.hitad.chromLev.Chrom` objects each representing a
        single chromosome data under certain resolution and replicate label
        will be pickled (using :py:mod:`cPickle`) under the folder named
        *cache*. The folder will be created if specified but don't exist.
        If None, the folder returned by :py:func:`tempfile.gettempdir` will
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
    
    def __init__(self, datasets, chroms=['#','X'], maxsize=4000000, npzpre=None, cache=None):
        
        self.chroms = set(chroms)
        data = datasets
        
        self._npzpre = npzpre
        if not self._npzpre is None:
            self._npzpre = os.path.abspath(os.path.expanduser(npzpre))
            for res in data:
                for rep in data[res]:
                    rl = '%dK' % (res//1000)
                    output = '.'.join([self._npzpre, rl, rep, 'npz'])
                    if os.path.exists(output):
                        log.error('The destination npz file will be overriden, reset npz prefix and run again ...')
                        log.error('Exit ...')
                        sys.exit(1)
        
        # We don't read data in memory at this point.
        # We only construct the mapping for loading convenience
        self.data = {}
        for res in data:
            for rep in data[res]:
                if data[res][rep].endswith('.npz'):
                    lib = np.load(data[res][rep])
                    for i in lib.files:
                        if ((not self.chroms) or (i.isdigit() and '#' in self.chroms)
                           or (i in self.chroms)):
                            if not i in self.data:
                                self.data[i] = {res:{rep:lib}}
                            else:
                                if res in self.data[i]:
                                    self.data[i][res][rep] = lib
                                else:
                                    self.data[i][res] = {rep:lib}
                else:
                    Map = self._scanFolder(data[res][rep])
                    for i in Map:
                        if not i in self.data:
                            self.data[i] = {res:{rep:Map[i]}}
                        else:
                            if res in self.data[i]:
                                self.data[i][res][rep] = Map[i]
                            else:
                                self.data[i][res] = {rep:Map[i]}
            
        
        if cache is None:
            self._cache = tempfile.gettempdir()
        else:
            self._cache = os.path.abspath(os.path.expanduser(cache))
            if not os.path.isdir(cache):
                os.makedirs(cache)
        
        self._intertype = np.dtype({'names':['bin1', 'bin2', 'IF'],
                                    'formats':[np.int, np.int, np.float]})
        
        # Before starting calling, we make cache data under the cache folder
        for chrom in self.data:
            log.debug('Chrom %s:', chrom)
            ms = self.data[chrom]
            for res in ms:
                for rep in ms[res]:
                    log.debug('  resolution: %d, %s', res, rep)
                    if type(ms[res][rep])==str:
                        tdata = np.loadtxt(ms[res][rep], dtype = self._intertype)
                    else:
                        tdata = ms[res][rep][chrom]
                    work = Chrom(chrom, res, tdata, rep, maxsize)
                    work.Label = rep
                    tl = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
                    kw = {'suffix':tl, 'dir':self._cache}
                    fd, tmpfil = tempfile.mkstemp(**kw)
                    os.close(fd)
                    with open(tmpfil, 'wb') as output:
                        log.debug('  Cache Chrom object into %s ...', tmpfil)
                        cPickle.dump(work, output, protocol = 2)
                    self.data[chrom][res][rep] = tmpfil
                    
                    time.sleep(3)
                    
                    if not self._npzpre is None:
                        rl = '%dK' % (res//1000)
                        output = '.'.join([self._npzpre, rl, rep, 'npz'])
                        if not os.path.exists(output):
                            Zip = zipfile.ZipFile(output, mode = 'w', allowZip64 = True)
                        else:
                            Zip = zipfile.ZipFile(output, mode = 'a', allowZip64 = True)
                        tl = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
                        fd, tmpfile = tempfile.mkstemp(suffix = '.'.join([tl, 'npy']))
                        os.close(fd)
                        log.debug('  Save data into %s for accelerating IO next time ...', output)
                        fname = '.'.join([chrom, 'npy'])
                        fid = open(tmpfile, 'wb')
                        try:
                            write_array(fid, tdata)
                            fid.close()
                            fid = None
                            Zip.write(tmpfile, arcname = fname)
                        finally:
                            if fid:
                                fid.close()
                        
                        os.remove(tmpfile)
                        Zip.close()
    
    def _extractChrLabel(self, filename):
        """
        Extract chromosome label from a file name.
        """
        # Full filename including path prefix
        _, interName = os.path.split(filename)
        regexp = 'chr(.*).txt'
        search_results = re.search(regexp, filename)
        label = search_results.group(1)
        
        # Remove leading zeroes.
        if label.isdigit():
            label = str(int(label))
    
        return label

    def _scanFolder(self, folder):
        """
        Create a map from chromosome labels to file names under the folder.
        """
        oriFiles = glob.glob(os.path.join(folder, 'chr*.txt'))
        
        # Read chromosome labels
        labels = []
        interFiles = [] # Depend on user's selection
        for i in oriFiles:
            label = self._extractChrLabel(i)
            if ((not self.chroms) or (label.isdigit() and '#' in self.chroms)
                or (label in self.chroms)):
                labels.append(label)
                interFiles.append(i)
        
        # Map from labels to files
        Map = dict(zip(labels, interFiles))
        
        return Map
        
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
                        log.debug('%s: completed', proc_name)
                        self.task_queue.task_done()
                        break
                    chrom, res, rep = current
                    log.debug('%s: Chrom %s (res %d, %s)', proc_name, chrom, res, rep)
                    worker(chrom, res, rep)
                    self.task_queue.task_done()
                
        def worker(chrom, res, rep):
            tmpfil = self.data[chrom][res][rep]
            with open(tmpfil, 'rb') as source:
                curChrom = cPickle.load(source)
            curChrom.callDomain()
            # Update cache data
            with open(tmpfil, 'wb') as output:
                cPickle.dump(curChrom, output, protocol = 2)
            
        cpu_core = min(cpu_core, total_cpu)
        
        log.debug('Spawn %d subprocesses ...', cpu_core)
        
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
            log.debug('Chrom %s ...', chrom)
            pool = {}
            for res in self.data[chrom]:
                reps = {}
                for rep in self.data[chrom][res]:
                    tmpfil = self.data[chrom][res][rep]
                    with open(tmpfil, 'rb') as source:
                        reps[rep] = cPickle.load(source)
                mrep = MultiReps(chrom, res, reps)
                mrep.callDomain()
                pool[res] = mrep
            minres = min(pool)
            self.Results.extend(pool[minres].mergedDomains)
    
    def outputDomain(self, filename):
        
        with open(filename, 'wb') as output:
            for d in self.Results:
                line = '%s\t%d\t%d\t%d\n' % tuple(d)
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
    