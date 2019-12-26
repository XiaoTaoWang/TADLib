# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 19:38:21 2018

@author: XiaoTao Wang
"""
import itertools, cooler
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

class Triangle(object):
    
    def __init__(self, uri, chrom, start, end, correct='weight', figsize=(7, 3.5)):

        self.clr = cooler.Cooler(uri)
        self.res = self.clr.binsize
        fig = plt.figure(figsize=figsize)
        self.fig = fig

        self.chrom = chrom
        self.start = start
        self.end = end

        M = self.clr.matrix(balance=correct, sparse=False).fetch((chrom, start, end))
        M[np.isnan(M)] = 0
        self.matrix = M

        # define my colormap (traditional w --> r)
        self.cmap = LinearSegmentedColormap.from_list('interaction',
                ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])

    
    def print_coordinate(self, pos):

        i_part = int(pos) // 1000000 # Integer Part
        d_part = (int(pos) % 1000000) // 1000 # Decimal Part
        
        if (i_part > 0) and (d_part > 0):
            return ''.join([str(i_part), 'M', str(d_part), 'K'])
        elif (i_part == 0):
            return ''.join([str(d_part), 'K'])
        else:
            return ''.join([str(i_part), 'M'])
    
    def matrix_plot(self, colormap='traditional', vmin=None, vmax=None, cbr_fontsize=9,
        nticks=4, label_size=9, remove_label=False, heatmap_pos=[0.1, 0.1, 0.8, 0.8],
        colorbar_pos=[0.08, 0.45, 0.02, 0.15], chrom_pos=[0.1, 0.08, 0.8, 0.015]):

        h_ax = self.fig.add_axes(heatmap_pos)
        c_ax = self.fig.add_axes(colorbar_pos)

        M = self.matrix
        n = M.shape[0]

        # Create the rotation matrix
        t = np.array([[1,0.5], [-1,0.5]])
        A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)

        if colormap=='traditional':
            cmap = self.cmap
        else:
            cmap = colormap
        
        # Plot the Heatmap ...
        x = A[:,1].reshape(n+1, n+1)
        y = A[:,0].reshape(n+1, n+1)
        y[y<0] = -y[y<0]

        if vmax is None:
            vmax = np.percentile(M[M.nonzero()], 95)
        if vmin is None:
            vmin = M.min()
        
        sc = h_ax.pcolormesh(x, y, np.flipud(M), vmin=vmin, vmax=vmax, cmap=cmap,
                        edgecolor='none', snap=True, linewidth=.001)
        
        # colorbar
        cbar = self.fig.colorbar(sc, cax=c_ax, ticks=[vmin, vmax], format='%.3g')
        c_ax.tick_params(labelsize=cbr_fontsize)


        # Hide the bottom part
        xmin = A[:,1].min()
        xmax = A[:,1].max()
        ymin = A[:,0].min()
        ymax = 0
        h_ax.fill([xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax], 'w', ec='none')
        h_ax.axis('off')

        # chromosome bar
        if not remove_label:
            chrom_ax = self.fig.add_axes(chrom_pos)
            chrom_ax.tick_params(axis='both', bottom=True, top=False, left=False,
                                right=False, labelbottom=True, labeltop=False,
                                labelleft=False, labelright=False)
            interval = (self.end - self.start) // self.res
            ticks = list(np.linspace(0, interval, nticks).astype(int))
            pos = list(np.linspace(self.start, self.end, nticks).astype(int))
            labels = [self.print_coordinate(p) for p in pos]
            chrom_ax.set_xticks(ticks)
            if len(ticks) < 7:
                chrom_ax.set_xticklabels(labels, fontsize=label_size)
            else:
                chrom_ax.set_xticklabels(labels, fontsize=label_size, rotation=15, ha='right')
            
            chrom_ax.set_xlabel('chr'+self.chrom.lstrip('chr'), fontsize=label_size+2)

            chrom_ax.set_xlim(ticks[0], ticks[-1])
            chrom_ax.set_ylim(0, 0.02)

        self.heatmap_ax = h_ax
        self.cbar_ax = c_ax
        self.chrom_ax = chrom_ax
        self.hx = x
        self.hy = y
    

    def plot_loops(self, loop_file, marker_size=50, marker_color='#111111', marker_type='o',
        marker_alpha=0.5):

        loopType = np.dtype({'names':['chr','start1','end1','start2','end2'],
                         'formats':['U5', np.int, np.int, np.int, np.int]})
        loops = np.loadtxt(loop_file, dtype=loopType, usecols=[0,1,2,4,5])
        loops = loops[(loops['chr']==self.chrom)]

        test_x = loops['start1']
        test_y = loops['end2']
        mask = (test_x >= self.start) & (test_y < self.end)
        loops = loops[mask]

        n = self.matrix.shape[0]

        # mark the loop loci
        Bool = np.zeros((n, n), dtype=bool)
        for xs, xe, ys, ye in zip(loops['start1'], loops['end1'], loops['start2'], loops['end2']):
            # Lodate the loop pixel at given resolution
            s_l = range(xs//self.res-1, int(np.ceil(xe/float(self.res)))+1)
            e_l = range(ys//self.res-1, int(np.ceil(ye/float(self.res)))+1)
            si, ei = None, None
            for i in s_l:
                for j in e_l:
                    st = i - self.start//self.res
                    et = j - self.start//self.res
                    if (st < n) and (et < n):
                        if si is None:
                            si, ei = st, et
                        else:
                            if self.matrix[st,et] > self.matrix[si,ei]:
                                si, ei = st, et
            if not si is None:
                Bool[si, ei] = 1
                #Bool[ei, si] = 1
        
        lx = self.hx[:-1,:-1][np.flipud(Bool)]
        ly = self.hy[:-1,:-1][np.flipud(Bool)] + 1
        if lx.size > 0:
            self.heatmap_ax.scatter(lx, ly, s=marker_size, c=marker_color, marker=marker_type,
                alpha=marker_alpha)
        
        self.heatmap_ax.set_xlim(self.hx.min(), self.hx.max())
        self.heatmap_ax.set_ylim(self.hy.min(), self.hy.max())

        self.loops = loops
    
    def plot_TAD(self, tad_fil, line_color='#60636A', linewidth=3, line_style='-'):

        tadtype = np.dtype({'names':['chr','start','end'],
                            'formats':['U5', np.int, np.int]})
        tads = np.loadtxt(tad_fil, dtype=tadtype, usecols=[0,1,2])
        tads = tads[(tads['chr']==self.chrom)]
        mask = (tads['end'] > self.start) & (tads['start'] < self.end)
        tads = tads[mask]

        n = self.matrix.shape[0]

        for s, e in zip(tads['start'], tads['end']):
            si = s // self.res - self.start // self.res
            ei = e // self.res - self.start // self.res
            if si < 0:
                si = 0
            if ei > n - 1:
                ei = n - 1
            
            if ei - si < 2:
                continue

            x = [self.hx[:-1, :-1][n-1-si, si],
                 self.hx[:-1, :-1][n-1-si, ei],
                 self.hx[:-1, :-1][n-1-ei, ei]]
            y = [self.hy[:-1, :-1][n-1-si, si] - 1,
                 self.hy[:-1, :-1][n-1-si, ei] + 1,
                 self.hy[:-1, :-1][n-1-ei, ei] - 1]
            self.heatmap_ax.plot(x, y, color=line_color, linestyle=line_style,
                linewidth=linewidth)
        
        self.heatmap_ax.set_xlim(self.hx.min(), self.hx.max())
        self.heatmap_ax.set_ylim(self.hy.min(), self.hy.max())

        self.tads = tads
    

    def outfig(self, outfile, dpi=200, bbox_inches='tight'):

        self.fig.savefig(outfile, dpi=dpi, bbox_inches=bbox_inches)
    
    def show(self):

        self.fig.show()
