# -*- coding: utf-8 -*-
"""
Created on Thu Nov 03 16:57:05 2016

@author: wxt
"""

class Graph(dict):
    """
    A collection of graph methods required in our hierarchical domain and
    boundary alignment algorithm.
    
    Parameters
    ----------
    vs : list
        List of graph vertices.
    es : dict
        The keys are tuples of vertex pairs representing edges, and the values
        are corresponding weights.
        
    """
    def __init__(self, vs=[], es={}):
    
        for v in vs:
            self.add_vertex(v)
        
        for e, s in es.items():
            self.add_edge(e, s)
        
    def add_vertex(self, v):
        """
        Add vertex *v* to the graph.
        
        Parameters
        ----------
        v : any immutable object that can be used as a dictionary key
            A single vertex representation.
        """
        self[v] = {}
    
    def add_edge(self, e, s):
        """
        Bidirectionally link two vertices of *e*, and assign the attribute
        *s* (weight) to the edge.
        
        Parameters
        ----------
        e : tuple
            A single edge representation.
        s : any object
            Edge attribute.
        """
        v, w = e
        self[v][w] = s
        self[w][v] = s
    
    def find_connected_components(self):
        """
        Find all connected components of the graph using depth-first search.
        """
        pool = set()
        components = []
        for v in self:
            if not v in pool:
                cache = set()
                self.DFS(v, cache)
                components.append(cache)
                pool.update(cache)
        return components
    
    def DFS(self, v, cache):
        """
        Depth-first search for all vertices reachable from *v*.
        
        Parameters
        ----------
        v : any immutable object that can be used as a dictionary key
            A single vertex in the graph.
        cache : set
            Cache the visited vertices.
        """
        if not v in cache:
            cache.add(v)
        
        for w in self[v]:
            if not w in cache:
                self.DFS(w, cache)