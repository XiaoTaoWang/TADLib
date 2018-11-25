# Created on Wed Sep 24 15:05:28 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University
from __future__ import division
import numpy as np
from scipy.spatial import ConvexHull

class Polygon(ConvexHull):
    """Main class for polygon creation, manipulation and calulation.
    
    The API is inherited from :class:`scipy.spatial.ConvexHull`, which
    constructs convex hull from a point set. So you should input coordinates
    of a point set to initiate an instance.
    
    More general, not all points are required. Vertices alone are enough
    to construct the convex hull.
    
    Also, this class is designed to operate on 2-D space, although
    **ConvexHull** is okay for higher dimensions.
    
    Parameters
    ----------
    points : array_like
        Coordinates of points to construct a convex hull from. All objects
        that can be converted to a ndarray shaped (npoints, 2) are okay.
    
    Attributes
    ----------
    Common attributes from **ConvexHull**:
    
    points : ndarray of double, shape (npoints, 2)
        Coordinates of input points. Converted ndarray.
    vertices : ndarray of ints, shape (nvertices,)
        Indices of points forming the vertices of the convex hull. The
        vertices are counter-clockwise ordered.
    simplices : ndarray of ints, shape (nfacet, 2)
        Indices of points forming the simplical facets of the convex hull.
    
    Customized attributes:
    
    anchors : ndarray of double, shape (nvertices,)
        Vertices of the convex hull, counter-clockwise ordered.
    
    area : float
        Area of the polygon.
    
    Examples
    --------
    Convex hull of a random set of points:
    
    >>> import numpy as np
    >>> from tadlib.calfea.polygon import Polygon
    >>> points = np.random.rand(20, 2) # 20 random points in 2-D space
    >>> P = Polygon(points)
    >>> print P.anchors
    [[ 0.71119071  0.216714  ]
     [ 0.92908323  0.77903127]
     [ 0.79530032  0.93279735]
     [ 0.02160083  0.98758203]
     [ 0.09094998  0.4423167 ]
     [ 0.16924963  0.05900692]
     [ 0.71119071  0.216714  ]]
    
    """
    
    def __init__(self, points):
        """A customized constructor.
        """
        ConvexHull.__init__(self, points)
        indices = self.vertices
        points = self.points
        # Customize
        self.anchors = points[indices]
    
    def calarea(self):
        """Calculate the polygon area.
        
        An attribute called **area** is assigned.
        
        See Also
        --------
        tadlib.calfea.polygon.shoelace : Twice the area of polygon
        
        """
        # Our Calculation
        # Call external function
        self.area = abs(shoelace(self.anchors)) / 2.0
    
    def close(self):
        """
        Close the polygon, i.e., the first point is also the last one.
        
        See Also
        --------
        tadlib.calfea.polygon.isinside : judge if points are inside a polygon
                                         or not.
        
        Notes
        -----
        Must be called before :py:meth:`tadlib.calfea.polygon.Polygon.isinside`.
        
        """
        first = self.anchors[0]
        last = self.anchors[-1]
        if not np.all(first==last):
            # Concatenate along first axis, dim>=2
            self.anchors = np.r_['0,2', self.anchors, self.anchors[0]]
    
    def isinside(self, query, zerolike=1e-12):
        """Compute point location relative to a polygon.
        
        Parameters
        ----------
        query : tuple or array_like
            A tuple indicates one point. (The length must be 2)
            Or you can input an array-like sequence. Any object that can be
            converted to a ndarray shaped (npoints, 2) is okay.
        
        zerolike : float
            A number used to approximate 0.

        Returns
        -------
        mindst : scalar or array_like
            If mindst < 0, point is outside the polygon.
            If mindst = 0, point in on a side of the polygon.
            If mindst > 0, point is inside the polygon.

        Notes
        -----
        Sloan's improved version of the Nordbeck and Rydstedt algorithm.
        
        Examples
        --------
        We start with :class:`tadlib.calfea.polygon.Polygon` construction:
        
        >>> import numpy as np
        >>> from tadlib.calfea.polygon import Polygon
        >>> points = np.random.rand(20, 2) # Used for constructing Polygon
        >>> P = Polygon(points)
        >>> check = np.random.rand(3, 2) # Another 3 random points
        >>> P.isinside(check)
        array([ 0.21145311,  0.09807244, -0.15341914])

        """
        # Indicator
        label = type(query)==tuple
        # Unify the interface
        query = np.array(query, ndmin=2)
        # Close the polygon
        self.close()
        
        # If snear = True, dist to nearest side < nearest vertex
        # If snear = False, dist to nearest vertex < nearest side
        snear = np.ma.masked_all(query.shape[0], dtype=bool)
        # Initialization
        mindst = np.ones(query.shape[0], dtype=float) * np.inf
        j = np.ma.masked_all(query.shape[0], dtype=int)
        # Number of sides/vertices defining the polygon
        n = len(self.anchors) - 1
        # Loop over each side of the polygon
        for i in range(n):
            d = np.ones(query.shape[0], dtype=float) * np.inf
            # Vertex (x1, y1), start of side
            start = self.anchors[i]
            # Vertex (x2, y2), end of side
            sec = self.anchors[i + 1] - start
            # Query has coordinates (qx, qy)
            devs = self.anchors[i] - query
            # Points on infinite line defined by
            #     x = x1 + t * (x1 - x2)
            #     y = y1 + t * (y1 - y2)
            # Find where normal passing through (qx, qy) intersects
            # infinite line
            t = -(devs[:,0] * sec[0] + devs[:,1] * sec[1]) / \
                 (sec[0] ** 2 + sec[1] ** 2)
            tlt0 = t < 0
            tle1 = (0 <= t) & (t <= 1)
            # Normal intersects side
            d[tle1] = ((devs[:,0][tle1] + t[tle1] * sec[0]) ** 2 + 
                       (devs[:,1][tle1] + t[tle1] * sec[1]) ** 2)
            # Normal does not intersects side
            # Point is closest to vertex (x1, y1)
            # Compute square of distance to this vertex
            d[tlt0] = devs[:,0][tlt0] ** 2 + devs[:,1][tlt0] ** 2
            # Store distances
            mask = d < mindst
            mindst[mask] = d[mask]
            j[mask] = i
            # Point is closer to (x1, y1) than any other vertex or side
            snear[mask & tlt0] = False
            # Point is closer to this side than to any other side or vertex
            snear[mask & tle1] = True
        if np.ma.count(snear) != snear.size:
            raise IndexError('Error when computing distances')
        mindst **= 0.5
        # Point is closer to its nearest vertex than its nearest side,
        # check if nearest vertex is concave.
        # If the nearest vertex is concave then point is inside the polygon,
        # else the point is outside the polygon.
        jo = j.copy()
        jo[j == 0] -= 1
        area = shoelace([self.anchors[j + 1], self.anchors[j],
                         self.anchors[jo - 1]])
        mindst[~snear] = np.copysign(mindst, area)[~snear]
        # Point is closer to its nearest side than to its nearest vertex,
        # check if point is to left or right of this side.
        # If point is to left of side it is inside polygon, else point is
        # outside polygon.
        area = shoelace([self.anchors[j], self.anchors[j + 1], query])
        mindst[snear] = np.copysign(mindst, area)[snear]
        # Point is on side of polygon
        mindst[np.fabs(mindst) < zerolike] = 0
        # If input values were scalar then the output should be too
        if label:
            mindst = float(mindst)
        
        return mindst
        
def shoelace(vertices):
    """        
    Calculate twice the area of polygon using Shoelace formula.
    
    Polygon is defined by vertices.
    
    Parameters
    ----------
    vertices : array_like
        Vertex coordinates in a 2-D space.
        Coordinates must be placed along the last axis. And data points are
        along the first axis.
    
    Returns
    -------
    area : float
        You can deduce the order of input vertices from the sign:
        area is positive if vertices are in counter-clockwise order.
        area is negative if vertices are in clockwise order.
        area is zero if all points are colinear.
    
    Notes
    -----
    This function can be also used to judge if all points in a data set are
    collinear. Collinear points as input for initializing Polygon instance
    will raise a QhullError.
    
    Examples
    --------
    Vertices of a square:

    Clockwise:    
    
    >>> from tadlib.calfea.polygon import shoelace
    >>> sq = [(0,0), (0,1), (1,1), (1,0)]
    >>> shoelace(sq)
    -2.0
    
    Counter-clockwise:
    
    >>> sq = [(0,0), (1,0), (1,1), (0,1)]
    >>> shoelace(sq)
    2.0
        
    """
    
    vertices = np.asfarray(vertices)
    # Rule for stacking multiple comma separated arrays
    rule = '0,' + str(len(vertices.shape))
    # Slip the array along the first axis
    slip_v = np.r_[rule, vertices[-1], vertices[:-1]]
    # Extract coordinates
    x = np.take(vertices, [0], axis=-1).reshape(vertices.shape[:-1])
    y = np.take(vertices, [1], axis=-1).reshape(vertices.shape[:-1])
    slip_x = np.take(slip_v, [0], axis=-1).reshape(vertices.shape[:-1])
    slip_y = np.take(slip_v, [1], axis=-1).reshape(vertices.shape[:-1])
    # Sholelace Foluma
    area = np.sum(y * slip_x - x * slip_y, axis=0)
    
    return area

def collinear(points):
    """Test whether all given points are collinear.
    
    Collinear points will trigger an error called **QhullError** when used
    to initialize a :class:`tadlib.calfea.polygon.Polygon` instance. However,
    other conditions may also trigger **QhullError**. Doing this test in
    advance can avoid this error and make things clearer.
    
    Parameters
    ----------
    points : array_like
        Coordinates of points to construct a polygon. Any object that can 
        be converted to a ndarray shaped (npoints, 2) is okay.
    
    Returns
    -------
    judge : bool
        True if the input points are collinear else False.
    
    See Also
    --------
    tadlib.calfea.polygon.shoelace
    
    Notes
    -----
    To judge if all points are collinear, we use a simple heuristic
    algorithm. We sort the points at first. Then test whether consecutive
    triples are collinear.
    
    Examples
    --------
    Trival but still effective:
    
    >>> from tadlib.calfea.polygon import collinear
    >>> line = [(2, 0.4), (2, 0.8), (2, 4), (2, 100)]
    >>> collinear(line)
    True
    
    """
    
    if len(points)<3:
        judge = True # Points are always collinear if N<3
    else:
        points = np.asfarray(points)
        # Indices for sorting
        idx = np.argsort(points, axis=0)
        mask = points[:,0]==points[0][0]
        if not np.all(mask):
            asort = points[idx[:,0]]
        else:
            asort = points[idx[:,1]]
        # The first point of the triples
        first = np.arange(len(points)-2)
        triples = [asort[first], asort[first + 1], asort[first + 2]]
        # Calculate triple areas
        areas = shoelace(triples)
        # Triples are collinear if area=0
        # If all sorted triples are collinear, then all points are collinear
        judge = np.all(areas==0)
    
    return judge
