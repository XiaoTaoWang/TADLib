TAD Data Parsing
================
.. autofunction:: tadlib.calfea.analyze.load_TAD

Intra-TAD Interaction Analysis
==============================
.. autoclass:: tadlib.calfea.analyze.Core
   :members:

Matrix Manipulating
===================
.. autofunction:: tadlib.calfea.analyze.manipulation

Polygon Creation and Operations
===============================
Here, polygon means convex polygon. In the simplest case, a polygon can be
created or presented by vertices. More ofen, you don't know the vertices but
want a polygon enclosing a set of points.

**ConvexHull** defined in **scipy.spatial** module uses the Qhull library to
compute convex hull for a finite set of points, but doesn't provide any
further polygon operations.

We fix these problems by customizing more methods for **ConvexHull**.

.. note:: Obviously, collinear point sets cannot be used to construct polygon.
   So a collinear test should be performed in advance.

.. autoclass:: tadlib.calfea.polygon.Polygon
   :members:

.. autofunction:: tadlib.calfea.polygon.shoelace

.. autofunction:: tadlib.calfea.polygon.collinear


