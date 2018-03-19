# Numerical steepest descent solver

Currently quite slow and unstable. Please let me know if you find ways to break it. Designed to solve integrals of the form:

![equation](http://latex.codecogs.com/gif.latex?%5Cint_%5Cgamma%20f%28z%29%5Cmathrm%7Be%7D%5E%7B%5Cmathrm%7Bi%7D%5Comega%20g%28z%29%7D%5Cmathrm%7Bd%7Dz)

*Infinite contours*

To compute an integral over an infinite path ![equation](http://latex.codecogs.com/gif.latex?%5Cgamma), from a valley with complex angle ``a`` to angle ``b``, type
```
[z, w] = NSD45( a, b, freq, Npts, G, 'analytic', true, 'settleRad', R, 'ainf', 'binf');
```
to obtain weighs ``w`` and nodes ``z``. The option ``'analytic', true`` assumes that the phase is analytic. ``Npts`` is proportional number of points weights and nodes, typically no more than ``Npts=15`` is sufficient. ``freq>0`` is the frequency parameter. ``R>0`` represents the radius of a ball, outside of which the SD can be assumed to be approximately straight lines, all the way to infinity. The array ``G`` should contain (vectorised) function handles corresponding to derivatives of the phase ![equation](http://latex.codecogs.com/gif.latex?g), that is

![equation](http://latex.codecogs.com/gif.latex?G%3A%3D%28g%2Cg%27%2Cg%27%27%2C%5Cldots%2Cg%5E%7B%28n%29%7D%29%2C)

where ![equation](http://latex.codecogs.com/gif.latex?n) is chosen to be at least one larger than the highest order stationary point of ![equation](http://latex.codecogs.com/gif.latex?g)
