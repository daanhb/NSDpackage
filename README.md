**Numerical steepest descent solver**

Currently quite slow and unstable. Please let me know if you find ways to break it.

To compute an integral over an infinite contour, from a valley with complex angle ``a`` to angle ``b``, type
```
[z, w] = NSD45( a, b, freq, Npts, G, 'analytic', true, 'settleRad', R, 'ainf', 'binf');
```
to obtain weighs ``w`` and nodes ``z``. The option ``'analytic', true`` assumes that the phase is analytic. ``Npts`` is proportional number of points weights and nodes, typically no more than ``Npts=15`` is sufficient. ``freq>0`` is the frequency parameter. ``R>0`` represents the radius of a ball, outside of which the SD can be assumed to be approximately straight lines, all the way to infinity.
