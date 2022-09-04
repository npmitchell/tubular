Common issues in TubULAR and suggestions
=======


- Data loads but only one slice is present
> Try toggling the swapZT parameter in the xp metadata.

- Data loads but each timepoint is loaded as multiple timestamps
> Try toggling the swapZT parameter in the xp metadata, since the Z dimension is likely interpreted as time T.

- The surfaces aren't coming out well.
> Try adjusting the tension and pressure values of the detectOptions. Typically setting these very nearly to zero is a good starting point.
> If that fails, check whether your initial guess for the level set is inside the surface you want to find. The initial guess could be a small sphere that is well within the surface you want to find. Subsequent timepoints are found by first using the previous timepoint's output as an initial guess. To shrink the previous timepoint's output before evolving it, adjust pre_pressure.

- Fast marching centerline extraction fails and/or the image output shows the path leaving the outside of the tubular surface. 
> First erase the existing (bad) centerline files from disk (txt files) and try re-running with options.preview = true. If the visual output doesn't tell you what's wrong right away, try reducing the subsampling factor (``options.res``) by which the volume is reduced for fast marching. You may also/alternatively need to adjust options.normal_step in tubi.alignAPDVCoords(), which places the A and P points connected by the centerline a distance options.normal_step from the surface.


- Fast marching centerline extraction takes too long to calculate.
> Increase the subsampling factor (``options.res``). These crude fast-marching centerlines will be superceded by centerlines computed from the conformal mapping to the plane later anyway, so there is no point wasting too much time on them. Their only purpose is to topologically constrain the conformal map minimizing Dirichlet energy (creaing uv/sphi coordinates), so as long as they are curves that live inside the mesh and span the system, they should work fine. 

- ``Undefined function 'smooth' for input arguments of type 'double'``
> Install Curve Fitting Toolbox or obtain a working version of the ``smooth`` function.

- ``Undefined function 'perform_front_propagation_3d' for input arguments of type 'double'.``
> This is a function that is inside gptoolbox, in the external folder. Make sure you run GPToolbox's external/toolbox_fast_marching/compile_mex.m successfully, run with MATLAB from within the parent directory (ie the current working directory should be something like ``tubular/external/gptoolbox/external/toolbox_fast_marching/``).  