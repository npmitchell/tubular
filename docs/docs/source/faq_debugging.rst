Common issues in TubULAR and suggestions
=======


- Data loads but only one slice is present
> Try toggling the swapZT parameter in the xp metadata.

- Data loads but each timepoint is loaded as multiple timestamps
> Try toggling the swapZT parameter in the xp metadata, since the Z dimension is likely interpreted as time T.

- Fast marching centerline extraction fails and/or the image output shows the path leaving the outside of the tubular surface. 
> First erase the existing (bad) centerline files from disk (txt files) and try re-running with options.preview = true. If the visual output doesn't tell you what's wrong right away, try reducing the subsampling factor (`options.res') by which the volume is reduced for fast marching. You may also/alternatively need to adjust options.normal_step in tubi.alignAPDVCoords(), which places the A and P points connected by the centerline a distance options.normal_step from the surface.


- Fast marching centerline extraction takes too long to calculate.
> Increase the subsampling factor (`options.res'). These crude fast-marching centerlines will be superceded by centerlines computed from the conformal mapping to the plane later anyway, so there is no point wasting too much time on them. Their only purpose is to topologically constrain the conformal map minimizing Dirichlet energy (creaing uv/sphi coordinates), so as long as they are curves that live inside the mesh and span the system, they should work fine. 
