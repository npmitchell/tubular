# TubULAR


[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://npmitchell.github.io/tubular/)

> Tube-like sUrface Lagrangian Analysis Resource (TubULAR): for analysis of in-plane and out-of-plane motions of 3D data with dynamic tube-like surfaces of tissues or interfaces.

Cite as: 
> Mitchell, N.P. & Cislo, D.J. "TubULAR: Tracking in toto deformations of dynamic tissues via constrained maps". bioRxiv 2022.04.19.488840 doi:10.1101/2022.04.19.488840


## Overview
![alt text](https://github.com/npmitchell/tubular/blob/gh-pages/docs/_images/fig_tubular_overview_v3.jpg?raw=true)
**TubULAR** transforms 3D volumentric data with a dynamic surface of interest into 
a tube-like, global parameterization based on the intrinsic geometry of the surface 
and then builds measurements of mtion on the dynamic surfaces.

There are four components followed in this package to work with dynamic surfaces: 
1. Obtaining the surfaces (using any method you like, demonstrated here with level sets methods which we incorporated into ImSAnE [1] or using level sets methods from MATLAB built-ins [2] or standalone python modules [3])
2. Visualizing the surfaces (using TexturePatch, included in TubULAR and invoked by TubULAR methods)
3. Parameterizing the entire surfaces over space and time for making measurements (TubULAR methods)
4. Characterizing and separating the different components of this motion (using DEC, referenced by TubULAR methods)
Components 1, 2, and 4 are built on standalone modules in this code package.

We have also incorporated TubULAR into ImSAnE, so that if you already use ImSaNE, you can 
simply update your version and use integralDetector (or morphsnakesDetector) as the detectorType
and tubularMeshWrapper as the fitterType. This yields ImSAnE-style embedding grids, surfaces of
interest, etc.

System Requirements
-------------------
This code is written in MATLAB and has been tested with MATLAB 2018-2021 on Linux Ubuntu 18.04. We expect no issues with Mac or Windows machines, but an installation of GPToolbox and of CGAL are required for full functionality. GPToolbox included along with the license in external/. 
Please also note that upon cloning TubULAR on your local machine, DECLab and TexturePatch, which are linked repositories, will not be pulled automatically with TubULAR. Those repositories are available here: http://github.com/DillonCislo/DEC and https://github.com/npmitchell/TexturePatch. Note that the example scripts expect these repositories to be populated in the folders called DECLab and TexturePatch.

Installation Guide
------------------
Installation instructions for GCAL and GPToolbox are found on the TubULAR documentation:
[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://npmitchell.github.io/tubular/)

The TubULAR codebase itself should download within a few seconds on a "normal" desktop computer.

Demos:
------
For demonstrations:
 example/example_timeseries_gut11Timepoints.m  analyzes a midgut surface, using ImSAnE and morphsnakes to compute surfaces first.
 example/example_zebrafish_heart.m  analyzes a zebrafish heart surface, using morphsnakes to compute surfaces first.
 example/example_singlecoil.m  demonstrates functionality on a synthetic dataset
 
The first demo has both raw data on FigShare: 
and the expected output on FigShare: 
The expected run time for the entire pipeline on a normal computer is a few hours, but many steps are optional. The surface texturepatching in 3D for the first example dataset is prohibitively slow, so it is commented out in the demo script.

Obtaining surfaces:
-------------------
Since different tools will suit different needs, we allow surfaces to be obtained by whatever
method you choose (see for ex, [1]). However, we include a default methods for computing 
these surfaces based on active contour methods. We use build-in MATLAB methods for finding level set surfaces [2]. Note that in some of the data for the publication, we used a similar, but slightly more technical, morphological snakes package in python to find surfaces [3].

Either method can be implemented within the ImSAnE environment [1] as a different detectorType, and TubULAR can utilize an ImSAnE class instance as an input. However, this is not required, since TubULAR also has a getMeshes() method.

If your surfaces aren't coming out well, try adjusting the tension and pressure values of the detectOptions. Typically setting these very nearly to zero is a good starting point.

Surface visualization:
----------------------
Tools for visualizing textured data on the surface are provided in TexturePatch.
The TexturePatch package is integrated into TubULAR so that as you move through a 
TubULAR pipeline, images of your data on the surface in 3D and in mapped 2D spaces are 
drawn on demand.

The TexturePatch package can be used as a standalone package as well.

Coordinate system acquisition:
------------------------------
The surfaces are given global parameterizations for mapping them into the plane and measuring
their dynamics -- both in-plane and out-of-plane. 
In ImSAnE, this step would be called Surface Fitting.

This is done by first ensuring the mesh is a topological cylinder. (If we begin with a 
topological sphere, then two endcaps are cut off based on designation of the two endpoints,
denoted A and P in the code for anterior and posterior.) Then we define coordinates based on 
an orbifold mapping of the surface to the plane, with either endcap being mapped to u=0 and u=1,
and the bulk being mapped to the periodic unit disk. The uv coordiantes of the mapped space provide
a coordinate net for our surface at each timepoint. If the surface is not warping too rapidly,
this uv coordinate mapping will provide a good enough series of images to track objects and follow
tissue motion in 2D. In contrast to [1], this package uses Orbifold-Tutte embeddings [4] for 
conformal maps, which leverage the boundary constraints from the cylindrical topology of the 
surface.

Further refinements of this coordinate system towards a Lagrangian coordinate system follow 
by examining the motion of the coordinate net in 3D space and minimizing the motion of the 
virtual ("pullback") coordinates with respect to their physical motion over time. 
This results in an "s,phi" coordinate system, with s being the integrated 
pathlength along the surface (along the longitudinal coordinate, which cooresponds to the u 
direction at the reference time t0), and phi being the circumferential coordinate (corresponding
to the v direction at the reference time t0). We call this parameterization a "surface Lagrangian"
parameterization, since deformation of the surface in 3D changes the mapping to the plane, but 
cell rearrangements within the surface that do not change the shape of the surface do not affect 
the mapping to the plane.

These "surface Lagrangian" coordinates can then be used to define truly Lagrangian coordinates
by automated analysis of the tissue motion on the surface.

Other pullback coordinates are also supported, as described in the docs. 

Surface velocity measurements:
------------------------------
The DEC package is a self-contained discrete exterior calculus toolkit included here. DEC is
used within TubULAR to compute the divergence and curl of the flow fields on the surface and 
their respective potential fields. Documentation shows how to use DEC -- either with or 
without a TubULAR object -- to compute Laplacians and other DEC metrics of tissue deformation.

Full Documentation
------------------
Click here for the full documentation:
[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://npmitchell.github.io/tubular/)




[1] Heemskirk & Streichan, Nature Methods (2015)

[2] T. Chan and L. Vese, IEEE Transactions on Image Processing 10, 266 (2001).

[3] P. Marquez-Neila, L. Baumela, and L. Alvarez, IEEE Transactions on Pattern Analysis and Machine Intelligence 36, 2 (2014).

[4] Aigerman, N. & Lipman, Y. Orbifold Tutte embeddings. ACM Transactions on Graphics 34, 190:1–190:12 (2015).



