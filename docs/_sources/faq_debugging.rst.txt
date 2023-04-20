Common issues in TubULAR and suggestions
========================================

- Almost every method of TubULAR has an ``overwrite`` option (which is false be default for any given method), and most have a ``preview`` option, defined in an optional struct that can be passed to the method. If you overwrite something upstream of other calculations, there are potential scenarios where this could lead to issues. For example, if you change the endpoints which are used to cut off the endcaps of the mesh but don't rerun ``tubi.alignMeshesAPDV()``, then the aligned meshes will not reflect your changes.
> Set ``options.overwrite=true`` and pass ``options`` to any relevant method that is downstream of what you overwrote. 

- Data loads but only one slice is present
> Try toggling the swapZT parameter in the xp metadata.

- Data loads but each timepoint is loaded as multiple timestamps
> Try toggling the swapZT parameter in the xp metadata, since the Z dimension is likely interpreted as time T.

- The surfaces aren't coming out well.
> Try adjusting the tension and pressure values of the detectOptions. Typically setting these very nearly to zero is a good starting point.
> If that fails, check whether your initial guess for the level set is inside the surface you want to find. The initial guess could be a small sphere that is well within the surface you want to find. Subsequent timepoints are found by first using the previous timepoint's output as an initial guess. To shrink the previous timepoint's output before evolving it, adjust pre_pressure.

- Cutting off the endcaps fails for a particular timepoint.
> This usually means something wrong with the mesh itself. You can inspect that mesh in MeshLab to look for issues. If nothing appears unusual, consider remaking the mesh for that one timepoint, maybe with the tiniest adjustment of parameters or a little more iLastik training just to make sure you get a slightly different mesh. You can also redo the endcap removal with slightly larger or slightly smaller distance thresholds to the endcap vertices. If all else fails, try a different way of identifying which part of the mesh to cut off. The default method is to cut off the largest contiguous mesh patch that is within a threshold distance from the endcap vertex. 

- Fast marching centerline extraction fails and/or the image output shows the path leaving the outside of the tubular surface. 
> First erase the existing (bad) centerline files from disk (txt files) and try re-running with options.preview = true. If the visual output doesn't tell you what's wrong right away, try reducing the subsampling factor (``options.res``) by which the volume is reduced for fast marching. You may also/alternatively need to adjust options.normal_step in tubi.alignAPDVCoords(), which places the A and P points connected by the centerline a distance options.normal_step from the surface.

- Fast marching centerline extraction takes too long to calculate.
> Increase the subsampling factor (``options.res``). These crude fast-marching centerlines will be superceded by centerlines computed from the conformal mapping to the plane later anyway, so there is no point wasting too much time on them. Their only purpose is to topologically constrain the conformal map minimizing Dirichlet energy (creaing uv/sphi coordinates), so as long as they are curves that live inside the mesh and span the system, they should work fine. 

- ``Undefined function 'smooth' for input arguments of type 'double'``
> Install Curve Fitting Toolbox or obtain a working version of the ``smooth`` function.

- ``Undefined function 'perform_front_propagation_3d' for input arguments of type 'double'.``
> This is a function that is inside gptoolbox, in the external folder. Make sure you run gptoolbox's external/toolbox_fast_marching/compile_mex.m successfully, run with MATLAB from within the parent directory (ie the current working directory should be something like ``tubular/external/gptoolbox/external/toolbox_fast_marching/``).  

- ``Could NOT find Matlab`` error while compiling gptoolbox
> As described on the gptoolbox `mex troubleshooting page <https://github.com/alecjacobson/gptoolbox/tree/master/mex>`_, cmake requires a hardcoded version mapping between MATLAB's numerical versions (e.g. 9.5) and their named versions (e.g. 2018b). Unfortunately this mapping goes stale a couple of times a year when new MATLAB versions are released. If you find that cmake can't figure out your MATLAB version, and you have a fairly recent version of MATLAB, it's likely that you need to update the mapping. To do so, update the ``MATLAB_VERSIONS_MAPPING`` variable in ``tubular/external/gptoolbox/mex/cmake/FindMATLAB.cmake`` and add your version number. You can check your MATLAB version using the ``version`` command. If you still run into an error like

.. code-block:: console

        CMake Error at ... (message):
        Could NOT find Matlab (missing: Matlab_INCLUDE_DIRS Matlab_MEX_LIBRARY
        Matlab_MEX_EXTENSION Matlab_ROOT_DIR MEX_COMPILER MX_LIBRARY ENG_LIBRARY)
        (found version "NOTFOUND")

Then you might consider hardcoding your local MATLAB directory into the cmake instructions. Navigate to the file ``tubular/external/gptoolbox/mex/CMakeLists.txt`` and add the modify the appropriate section to read

.. code-block:: bash

        # Find matlab
        if(MATLAB_PROXY)
          set(Matlab_ROOT_DIR "${GPTOOLBOX_MEX_ROOT}/external/matlab")
          gptoolbox_download_matlab()
        endif()

        set(Matlab_ROOT_DIR "your/local/MATLAB/R?????")
        find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)

where you should replace ``your/local/MATLAB/R?????`` with the correct path to your local MATLAB directory.
