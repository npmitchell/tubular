function [adIDx, pdIDx] = aux_adjust_dIDx(mesh, cylmesh, t, dpFile,...
    ADBase, PDBase, cylinderMeshCleanDir, ...
    cylinderMeshCleanBase, outadIDxfn, outpdIDxfn, timePoints) 
%[adIDx, pdIDx] = aux_adjust_dIDx(mesh, cylmesh, t, dpFile, ...
%       ADBase, PDBase, cylinderMeshCleanDir, cylinderMeshCleanBase, ...
%       outadIDxfn, outpdIDxfn, timePoints) 
%
% Auxilliary function for adjusting adIDx and pdIDx in
% Generate_Axisymmetric_Pullbacks_Orbifold.m script
% 
% The anterior and posterior "dorsal" points (ie ad and pd) are where the
% cutpath of the cylinderCutMesh starts and ends, respectively.
% 
% Parameters
% ----------
% mesh: cylinder cut mesh with Cleaned Ears (cleanCylCutMesh)
% cylmesh: cylinder cut mesh before ear cleaning
%
% Returns
% -------
% adIDx : anterior dorsal point for cutting the anterior endcap 
% pdIDx : posteriod dorsal point for cutting the posterior endcap 

% Load the AD/PD vertex IDs
disp('Loading ADPD vertex IDs...')
if t == timePoints(1)
    disp(['reading h5 file: ' dpFile])
    adIDx = h5read( dpFile, sprintf( ADBase, t ) );
    pdIDx = h5read( dpFile, sprintf( PDBase, t ) );

    ad3D = cylmesh.v( adIDx, : );
    pd3D = cylmesh.v( pdIDx, : );
else
    currtidx = find(timePoints == t) ;
    prevtp = timePoints(currtidx - 1) ;
    % Load previous mesh and previous adIDx, pdIDx
    prevcylmeshfn = fullfile(cylinderMeshCleanDir, ...
        sprintf( cylinderMeshCleanBase, prevtp )) ;
    prevmesh = read_ply_mod( prevcylmeshfn ); 
    prevadIDx = h5read(outadIDxfn, ['/' sprintf('%06d', prevtp) ]) ;
    % read previous pdIDx with new indices
    prevpdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', prevtp) ]) ;
    ad3D = prevmesh.v(prevadIDx, :) ;
    pd3D = prevmesh.v(prevpdIDx, :) ;
end

trngln = triangulation(mesh.f, mesh.v) ;
boundary = trngln.freeBoundary ;
adIDx = boundary(pointMatch( ad3D, mesh.v(boundary(:, 1), :) ), 1);
pdIDx = boundary(pointMatch( pd3D, mesh.v(boundary(:, 1), :) ), 1);

