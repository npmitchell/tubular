function [length_lobes, area_lobes, volume_lobes] = measureLobeDynamics(QS, options)
%MEASURELOBEDYNAMICS(QS, options)
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   The current QuapSlap object whose lobes to plot
% options : struct with fields
%   overwrite : bool, default=false
%       overwrite saved lobes on disk if previously computed
%   save_ims : bool, default=true
%       save images of the lobes evolving over time
% 
% Returns
% -------
% length_lobes : #timepoints x 1 float
% area_lobes : #timepoints x 1 float 
% volume_lobes : #timepoints x 1 float
%
% NPMitchell 2020

%% Default options
t0 = QS.t0set() ;

%% Unpack options
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'save_ims') 
    save_ims = options.save_ims ;
else
    save_ims = true ;
end

%% Load or compute lobe dynamics

lobe_dynamics_fn = QS.fileName.lobeDynamics ;
if exist(lobe_dynamics_fn, 'file') && ~overwrite
    % Load length, surface area, and volume dynamics for each lobe
    disp('Loading lobe length, area, and volume...')
    load(lobe_dynamics_fn, 'length_lobes', 'area_lobes', 'volume_lobes')
else
    % Unpack QS
    timePoints = QS.xp.fileMeta.timePoints ;
    lobeDir = QS.dir.lobe ;
    spcutMeshBase = QS.fullFileBase.spcutMesh ;
    nU = QS.nU ;
    nV = QS.nV ;
    [rot, trans] = QS.getRotTrans() ;
    resolution = QS.APDV.resolution ;
    [~, ~, ~, xyzlim] = QS.getXYZLims() ;
    
    % Define colors for plotting
    colors = QS.plotting.colors ;
    
    % Load fold positions
    load(QS.fileName.fold, 'folds', 'ssfold', 'ssmax')    
    
    [length_lobes, area_lobes, volume_lobes] = ...
        aux_compute_lobe_dynamics(folds, ssfold, ssmax, lobeDir, ...
                timePoints, t0, QS.timeInterval, QS.timeUnits, QS.spaceUnits, ...
                spcutMeshBase, nV, nU, rot, trans, ...
                resolution, QS.flipy, ...
                xyzlim, colors, save_ims, overwrite) ;
            
            
    % % For debugging -----------------------------------------------------
    % tp = 200 ;
    % load(sprintf(spcutMeshBase, tp), 'spcutMesh') ;
    % % Load the centerline too
    % % fn = sprintf(clineDVhoopBase, tp) ;
    % % load(fn, 'avgpts')
    % avgpts = spcutMesh.avgpts ;
    % mcline = spcutMesh.mcline ;
    % xyz = QS.xyz2APDV(spcutMesh.v) ;
    % figure; plot3(avgpts(:, 1), avgpts(:, 2), avgpts(:, 3))
    % hold on;
    % plot3(spcutMesh.mcline(:, 1), spcutMesh.mcline(:, 2), spcutMesh.mcline(:, 3))
    % trisurf(triangulation(spcutMesh.f, xyz), 'edgecolor', 'none', 'faceAlpha', 0.1)
    % %--------------------------------------------------------------------
    
    % Save surface area and volume dynamics for each lobe
    save(lobe_dynamics_fn, 'length_lobes', 'area_lobes', 'volume_lobes')
end
clearvars fold_ofn

