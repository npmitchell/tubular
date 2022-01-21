function generateSPCutMeshSmStack(QS, spcutMeshSmStackOptions)
% generateSPCutMeshSmStack(QS, spcutMeshSmStackOptions)
%   Create TIFF stacks of normally evolved smoothed meshes in (s, phi)
%   coordinate system.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% spcutMeshSmStackOptions: struct with fields
%   - n_outward : int 
%       number of steps in positive normal direction to sample
%   - n_inward : int
%       number of steps in negative normal direction to sample
%   - overwrite : bool
%       whether to overwrite spcutMeshSmStack on disk
%
% Returns 
% -------
%
% NPMitchell 2020

% Unpack options
n_outward = 10 ;
n_inward = 10 ;
layer_spacing = 1 ;
smoothIter = 0 ;
preSmoothIter = 0 ;
overwrite = false ;
imSize = 500 ;
if isfield(spcutMeshSmStackOptions, 'overwrite')
    overwrite = spcutMeshSmStackOptions.overwrite ;
end
if isfield(spcutMeshSmStackOptions, 'n_outward')
    n_outward = spcutMeshSmStackOptions.n_outward ;
end
if isfield(spcutMeshSmStackOptions, 'n_inward')
    n_inward = spcutMeshSmStackOptions.n_inward ;
end
if isfield(spcutMeshSmStackOptions, 'layer_spacing')
    layer_spacing = spcutMeshSmStackOptions.layer_spacing ;
end
if isfield(spcutMeshSmStackOptions, 'smoothIter')
    smoothIter = spcutMeshSmStackOptions.smoothIter ;
end
if isfield(spcutMeshSmStackOptions, 'preSmoothIter')
    preSmoothIter = spcutMeshSmStackOptions.preSmoothIter ;
end
if isfield(spcutMeshSmStackOptions, 'imSize')
    imSize = spcutMeshSmStackOptions.imSize;
end

% Unpack QS
spcutMeshSmBase = QS.fullFileBase.spcutMeshSm ;
fileNameBase = QS.fileBase.name ;

for qq = 1:length(QS.xp.fileMeta.timePoints)
    tt = QS.xp.fileMeta.timePoints(qq) ;
    disp(['t = ' num2str(tt)])
    QS.setTime(tt)
    
    % Load time-smoothed mesh
    load(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;
    
    %--------------------------------------------------------------
    % Generate Output Image File
    %--------------------------------------------------------------
    spacingstr = strrep(sprintf('%0.2fum', layer_spacing * QS.APDV.resolution), '.', 'p') ;
    imfn_spsm = sprintf( fullfile( QS.dir.im_r_sme_stack, ...
        [fileNameBase, '_%02d_%02d_' spacingstr '.tif']), tt, n_outward, n_inward ) ;
    
    if ~exist(imfn_spsm, 'file') || overwrite
        % Load 3D data for coloring mesh pullback
        QS.getCurrentData() ;
        IV = QS.currentData.IV ;
        
        fprintf(['Generating SP output image for sm mesh: ' imfn_spsm]);
        % Assigning field spcutMesh.u to be [s, phi] (ringpath
        % and azimuthal angle)
        Options.numLayers = [n_outward, n_inward] ;
        Options.layerSpacing = layer_spacing ;
        Options.smoothIter = smoothIter ;
        Options.preSmoothIter = preSmoothIter ;
        Options.yLim = [-0.5, 1.5] ;
        Options.imSize = imSize ;
        
        % Note that we pass a_fixed * 0.5 since the image is extended by a
        % factor of two
        % Compute relaxed aspect ratio
        tmp = spcutMeshSm.u ;
        tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
        arspsm = minimizeIsoarealAffineEnergy( spcutMeshSm.f, spcutMeshSm.v, tmp );
                
        % aux_generate_orbifold( spcutMeshSm, QS.a_fixed * 0.5, IV, imfn_spsm, Options)
        aux_generate_orbifold( spcutMeshSm, arspsm * 0.5, IV, imfn_spsm, Options)
        % error('stopping here')
    end
    clear Options
end
