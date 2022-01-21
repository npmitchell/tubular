function generateSPCutMeshStack(QS, spcutMeshStackOptions)
% generateSPCutMeshSmStack(QS, spcutMeshStackOptions)
%   Create TIFF stacks of normally evolved smoothed meshes in (s, phi)
%   coordinate system.
%   Same as generateSPCutMeshSmStack() but for non-smoothed (in time) 
%   pullback images. 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% spcutMeshStackOptions: struct with fields
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
if isfield(spcutMeshStackOptions, 'overwrite')
    overwrite = spcutMeshStackOptions.overwrite ;
end
if isfield(spcutMeshStackOptions, 'n_outward')
    n_outward = spcutMeshStackOptions.n_outward ;
end
if isfield(spcutMeshStackOptions, 'n_inward')
    n_inward = spcutMeshStackOptions.n_inward ;
end
if isfield(spcutMeshStackOptions, 'layer_spacing')
    layer_spacing = spcutMeshStackOptions.layer_spacing ;
end
if isfield(spcutMeshStackOptions, 'smoothIter')
    smoothIter = spcutMeshStackOptions.smoothIter ;
end
if isfield(spcutMeshStackOptions, 'preSmoothIter')
    preSmoothIter = spcutMeshStackOptions.preSmoothIter ;
end
if isfield(spcutMeshStackOptions, 'imSize')
    imSize = spcutMeshStackOptions.imSize;
end

% Unpack QS
spcutMeshBase = QS.fullFileBase.spcutMesh ;
fileNameBase = QS.fileBase.name ;

for qq = 1:length(QS.xp.fileMeta.timePoints)
    tt = QS.xp.fileMeta.timePoints(qq) ;
    disp(['t = ' num2str(tt)])
    QS.setTime(tt)
    
    % Load time-smoothed mesh
    load(sprintf(spcutMeshBase, tt), 'spcutMesh') ;
    
    %--------------------------------------------------------------
    % Generate Output Image File
    %--------------------------------------------------------------
    spacingstr = strrep(sprintf('%0.2fum', layer_spacing * QS.APDV.resolution), '.', 'p') ;
    if contains(fileNameBase, '%')
        imfn_spsm = sprintf( fullfile( QS.dir.im_re_stack, ...
            [fileNameBase, '_%02d_%02d_' spacingstr '.tif']), tt, n_outward, n_inward ) ;
    else
        disp('No timestamp in image name --> must be fixed sample')
        try 
            assert(length(QS.xp.fileMeta.timePoints) == 1)
        catch
            disp('No timestamp in image name, but more than one timepoint!')
        end
        imfn_spsm = sprintf( fullfile( QS.dir.im_re_stack, ...
            [fileNameBase, '_%02d_%02d_' spacingstr '.tif']), n_outward, n_inward ) ;
    end
    
    if ~exist(imfn_spsm, 'file') || overwrite
        % Load 3D data for coloring mesh pullback
        % Note that axisOrder is applying upon invoking getCurrentData()
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
        spcutMesh.u = spcutMesh.sphi ;
        spcutMesh = rmfield(spcutMesh, 'uv') ;
        spcutMesh = rmfield(spcutMesh, 'sphi') ;
        tmp = spcutMesh.u ;
        tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
        arspsm = minimizeIsoarealAffineEnergy( spcutMesh.f, spcutMesh.v, tmp );
                
        % aux_generate_orbifold( spcutMeshSm, QS.a_fixed * 0.5, IV, imfn_spsm, Options)
        aux_generate_orbifold( spcutMesh, arspsm * 0.5, IV, imfn_spsm, Options)
        % error('stopping here')
    end
    clear Options
end
