function results = computePCAoverTime(tubi, options)
% results = computePCAoverTime(tubi, options)
% Compute Principal Component Analysis of the sequence of deformations over
% time.
% 
% Parameters
% ----------
% tubi : TubULAR class instance
% options : struct with fields
%   overwrite : bool (default = false)
%       overwrite previous results on disk
%   overwriteImages : bool (default = false)
%       overwrite figures of previous results on disk
%   convert_to_period : bool (default=false)
%       convert all timestamps to period
%   T : numeric (default=NaN)
%       number of timepoints per period if data is periodic in time
%   t0 : int (default=tubi.t0set())
%       timestamp to mark as t=0
%   NmodesToView : int (default= 5)
%       how many modes to plot individually on the reference surface
%       configuration
%   nTimePoints2RmEnds : int (default= 0)
%       how many timepoints to remove from the beginning and end of the
%       timeseries. This is useful if the first and last few datapoints are
%       noisy
%   drawArrowsInPCAPlane : bool (default= false)
%       a plotting option, to draw arrows showing trajectory from 
%       timepoint to timepoint in 2D PCA projective planes
%   drawArrowsInPCA3D : bool (default=false)
%       a plotting option, to draw arrows showing trajectory from 
%       timepoint to timepoint in PCA space
%   meshChoice : cell of strings or single string 
%       (default={'Lagrangian', 'sphi'})
%       the mesh style to use for computing PIV. 
%   axPairs : optional int array (default=[1,2;1,3;2,3])
%       pairs of mode numbers to plot in 2d projective planes (modes are 
%       ordered by rank)
%   plotArrowsOnModes : bool
%       plot arrows on the surfaces for mode decompositions
%   nArrows : int
%       subsampling factor for quiverplot on modes, only used if
%       plotArrowsOnModes == true 
%
% Returns
% -------
% results : cell of structs with fields
%   
%   
% Saves to disk
% -------------
% saves fullfile(tubi.dir.PCAoverTime, meshChoice, sprintf('pcaResults_%s.mat', pcaType)) ;
%   where meshChoice is either 'Lagrangian' or 'sphi'
% saves figures to fullfile(tubi.dir.PCAoverTime, meshChoice)
%
% NPMitchell 2022

%% Default Options
% overwrite previous results if on disk
overwrite = false ;
overwriteImages = false ;
% Convert time units to fractions of a temporal period
convert_to_period = false;
% Period of the periodicity in the data, if there is any. Would adjust
% labels to represent time in units of period
T = NaN ;
% The timepoint at which we set t=0.
t0 = tubi.t0set() ;
% How many of the most significant modes to plot individually (one-by-one)
% on the surface in 3D
NmodesToView = 5 ;
% How many timepoints to cut off the ends, should be set to half-width of
% temporal tripulse filter
nTimePoints2RmEnds = 0 ;
% draw arrows showing trajectory from timepoint to timepoint in 2D PCA
% projective planes
drawArrowsInPCAPlane = false ;
% draw arrows showing trajectory from timepoint to timepoint in PCA space
drawArrowsInPCA3D = false ;
% How to construct the mesh on which we compute the fields to be PCA'd
meshChoice = {'Lagrangian', 'sphi'} ;
% On which fields do we compute PCA
pcaTypes = {'v3d', 'vnVector', 'vt', 'vnScalar', 'divv', '2Hvn', 'gdot'};
% pairs of mode numbers to plot in 2d projections (modes are ordered by rank)
axPairs = [1,2;1,3;2,3] ;
% smoothing parameters for sphi measurement
lambda = tubi.smoothing.lambda ;
lambda_mesh = tubi.smoothing.lambda_mesh ;
lambda_err = tubi.smoothing.lambda_err ;
nmodes = tubi.smoothing.nmodes ;
zwidth = tubi.smoothing.zwidth ;
% Plotting: draw quiverplot arrows on modes
plotArrowsOnModes = false ;
% Plotting: subsampling of quiver plot arrows
nArrows = 250 ;


%% Overwrite default options with supplied options
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
else
    overwriteImages = overwrite ;
end

if isfield(options, 'pcaTypes')
    pcaTypes = options.pcaTypes ;
end

if isfield(options, 'convert_to_period')
    convert_to_period = options.convert_to_period ;
end
if isfield(options, 'T')
    T = options.T ;
end
if isfield(options, 't0')
    t0 = options.t0 ;
end
if isfield(options, 'NmodesToView')
    NmodesToView = options.NmodesToView ;
end
if isfield(options, 'nTimePoints2RmEnds')
    nTimePoints2RmEnds = options.nTimePoints2RmEnds ;
end
if isfield(options, 'drawArrowsInPCAPlane')
    drawArrowsInPCAPlane = options.drawArrowsInPCAPlane ;
end
if isfield(options, 'drawArrowsInPCA3D')
    drawArrowsInPCA3D = options.drawArrowsInPCA3D ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'lambda_err')
    lambda_err = options.lambda_err ;
end
if isfield(options, 'qsub')
    nArrows = options.nArrows ;
end

if isfield(options, 'plotArrowsOnModes')
    plotArrowsOnModes = options.plotArrowsOnModes ;
end


if isfield(options, 'meshChoice')
    meshChoice = options.meshChoice ;
    if isa(meshChoice, 'char')
        meshChoice = {meshChoice} ;
    end
elseif isfield(options, 'meshStyles')
    meshChoice = options.meshStyles ;
    if isa(meshChoice, 'char')
        meshChoice = {meshChoice} ;
    end
end

if isfield(options, 'axPairs')
    axPairs = options.axPairs ;
end


%% Construct PCA for each kind of meshChoice
if nargout > 0
    results = {} ;
end
dmykk = 1 ; % dummy variable for filling in results
for meshChoiceID = 1:numel(meshChoice)
    meshStyle = meshChoice{meshChoiceID} ;
    outputDirRoot = fullfile(tubi.dir.PCAoverTime, meshStyle) ;
    
    if ~exist(outputDirRoot, 'dir')
        mkdir(outputDirRoot)
    end
    
    if strcmpi(meshStyle, 'lagrangian')
        %% Obtain pullback pathlines that define Lagrangian mesh
        t0idx = tubi.xp.tIdx(t0) ;
        pathlines = tubi.getPullbackPathlines(t0, 'vertexPathlines', ...
            'vertexPathlines3D') ;
        % make allV a #vx3x#T array
        % for tidx = 1:numel(tubi.xp.fileMeta.timePoints)
        ntps = length(tubi.xp.fileMeta.timePoints) ;
        % Note permutation to put timepoints last, space indices first
        allVx = reshape(permute(pathlines.vertices3d.vXrs, [2,3,1]), [tubi.nU*tubi.nV, 1, ntps]);
        allVy = reshape(permute(pathlines.vertices3d.vYrs, [2,3,1]), [tubi.nU*tubi.nV, 1, ntps]);
        allVz = reshape(permute(pathlines.vertices3d.vZrs, [2,3,1]), [tubi.nU*tubi.nV, 1, ntps]);
        % Unroll space indices, dimension indices are second, time last
        % allV = cat(2, allVx, allVy, allVz) ;
        
        % Note permutation to put timepoints last, space indices first
        allVxc = allVx(1:tubi.nU*(tubi.nV-1), 1, :) ;
        allVyc = allVy(1:tubi.nU*(tubi.nV-1), 1, :) ;
        allVzc = allVz(1:tubi.nU*(tubi.nV-1), 1, :) ;
        % Unroll space indices, dimension indices are second, time last
        allVc = cat(2, allVxc, allVyc, allVzc) ;

        % face list is static in Lagrangian frame  
        refMeshC = glueRectCylinderCutMeshSeam(pathlines.refMesh) ;
        Fc = refMeshC.f ;

        % allF = repmat(F, [1, 1, ntps]) ;
        allFNc = zeros(size(Fc,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
        for tidx = 1:numel(tubi.xp.fileMeta.timePoints)
            Vc = squeeze(allVc(:,:,tidx)) ;

            % Generate coordinate frame field on each face
            TRc = triangulation(Fc,Vc);
            allFNc(:,:,tidx) = TRc.faceNormal;
        end
        
        if any(contains(pcaTypes, 'vnScalar'))
            % As of now, getPathlineVelocities() does not return normal
            % velocities on vertices, only on advected PIV eval coords.
            % outStruct = measurePathlineVelocities(tubi) ;   
            % So instead do:
            vnScalar = zeros(ntps-1, tubi.nU * tubi.nV) ;
            MKDir = sprintf(tubi.dir.metricKinematics.pathline.measurements, t0) ;
            for tidx = 1:numel(tubi.xp.fileMeta.timePoints)-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('veln_pathline%04d_%06d.mat', t0, tp)) ;
                load(fn, 'veln')
                vnScalar(tidx, :) = veln(:) ;
            end
            vnScalar = vnScalar(:, 1:tubi.nU*(tubi.nV-1)) ;            
        end
        if any(contains(pcaTypes, 'H2vn'))
            H2Vn = zeros(ntps-1, tubi.nU * tubi.nV) ;
            MKDir = sprintf(tubi.dir.metricKinematics.pathline.measurements, t0) ;
            for tidx = 1:numel(tubi.xp.fileMeta.timePoints)-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('H2vn_pathline%04d_%06d.mat', t0, tp)) ;
                load(fn, 'H2vn')
                H2Vn(tidx, :) = H2vn(:) ;
            end            
            H2Vn = H2Vn(:, 1:tubi.nU*(tubi.nV-1)) ;
        end
        if any(contains(pcaTypes, 'divv'))
            divVt = zeros(ntps-1, tubi.nU * tubi.nV) ;
            MKDir = sprintf(tubi.dir.metricKinematics.pathline.measurements, t0) ;
            for tidx = 1:numel(tubi.xp.fileMeta.timePoints)-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('divv_pathline%04d_%06d.mat', t0, tp)) ;
                load(fn, 'divv')
                divVt(tidx, :) = divv(:) ;
            end            
            divVt = divVt(:, 1:tubi.nU*(tubi.nV-1)) ;
        end
        if any(contains(pcaTypes, 'gdot'))
            gdotS = zeros(ntps-1, tubi.nU * tubi.nV) ;
            MKDir = sprintf(tubi.dir.metricKinematics.pathline.measurements, t0) ;
            for tidx = 1:numel(tubi.xp.fileMeta.timePoints)-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('gdot_pathline%04d_%06d.mat', t0, tp)) ;
                load(fn, 'gdot')
                gdotS(tidx, :) = gdot(:) ;
            end            
            gdotS = gdotS(:, 1:tubi.nU*(tubi.nV-1)) ;
        end
        tubi.clearTime() ;
        
    elseif strcmpi(meshStyle, 'sphi')

        %% Load nearly-Lagrangian (s,phi) Meshes ==========================
        close all; clc;

        % time handling
        ntps = numel(tubi.xp.fileMeta.timePoints) ;
        t0idx = tubi.xp.tIdx(t0) ;
        
        % Collate mesh vertices and face normals
        allVc = [];
        allFNc = [];
        for tidx = 1:ntps
            progressbar(tidx, ntps);

            % Load data
            t = tubi.xp.fileMeta.timePoints(tidx);
            tubi.setTime(t); % Set fit time to the current time
            
            % Use closed mesh
            spcutMeshSmRSC = tubi.getCurrentSPCutMeshSmRSC; % The Lagrangian mesh
            Fc = spcutMeshSmRSC.f; % FACE ORDER HAS VERTEX NORMALS POINTING INWARD
            Vc = spcutMeshSmRSC.v; % ASK NOAH ABOUT THIS
            
            % Use open periodic BC mesh
            % spcutMeshSmRS = tubi.getCurrentSPCutMeshSmRS; % The Lagrangian mesh
            % F = spcutMeshSmRS.f; % FACE ORDER HAS VERTEX NORMALS POINTING INWARD
            % V = spcutMeshSmRS.v; % ASK NOAH ABOUT THIS

            % Generate coordinate frame field on each face
            % TR = triangulation(F,V);
            TRc = triangulation(Fc,Vc);
            if isempty(allVc)

                % allV = zeros(size(V,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
                % allFN = zeros(size(F,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
                allVc = zeros(size(Vc,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
                allFNc = zeros(size(Fc,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
            end

            % allV(:,:,tidx) = V;
            % allFN(:,:,tidx) = TR.faceNormal;
            allVc(:,:,tidx) = Vc;
            allFNc(:,:,tidx) = TRc.faceNormal;
            
            % % Check that face normals are the same whether TR or TRc
            % tris = find(any(abs(TRc.faceNormal - TR.faceNormal), 2)) ;
            % trisurf(TR, 'facecolor', 'none'); hold on;
            % trisurf(TRc, 'edgecolor', 'none')
            % COM = barycenter(V, F);
            % plot3(COM(tris, 1), COM(tris, 2), COM(tris, 3), 'ro')
            % diffFN = vecnorm((TRc.faceNormal(tris) - TR.faceNormal(tris)), 2, 2) ./ ...
            %     vecnorm(TRc.faceNormal(tris), 2, 2) ;
            % assert(max(diffFN(:)) < 1e-10)

        end
        
        
        MKDir = fullfile(tubi.dir.metricKinematics.root, ...
            strrep(sprintf('lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d', ...
            lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'), ...
            'measurements');
        if any(contains(pcaTypes, 'vnScalar'))
            % As of now, getPathlineVelocities() does not return normal
            % velocities on vertices, only on advected PIV eval coords.
            % outStruct = measurePathlineVelocities(tubi) ;   
            % So instead do:
            vnScalar = zeros(ntps-1, tubi.nU * tubi.nV) ;
            for tidx = 1:ntps-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('veln_vertices_%06d.mat', tp)) ;
                load(fn, 'veln_filt')
                vnScalar(tidx, :) = veln_filt(:) ;
            end
            vnScalar = vnScalar(:, 1:tubi.nU*(tubi.nV-1)) ;            
        end
        if any(contains(pcaTypes, 'H2vn'))
            H2Vn = zeros(ntps-1, tubi.nU * tubi.nV) ;
            for tidx = 1:ntps-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('H2vn_vertices_%06d.mat', tp)) ;
                load(fn, 'H2vn_filt')
                H2Vn(tidx, :) = H2vn_filt(:) ;
            end            
            H2Vn = H2Vn(:, 1:tubi.nU*(tubi.nV-1)) ;
        end
        if any(contains(pcaTypes, 'divv'))
            divVt = zeros(ntps-1, tubi.nU * tubi.nV) ;
            for tidx = 1:ntps-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('divv_vertices_%06d.mat', tp)) ;
                load(fn, 'divv_filt')
                divVt(tidx, :) = divv_filt(:) ;
            end            
            divVt = divVt(:, 1:tubi.nU*(tubi.nV-1)) ;
        end
        if any(contains(pcaTypes, 'gdot'))
            gdotS = zeros(ntps-1, tubi.nU * tubi.nV) ;
            for tidx = 1:ntps-1
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                tubi.setTime(tp) ;
                fn = fullfile(MKDir, sprintf('gdot_vertices_%06d.mat', tp)) ;
                load(fn, 'gdot_filt')
                gdotS(tidx, :) = gdot_filt(:) ;
            end            
            gdotS = gdotS(:, 1:tubi.nU*(tubi.nV-1)) ;
        end
        tubi.clearTime() ;
    else
        error(['Could not recognize meshStyle = ' meshStyle ...
            '. Try again with options.meshChoice=Lagrangian or sphi'])
    end

    % View Mesh ---------------------------------------------------------------
    % close all;
    % 
    % % Choose the time point to view --> t0idx
    % V = allV(:,:,t0idx);
    % fn = allFN(:,:,t0idx);
    % 
    % COM = barycenter(V, F);
    % ssf = false( size(F,1), 1 );
    % sst = false( size(F,1), 1 );
    % sst(1:4:end) = true;
    % % sst = true( size(P, 1), 1 );
    % 
    % trisurf(triangulation(F,V), ...
    %     'FaceColor', 0.8 * [1 1 1], 'EdgeColor', 'none');
    % 
    % hold on
    % 
    % quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
    %     COM([ssf; ssf; sst]), ...
    %     fn([sst; ssf; ssf]), fn([ssf; sst; ssf]), ...
    %     fn([ssf; ssf; sst]), ...
    %     1, 'LineWidth', 2, 'Color', 'k' );
    % 
    % hold off
    % 
    % axis equal

    clear spcutMeshSmRSC spcutMeshSmRS
    clear tidx V u TR 

    %% Visualize Velocity PCA Modes ===========================================
    close all; clc;

    % The full face-based 3D velocity field (#t x #F x 3)
    fprintf('First computing face-based velocity field over time... ');
    v3d = zeros(ntps-1, size(Fc,1), 3) ;    
    for tidx = 1:numel(tubi.xp.fileMeta.timePoints) -1 
        V0c = squeeze(allVc(:,:,tidx)) ; % vertices at t0
        V1c = squeeze(allVc(:,:,tidx+1)) ; % vertices at t1

        % velocities on vertices (Lagrangian frame) in units 
        % of tubi.spaceUnits / [tubi.timeUnits * (timepoint1-timepoint0)]
        vel = V1c - V0c ;

        % Push the velocities to the face barycenters
        [V2Fc, F2Vc] = meshAveragingOperators(Fc, V0c) ;
        velF = V2Fc * vel ;

        % Collimate velocities on barycenters
        v3d(tidx, :, :) = velF ;
    end

    %% Get subsampling for quiver plots
    if plotArrowsOnModes
        disp('Load/compute sampling for vector field')
        fn = fullfile(tubi.dir.PCAoverTime, 'quiver_subsampling.mat') ;
        if exist(fn, 'file')
            load(fn, 'sampleIDx_faces') 
        else
            % get subsampling
            COM = barycenter(V0c, Fc);
            % find farthest vertices
            sampleIDx_vertices = farthestPointVertexSampling(nArrows, V0c, Fc, ...
                [1], struct('preview', true)) ;
            % pointmatch those vertices to face barycenter sampling
            sampleIDx_faces = pointMatch(V0c(sampleIDx_vertices,:), COM) ;
            save(fn, 'sampleIDx_faces', 'sampleIDx_vertices')
        end
        ssf = false( size(Fc,1), 1 );
        sst = false( size(Fc,1), 1 );
        sst(sampleIDx_faces) = true;
    end
    
    % THIS IS ONLY APPROXIMATE SINCE (s,phi) is not fully Lagrangian
    % v3d = tubi.getVelocityAverage().vf;
    fprintf('Done\n');

    %% Load velocities and Perform PCA -----------------------------------------
    fprintf('Performing PCA... ');
    tic
    % The normal component of the velocity field
    FN1 = permute(allFNc, [3 1 2]); 
    FN1(end, :, :) = [];
    vn = sum( FN1 .* v3d, 3 ) .* FN1;

    % The tangential component of each velocity field
    vt = v3d - vn;

    % Consider doing PCA on each kind of velocity field
    for pcaTypeID = 1:length(pcaTypes)
        pcaType = pcaTypes{pcaTypeID} ;
        fnResult = fullfile(outputDirRoot, sprintf('pcaResults_%s.mat', pcaType)) ;
        outputDirPCA = fullfile(outputDirRoot, pcaType) ;
        if ~exist(outputDirPCA, 'dir')
            mkdir(outputDirPCA)
        end
        
        % load from disk or recompute
        if exist(fnResult, 'file') && ~overwrite
            disp(['Found result on disk. Loading: ' fnResult])
            if nargout > 0
                load(fnResult, 'pcaModeError', 'coeff', ...
                'pcaType', 'convert_to_period', 'T', 't0', ...
                'NmodesToView', 'nTimePoints2RmEnds', 'meshChoice') ;
                results{dmykk} = struct('pcaModeError', pcaModeError, ...
                    'coeff', coeff, ...
                    'pcaType', pcaType, ...
                    'convert_to_period', convert_to_period, ...
                    'T', T, 't0', t0, ...
                    'NmodesToView', NmodesToView, ...
                    'nTimePoints2RmEnds', nTimePoints2RmEnds, ...
                    'meshChoice', meshChoice) ;
                dmykk = dmykk + 1 ; 
            end
        else
            disp('Did not find result on disk...')
            pcaType = pcaTypes{pcaTypeID} ;

            % Choose which component of the velocity field to analyze
            if strcmpi(pcaType, 'v3d')
                v_pca = v3d;
                titleString0 = '3D velocity $$\bf v$$';
                saveStrPrefix = 'PCA_Full_Velocity';
            elseif strcmpi(pcaType, 'vnScalar')
                v_pca = vnScalar;
                titleString0 = 'normal velocity $$v_n$$';
                saveStrPrefix = 'PCA_Normal_Velocity_Scalar';
            elseif strcmpi(pcaType, 'vnVector')
                v_pca = vn ;
                titleString0 = 'normal velocity, $$v_n \hat{n}$$';
                saveStrPrefix = 'PCA_Normal_Velocity_Vector';
            elseif strcmpi(pcaType, 'vt')
                v_pca = vt;
                titleString0 = ['tangential velocity ' ...
                    '$${\bf v}_\parallel$$'];
                saveStrPrefix = 'PCA_Tangent_Velocity';
            elseif strcmpi(pcaType, 'H2vn')
                v_pca = H2Vn;
                titleString0 = 'out-of-plane deformation $$2H v_n$$';
                saveStrPrefix = 'PCA_H2vn';
            elseif strcmpi(pcaType, 'divv')
                v_pca = divVt;
                titleString0 = 'in-plane divergence $$\nabla \cdot {\bf v}_\parallel$$';
                saveStrPrefix = 'PCA_divv';
            elseif strcmpi(pcaType, 'gdot')
                v_pca = gdotS;
                titleString0 = 'area rate of change $$\mathrm{Tr}[g^{-1} \dot{g}] / 2$$';
                saveStrPrefix = 'PCA_gdot';
            else
                error('Invalid velocity type for PCA analysis');
            end

            if convert_to_period
                v_pca = v_pca * T;
            end

            % Cat into a large array
            if numel(size(v_pca)) == 3
                % Note we convert v_pca into a #timepoints x 3*#faces array
                v_pca = [ v_pca(:,:,1), v_pca(:,:,2), v_pca(:,:,3) ];
            else
                disp('v_pca is 2D, meaning it is a scalar field.')
            end
            
            % Clip ends to account for poor temporal smoothing
            if nTimePoints2RmEnds > 0
                v_pca(1:nTimePoints2RmEnds, :) = [];
                v_pca(end-nTimePoints2RmEnds+1:end, :) = [];
            end
            if convert_to_period
                v_pca = v_pca * T; 
            end
            
            % Perform PCA on the Cartesian components of the velocity
            coeff = pca(v_pca);

            toc
            fprintf(' Done\n');

            %% Generate Visualization -------------------------------------

            % Extract Lagrangian mesh
            V0c = allVc(:,:,t0idx);
            DEC = DiscreteExteriorCalculus(Fc, V0c);
            [V2Fc, F2Vc] = meshAveragingOperators(Fc, V0c);
            VNc = per_vertex_normals(V0c, Fc, 'Weighting', 'angle');
            
            
            % Chose the mode that will be visualized
            displayTypes = {'normal', 'displacement', 'rotational', 'dilatational'} ;
            for dispID = 1:length(displayTypes)

                plotType = displayTypes{dispID} ;
                for modeID = 1:NmodesToView

                    if (modeID == 0)
                        error('Broken handling');
                        vmode = 0 .* coeff(:,1);
                    else
                        vmode = coeff(:, modeID);
                    end

                    % Extract visualization parameters
                    [~, ~, ~, xyzlim] = tubi.getXYZLims();
                    xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
                    ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
                    zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;

                    % String for naming this mode figure
                    saveStrMode = sprintf([saveStrPrefix, '_mode_%d'], modeID);
                    
                    clearvars divV rotV harmV scalarP vectorP VM
                    
                    % Treat the visualization differently if Vector/Scalar
                    isVectorField = numel(vmode) == size(Fc,1)*3 ;
                    if isVectorField
                        vmode = reshape(vmode, size(Fc,1), 3);
                        vmode_F = vmode;
                        vmode = F2Vc * vmode;
                        % Take components of this 3D vector field
                        [divV, rotV, harmV, scalarP, vectorP] = ...
                            DEC.helmholtzHodgeDecomposition(vmode_F);
                        
                        FNc = faceNormal(triangulation(Fc, V0c));
                        normV = dot(harmV, FNc, 2) .* FNc;
                        harmV = harmV - normV;

                        % Scale PCA compnents to produce visible displacements
                        vscale = 400;
                        VMc = V0c + vscale * vmode;

                        saveStr = [saveStrMode, '_', plotType, '.png'];
                        outputFigFn = fullfile(outputDirPCA, saveStr);
                        if ~exist(outputFigFn, 'file') || overwriteImages

                            % fig = figure('Visible', 'off',  'units', 'centimeters', ...
                            %     'position', [0,0,xwidth,ywidth]) ;
                            fig = figure('Visible', 'off',  'units', 'centimeters') ;
                            for viewAngleID = 1:4
                                subplot(2, 2, viewAngleID)

                                quiver_lw = 0.75;
                                if strcmpi(plotType, 'displacement')

                                    cmap = brewermap(256, '*RdBu');
                                    crange = max(sqrt(sum(vmode.^2, 2))) * [-1 1];
                                    vcolors = mapValueToColor( ...
                                        sqrt(sum(vmode.^2, 2)) .* sign(dot(vmode, VNc, 2)), ...
                                        crange, cmap);

                                    patch( 'Faces', Fc, 'Vertices', VMc, 'FaceVertexCData', vcolors, ...
                                        'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                        'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                        'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                        'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                    titleString = [titleString0 ' - ' sprintf('mode %d', modeID)];

                                    if convert_to_period
                                        labelString = ...
                                            ['sgn$$[\vec{\bf V} \cdot \hat{n}_0] ||\vec{\bf V}||$$' ...
                                            ' [' tubi.spaceUnits '/T]'];
                                    else
                                        labelString = ...
                                            ['sgn$$[\vec{\bf V} \cdot \hat{n}_0] ||\vec{\bf V}||$$' ...
                                            ' [' tubi.spaceUnits '/' tubi.timeUnits ']'];
                                    end

                                elseif strcmpi(plotType, 'rotational')

                                    vectorPV = F2Vc * vectorP;
                                    cmap = brewermap(256, '*RdBu');
                                    crange = max(abs(vectorPV)) * [-1 1];

                                    plotRotV = normalizerow(rotV);

                                    patch( 'Faces', Fc, 'Vertices', V0c, ...
                                        'FaceVertexCData', vectorPV, ...
                                        'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                        'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                        'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                        'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                    hold on

                                    if plotArrowsOnModes
                                        quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                            COM([ssf; ssf; sst]), ...
                                            plotRotV([sst; ssf; ssf]), plotRotV([ssf; sst; ssf]), ...
                                            plotRotV([ssf; ssf; sst]), ...
                                            1, 'LineWidth', quiver_lw, 'Color', 'k' );

                                        hold off
                                    end

                                    titleString = [titleString0 ' - ' ...
                                        sprintf('rotational part, mode %d', modeID)];

                                elseif strcmpi(plotType, 'dilatational')

                                    cmap = brewermap(256, '*RdBu');
                                    crange = max(abs(scalarP)) * [-1 1];

                                    plotDivV = normalizerow(divV);

                                    patch( 'Faces', Fc, 'Vertices', V0c, ...
                                        'FaceVertexCData', scalarP, ...
                                        'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                        'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                        'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                        'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                    if plotArrowsOnModes
                                        hold on
                                        quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                            COM([ssf; ssf; sst]), ...
                                            plotDivV([sst; ssf; ssf]), plotDivV([ssf; sst; ssf]), ...
                                            plotDivV([ssf; ssf; sst]), ...
                                            1, 'LineWidth', quiver_lw, 'Color', 'k' );
                                        hold off
                                    end
                                    
                                    titleString = [titleString0 ' - ' ...
                                        sprintf('dilatational part, mode %d', modeID)];

                                elseif strcmpi(plotType, 'harmonic')

                                    cmap = brewermap(512, '*RdBu');
                                    cmap = cmap((1:256).', :);

                                    absHVP = sqrt(sum((F2Vc * harmV).^2, 2));
                                    crange = 0.01 * max(absHVP) * [0 1];

                                    plotHarmV = normalizerow(harmV);

                                    patch( 'Faces', Fc, 'Vertices', V0c, ...
                                        'FaceVertexCData', absHVP, ...
                                        'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                        'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                        'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                        'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
                                    
                                    if plotArrowsOnModes
                                        hold on
                                        quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                            COM([ssf; ssf; sst]), ...
                                            plotHarmV([sst; ssf; ssf]), plotHarmV([ssf; sst; ssf]), ...
                                            plotHarmV([ssf; ssf; sst]), ...
                                            1, 'LineWidth', quiver_lw, 'Color', 'k' );

                                        hold off
                                    end

                                    titleString = [titleString0 ' - ' ...
                                        sprintf('harmonic part, mode %d', modeID)];

                                    labelString = '$$||V_{harm}||$$';

                                elseif strcmpi(plotType, 'normal')

                                    cmap = brewermap(256, '*RdBu');
                                    crange = max(sqrt(sum((F2Vc*normV).^2, 2))) * [-1 1];
                                    vcolors = mapValueToColor( ...
                                        sqrt(sum((F2Vc*normV).^2, 2)) .* sign(dot(F2Vc * normV, VNc, 2)), ...
                                        crange, cmap);

                                    plotNormV = normalizerow(normV);

                                    patch( 'Faces', Fc, 'Vertices', V0c, ...
                                        'FaceVertexCData', vcolors, ...
                                        'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                        'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                        'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                        'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                    if plotArrowsOnModes
                                        hold on;
                                        quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                            COM([ssf; ssf; sst]), ...
                                            plotNormV([sst; ssf; ssf]), plotNormV([ssf; sst; ssf]), ...
                                            plotNormV([ssf; ssf; sst]), ...
                                            1, 'LineWidth', quiver_lw, 'Color', 'k' );

                                        hold off
                                    end

                                    titleString = [titleString0 ' - ' ...
                                        sprintf('normal part, mode %d', modeID)];
                                else
                                    error('Invalid plot type');
                                end

                                axis equal

                                scaleFactor = 1.1;
                                xlim(scaleFactor * [xmin, xmax]);
                                ylim(scaleFactor * [ymin, ymax]);
                                zlim(scaleFactor * [zmin, zmax]);

                                % camlight
                                xlabel(['AP Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
                                ylabel(['Lateral Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
                                zlabel(['DV Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');

                                colormap(cmap);
                                set(gca, 'Clim', crange);

                                % cb = colorbar;
                                % cb.Label.String = labelString;
                                % cb.Label.Interpreter = 'latex';

                                set(gcf, 'Color', [1 1 1]);
                                axis off

                                if viewAngleID == 1
                                    if contains(titleString, '$')
                                        sgtitle(titleString, 'Interpreter', 'latex');
                                    else
                                        sgtitle(titleString, 'FontWeight', 'normal');
                                    end
                                end

                                set(gca, 'FontSize', 6);
                                set(gca, 'FontWeight', 'normal');

                                view([0,(viewAngleID-1)*90])
                            end

                            % Resize Figure for Paper -------------------------------------------------
                            set(fig, 'Units', 'centimeters');

                            % ratio = fig.Position(4) ./ fig.Position(3);
                            % fig.Position(3) = 4;
                            % fig.Position(4) = ratio * fig.Position(3);

                            set(fig, 'PaperPositionMode', 'auto');
                            set(fig, 'position', [0, 0, 16, 10])
                            set(fig.Children, 'FontName', 'Helvetica');

                            % Save figure
                            disp(['Saving figure: ' outputFigFn])
                            export_fig(outputFigFn, '-png', '-r500');
                            close all
                        end
                    else
                        disp('PCA field is scalar field')
                        
                        % Plot displacement of vertices from reference mesh
                        if strcmpi(pcaType, 'vnScalar') && strcmpi(plotType, 'displacement')
                            % Scale PCA compnents to produce visible displacements
                            vscale = 400;
                            
                            VMc = V0c + vscale * vmode .* VNc;
                            saveStr = [saveStrMode, '_', plotType, '.png'];
                            outputFigFn = fullfile(outputDirPCA, saveStr);
                            if ~exist(outputFigFn, 'file') || overwriteImages
                                
                                fig = figure('Visible', 'off',  'units', 'centimeters') ;
                                for viewAngleID = 1:4
                                    subplot(2, 2, viewAngleID)
                               
                                    cmap = brewermap(256, '*RdBu');
                                    crange = max(sqrt(sum(vmode.^2, 2))) * [-1 1];
                                    vcolors = mapValueToColor( vmode, ...
                                        crange, cmap);

                                    patch( 'Faces', Fc, 'Vertices', VMc, 'FaceVertexCData', vcolors, ...
                                        'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                        'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                        'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                        'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                    titleString = [titleString0 ' - ' sprintf('mode %d', modeID)];

                                    if convert_to_period
                                        labelString = ...
                                            ['v_n' ...
                                            ' [' tubi.spaceUnits '/T]'];
                                    else
                                        labelString = ...
                                            ['V_n' ...
                                            ' [' tubi.spaceUnits '/' tubi.timeUnits ']'];
                                    end
                                    axis equal

                                    scaleFactor = 1.1;
                                    xlim(scaleFactor * [xmin, xmax]);
                                    ylim(scaleFactor * [ymin, ymax]);
                                    zlim(scaleFactor * [zmin, zmax]);

                                    xlabel(['AP Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
                                    ylabel(['Lateral Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
                                    zlabel(['DV Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');

                                    colormap(cmap);
                                    set(gca, 'Clim', crange);

                                    % cb = colorbar;
                                    % cb.Label.String = labelString;
                                    % cb.Label.Interpreter = 'latex';

                                    set(gcf, 'Color', [1 1 1]);
                                    axis off

                                    if viewAngleID == 1
                                        if contains(titleString, '$')
                                            sgtitle(titleString, 'Interpreter', 'latex');
                                        else
                                            sgtitle(titleString, 'FontWeight', 'normal');
                                        end
                                    end

                                    set(gca, 'FontSize', 6);
                                    set(gca, 'FontWeight', 'normal');

                                    view([0,(viewAngleID-1)*90])
                                end

                                % Resize Figure for Paper -----------------
                                set(fig, 'Units', 'centimeters');
                                set(fig, 'PaperPositionMode', 'auto');
                                set(fig, 'position', [0, 0, 16, 10])
                                set(fig.Children, 'FontName', 'Helvetica');

                                % Save figure
                                disp(['Saving figure: ' outputFigFn])
                                export_fig(outputFigFn, '-png', '-r500');
                                close all
                            end
                            
                        elseif strcmpi(plotType, 'normal')
                            
                            % Save image of scalar field alone
                            saveStr = [saveStrMode '.png'];
                            outputFigFn = fullfile(outputDirPCA, saveStr);
                            if ~exist(outputFigFn, 'file') || overwriteImages

                                % fig = figure('Visible', 'off',  'units', 'centimeters', ...
                                %     'position', [0,0,xwidth,ywidth]) ;
                                fig = figure('Visible', 'off',  'units', 'centimeters') ;
                                for viewAngleID = 1:4
                                    subplot(2, 2, viewAngleID)
                               
                                    cmap = brewermap(256, '*RdBu');
                                    crange = max(sqrt(sum(vmode.^2, 2))) * [-1 1];
                                    vcolors = mapValueToColor( vmode, ...
                                        crange, cmap);

                                    patch( 'Faces', Fc, 'Vertices', V0c, 'FaceVertexCData', vcolors, ...
                                        'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                        'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                        'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                        'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                    titleString = [titleString0 ' - ' sprintf('mode %d', modeID)];

                                    scaleFactor = 1.1;
                                    xlim(scaleFactor * [xmin, xmax]);
                                    ylim(scaleFactor * [ymin, ymax]);
                                    zlim(scaleFactor * [zmin, zmax]);
                                    
                                    xlabel(['AP Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
                                    ylabel(['Lateral Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
                                    zlabel(['DV Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');

                                    colormap(cmap);
                                    set(gca, 'Clim', crange);
                                    set(gcf, 'Color', [1 1 1]);
                                    axis off

                                    if viewAngleID == 1
                                        if contains(titleString, '$')
                                            sgtitle(titleString, 'Interpreter', 'Latex');
                                        else
                                            sgtitle(titleString, 'FontWeight', 'normal');
                                        end
                                    end
                                    set(gca, 'FontSize', 6);
                                    set(gca, 'FontWeight', 'normal');

                                    view([0,(viewAngleID-1)*90])
                                end

                                % Resize Figure for Paper -----------------
                                set(fig, 'Units', 'centimeters');
                                set(fig, 'PaperPositionMode', 'auto');
                                set(fig, 'position', [0, 0, 16, 10])
                                set(fig.Children, 'FontName', 'Helvetica');

                                % Save figure
                                disp(['Saving figure: ' outputFigFn])
                                export_fig(outputFigFn, '-png', '-r500');
                                close all
                            end
                        end
                    end
                end
            end

            % clear FN1 vn vt pcaType v_pca 
            % clear V0 V2F F2V VN DEC modeID titleString
            % clear vmode VM lapVMode xyzlim xmin ymin zmin xmax ymax zmax
            % clear cmap vrange vcolors scaleFactor

            %% Generate 2D PCA Mode Comparison Plot ===================================
            close all; clc;

            % Generate Visualization --------------------------------------------------
            % color by timestamp in proper units (tubi.timeInterval)
            tps = tubi.xp.fileMeta.timePoints(1:end-1) - tubi.t0;
            % Remove timepoints that we removed from the PCA
            tps(1:nTimePoints2RmEnds) = []; 
            tps((end-nTimePoints2RmEnds+1):end) = [];
            
            
            % cmap = brewermap(512, '*RdBu'); cmap = cmap((1:256).', :);
            % cmap = twilight_shifted(512); 
            % cmap = cmap((1:256).', :);
            cmap = lapaz ;
            % cmap = plasma(); 
            % cmap = brewermap(256, 'OrRd');

            if convert_to_period
                trange = [0 max(tps)]/T * tubi.timeInterval;
                timeColors = mapValueToColor(tps.'/T * tubi.timeInterval, ...
                    trange, cmap);
            else
                trange = [0 max(tps)] * tubi.timeInterval;
                timeColors = mapValueToColor(tps.' * tubi.timeInterval, ...
                    trange, cmap);
            end

            
            % Choose the modes that will be visualized
            for axPairID = 1:length(axPairs)
                axID1 = axPairs(axPairID, 1) ;
                axID2 = axPairs(axPairID, 2) ;

                saveStrPair = sprintf([saveStrPrefix '_m%dm%d.pdf'], axID1, axID2);
                outputFigFn = fullfile(outputDirPCA, saveStrPair);
                if ~exist(outputFigFn, 'file') || overwriteImages
                    % The projection of the velocities along the component axes
                    v1 = dot(v_pca.', repmat(coeff(:, axID1), 1, size(v_pca,1)), 1).';
                    v2 = dot(v_pca.', repmat(coeff(:, axID2), 1, size(v_pca,1)), 1).';

                    % Scaled separation vectors (used to visualize flow field);
                    vScale = 30;
                    if convert_to_period
                        vScale = vScale * T; 
                    end
                    vvec = [v1, v2]; vvec = vvec(2:end, :)-vvec(1:(size(vvec,1)-1), :);
                    vvec = vScale .* vvec ./ sqrt(sum(vvec.^2, 2));

                    % Clip the ends of the velocity vectors to account for poor temporal
                    % smoothing
                    % v1(1:2) = []; v1((end-1):end) = [];
                    % v2(1:2) = []; v2((end-1):end) = [];
                    % vvec(1:2, :) = []; vvec(end,:) = [];

                    % Axis limits
                    fig = figure('Visible', 'off',  'units', 'centimeters') ;
                    hold on

                    if drawArrowsInPCAPlane
                        quiver_lw = 0.75;
                        quiver(v1(1:size(vvec,1)), v2(1:size(vvec,1)), ...
                            vvec(:,1), vvec(:,2), 0, 'Color', 'k', 'LineWidth', quiver_lw );
                    end

                    scatter(v1, v2, [], timeColors, 'filled');
                    hold off
                    axis equal

                    % xlim(xLim); ylim(yLim);
                    axis square

                    if convert_to_period
                        xlabel(['mode ' num2str(axID1) ' [' tubi.spaceUnits '/T]']);
                        ylabel(['mode ' num2str(axID2) ' [' tubi.spaceUnits '/T]']);

                    else
                        xlabel(['principal component ' num2str(axID1) ...
                            ' (' tubi.spaceUnits '/' tubi.timeUnits ')'], 'interpreter', 'latex');
                        ylabel(['principal component ' num2str(axID2) ...
                            ' (' tubi.spaceUnits '/' tubi.timeUnits ')'], 'interpreter', 'latex');
                    end

                    cb = colorbar;
                    set(gca, 'Clim', trange);
                    colormap(cmap);
                    if convert_to_period
                        cb.Label.String = 'time [T]';
                    else
                        cb.Label.String = ['time [' tubi.timeUnits ']'];
                    end
                    % cb.Label.Interpreter = 'latex';
                    % cb.Ticks = [0.2 1 2];
                    cb.Position(1) = 0.8;
                    cb.Position(2) = 0.25;
                    cb.Position(3) = 0.045;
                    cb.Position(4) = 0.6;
                    cb.Label.Position(1) = 2.0;
                    set(gcf, 'Color', [1 1 1]);
                    if contains(titleString0, '$')
                        sgtitle(titleString0, 'interpreter', 'latex') 
                    else
                        sgtitle(titleString0) 
                    end
                    
                    % Resize Figure for Paper -----------------------------
                    set(fig, 'Units', 'centimeters');
                    fig.Position(3:4) = [16,10] ;
                    set(gca, 'FontWeight', 'normal');
                    set(fig, 'PaperPositionMode', 'auto');
                    set(fig.Children, 'FontName', 'Helvetica');
                    export_fig(outputFigFn, '-pdf', '-r200');
                end
            end

            %% Generate 3D PCA Mode Comparison Plot ===================================
            % Compares the 3 most important modes produced by PCA
            close all; clc;


            % Generate Visualization --------------------------------------------------
            outputFigFn = fullfile(outputDirPCA, [saveStrPrefix '_3D_trajectory.png']);
            if ~exist(outputFigFn, 'file') || overwriteImages 
                % The projection of the velocities along the component axes
                v1 = dot(v_pca.', repmat(coeff(:, 1), 1, size(v_pca,1)), 1).';
                v2 = dot(v_pca.', repmat(coeff(:, 2), 1, size(v_pca,1)), 1).';
                v3 = dot(v_pca.', repmat(coeff(:, 3), 1, size(v_pca,1)), 1).';

                % Scaled separation vectors (used to visualize flow field);
                vScale = 45;
                if convert_to_period, vScale = vScale * T; end
                vvec = [v1, v2, v3]; vvec = vvec(2:end, :)-vvec(1:(size(vvec,1)-1), :);
                vvec = vScale .* vvec ./ sqrt(sum(vvec.^2, 2));

                % Clip the ends of the velocity vectors to account for poor temporal
                % smoothing
                % v1(1:2) = []; v1((end-1):end) = [];
                % v2(1:2) = []; v2((end-1):end) = [];
                % v3(1:2) = []; v3((end-1):end) = [];
                % vvec(1:2, :) = []; vvec(end, :) = []; 

                % Axis limits
                scaleFactor = 1.5;
                xLim = scaleFactor * [-1 1] * max(abs(v1)); % (max(v1)-min(v1))/2;
                yLim = scaleFactor * [-1 1] * max(abs(v2)); % (max(v2)-min(v2))/2;
                zLim = scaleFactor * [-1 1] * max(abs(v3)); % (max(v3)-min(v3))/2;
                % if convert_to_period
                %     xLim = T * xLim;
                %     yLim = T * yLim;
                %     zLim = T * zLim;
                % end

                [x,y,z] = sphere(15);
                % x = x(:); y = y(:); z = z(:);
                r = 0.025 * max(max(max(xLim), max(yLim)), max(zLim)) ;
                if convert_to_period, r = T * r; end

                quiver_lw = 0.75;
                fig = figure('Visible', 'off',  'units', 'centimeters') ;
                hold on

                if  drawArrowsInPCA3D 
                    quiver3(v1(1:size(vvec,1)), v2(1:size(vvec,1)), v3(1:size(vvec,1)), ...
                        vvec(:,1), vvec(:,2), vvec(:,3), ...
                        0, 'Color', 'k', 'LineWidth', quiver_lw );
                end

                for tidx = 1:numel(v1)

                    surf(v1(tidx) + r.*x, v2(tidx) + r.*y, v3(tidx) + r.*z, ...
                        'FaceColor', timeColors(tidx,:), 'EdgeColor', 'none', ...
                        'FaceLighting', 'gouraud');

                    surf(v1(tidx) + r.*x, v2(tidx) + r.*y, zLim(1) + 0.*z, ...
                        'FaceColor', 0.6 * [1 1 1], 'EdgeColor', 'none', ...
                        'FaceLighting', 'none');

                    surf(v1(tidx) + r.*x, yLim(2) + 0.*y, v3(tidx) + r.*z, ...
                        'FaceColor', 0.6 * [1 1 1], 'EdgeColor', 'none', ...
                        'FaceLighting', 'none');

                    surf(xLim(2) + 0.*x, v2(tidx) + r.*y, v3(tidx) + r.*z, ...
                        'FaceColor', 0.6 * [1 1 1], 'EdgeColor', 'none', ...
                        'FaceLighting', 'none');

                end

                if drawArrowsInPCA3D 
                    quiver3(v1(1:size(vvec,1)), v2(1:size(vvec,1)), zLim(1) * ones(numel(v3)-1,1), ...
                        vvec(:,1), vvec(:,2), 0*vvec(:,3), ...
                        0, 'Color', 0.6 * [1 1 1], 'LineWidth', quiver_lw );

                    quiver3(v1(1:size(vvec,1)), yLim(2) * ones(numel(v2)-1,1), v3(1:size(vvec,1)), ...
                        vvec(:,1), 0*vvec(:,2), vvec(:,3), ...
                        0, 'Color', 0.6 * [1 1 1], 'LineWidth', quiver_lw );

                    quiver3(xLim(2) * ones(numel(v1)-1,1), v2(1:size(vvec,1)), v3(1:size(vvec,1)), ...
                        0*vvec(:,1), vvec(:,2), vvec(:,3), ...
                        0, 'Color', 0.6 * [1 1 1], 'LineWidth', quiver_lw );
                end

                hold off
                axis equal
                view([-45 30]);
                xlim(xLim); ylim(yLim); zlim(zLim);

                camlight
                if convert_to_period
                    xlabel(['mode ' num2str(1)]);
                    ylabel(['mode ' num2str(2)]);
                    zlabel(['mode ' num2str(3)]);
                else
                    xlabel(['mode ' num2str(1) ...
                        ' [' tubi.spaceUnits '/' tubi.timeUnits ']'], 'interpreter', 'latex');
                    ylabel(['mode ' num2str(2) ...
                        ' [' tubi.spaceUnits '/' tubi.timeUnits ']'], 'interpreter', 'latex');
                    zlabel(['mode ' num2str(3) ...
                        ' [' tubi.spaceUnits '/' tubi.timeUnits ']'], 'interpreter', 'latex');
                end

                % Get this code from the MathWorks 
                % axesLabelsAlign3D();

                grid off
                box on

                set(gcf, 'Color', [1 1 1]);

                cb = colorbar;
                cb.Location = 'eastoutside';
                set(gca, 'Clim', trange);
                colormap(cmap);
                if convert_to_period
                    cb.Label.String = 'time [T]';
                else
                    cb.Label.String = ['Time [' tubi.timeUnits ']'];
                end
                % cb.Label.Interpreter = 'latex';
                % cb.Ticks = [0.2 1 2];
                cb.Position(1) = 0.85;
                cb.Position(2) = 0.4;
                cb.Position(4) = 0.4;
                cb.Label.Position(1) = 2.5;

                axPos = get(gca, 'Position');
                axRatio = axPos(4) / axPos(3);
                axPos(1) = 0.175;
                axPos(2) = 0.25;
                axPos(3) = 0.6;
                axPos(4) = axRatio * axPos(3);
                set(gca, 'Position', axPos);

                if contains(titleString0, '$')
                    sgtitle(titleString0, 'interpreter', 'latex');
                else
                    sgtitle(titleString0);
                end
                
                set(gcf, 'Color', [1 1 1]);

                % Resize Figure for Paper -------------------------------------------------
                set(fig, 'Units', 'centimeters');
                % ratio = fig.Position(4) ./ fig.Position(3);
                % fig.Position(3) = 5.5;
                % fig.Position(4) = ratio * fig.Position(3);

                fig.Position(3:4) = [16, 10] ;
                % set(gca, 'FontSize', 5);
                % set(gca, 'FontWeight', 'normal');

                set(fig, 'PaperPositionMode', 'auto');
                set(fig.Children, 'FontName', 'Helvetica');

                % export_fig(saveStr, '-pdf', '-r200');
                export_fig(outputFigFn, '-png', '-r500');
            end
            clear x y z r

            %% Determine the Error Associated with Ignoring PCA Modes =================
            close all; clc;

            % Perform PCA on the Cartesian components of the velocity
            % [coeff,score,latent,tsquared,explained,mu] = pca(v_pca);

            % This lets you get all the modes but is ridiculously slow
            % coeff = pca(v_pca, 'Economy', false, 'NumComponents', 50);


            % Determine the Error Associated with Removing Each Mode Individually -----

            pcaModeError = zeros(size(coeff,2), 1);
            for i = 1:numel(pcaModeError)
                if mod(i, 20) == 0
                    disp(['Checking for contribution from mode ' num2str(i)])
                end

                % The projection of the velocities along the current component axis
                vi = dot(v_pca.', repmat(coeff(:, i), 1, size(v_pca,1)), 1) .* ...
                    repmat(coeff(:, i), 1, size(v_pca,1));

                curErr = sum(sum((vi.').^2, 2) ./ sum(v_pca.^2, 2), 1) ./ size(v_pca,1);

                pcaModeError(i) = curErr;

            end

            % Generate Visualization --------------------------------------------------

            plot_lw = 0.75;
            scatter_size = 2;

            fig = figure('Visible', 'off',  'units', 'centimeters') ;

            stem(1:numel(pcaModeError), pcaModeError.', ...
                'LineStyle', '-', ...
                'MarkerFaceColor', tubi.plotting.colors(1,:), ...
                'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
                'LineWidth', plot_lw, ...
                'MarkerSize', scatter_size);

            xlim([1 numel(pcaModeError)]);

            xlabel('mode number');
            ylabel('mode contribution');

            set(gcf, 'Color', [1 1 1]);

            % Resize Figure for Paper -------------------------------------------------
            set(fig, 'Units', 'centimeters');

            ratio = fig.Position(4) ./ fig.Position(3);
            fig.Position(3) = 3.5;
            fig.Position(4) = ratio * fig.Position(3);

            set(gca, 'FontSize', 5);
            set(gca, 'FontWeight', 'normal');

            set(fig, 'PaperPositionMode', 'auto');
            set(fig.Children, 'FontName', 'Helvetica');

            saveStr = ['Mode_Contribution_' pcaType '.pdf'];
            saveStr = fullfile(outputDirPCA, saveStr);
            export_fig(saveStr, '-pdf', '-r200');

            %% Save everything to disk
            disp(['Saving results to disk: ' fnResult])
            save(fnResult, 'pcaModeError', 'coeff', ...
                'pcaType', 'convert_to_period', 'T', 't0', ...
                'NmodesToView', 'nTimePoints2RmEnds', 'meshChoice') ;
            if nargout > 0
                results{dmykk} = struct('pcaModeError', pcaModeError, ...
                    'coeff', coeff, ...
                    'pcaType', pcaType, ...
                    'convert_to_period', convert_to_period, ...
                    'T', T, 't0', t0, ...
                    'NmodesToView', NmodesToView, ...
                    'nTimePoints2RmEnds', nTimePoints2RmEnds, ...
                    'meshChoice', meshChoice) ;
                dmykk = dmykk + 1 ; 
            end
        end
    end
end
