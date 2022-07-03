function results = computeLBSoverTime(tubi, options)
% results = computeLBSoverTime(tubi, options)
% Compute Laplace-Beltrami spectral decomposition of the sequence of
% deformations over time.
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
%   decompT : int (defualt = t0)
%       The time for which the mesh Laplacian eigenbasis will be
%       constructed
%   NLBSModes : int (default = 200)
%       The number of Laplacien eigenvectors to calculate
%   NmodesToView : int (default= 5)
%       how many modes to plot individually on the reference surface
%       configuration
%   nTimePoints2RmEnds : int (default= 0)
%       how many timepoints to remove from the beginning and end of the
%       timeseries. This is useful if the first and last few datapoints are
%       noisy
%   drawArrowsInLBSPlane : bool (default= false)
%       a plotting option, to draw arrows showing trajectory from 
%       timepoint to timepoint in 2D LBS projective planes
%   drawArrowsInLBS3D : bool (default=false)
%       a plotting option, to draw arrows showing trajectory from 
%       timepoint to timepoint in LBS space
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
% saves fullfile(tubi.dir.LBSoverTime, meshChoice, sprintf('lbsResults_%s.mat', lbsType)) ;
%   where meshChoice is either 'Lagrangian' or 'sphi'
% saves figures to fullfile(tubi.dir.LBSoverTime, meshChoice)
%
% Dillon Cislo + NPMitchell 2022

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
% The decomposition time point
decompT = tubi.t0set() ;
% The number of Laplacian eigenvectors to calculate
NLBSModes = 200;
% How many of the most significant modes to plot individually (one-by-one)
% on the surface in 3D
NmodesToView = 5 ;
% How many timepoints to cut off the ends, should be set to half-width of
% temporal tripulse filter
nTimePoints2RmEnds = 0 ;
% draw arrows showing trajectory from timepoint to timepoint in 2D LBS
% projective planes
drawArrowsInLBSPlane = false ;
% draw arrows showing trajectory from timepoint to timepoint in LBS space
drawArrowsInLBS3D = false ;
% How to construct the mesh on which we compute the fields to be LBS'd
meshChoice = {'Lagrangian', 'sphi'} ;
% On which fields do we compute LBS
lbsTypes = {'v3d', 'vnVector', 'vt', 'vnScalar', 'divv', '2Hvn', 'gdot'};
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

if isfield(options, 'lbsTypes')
    lbsTypes = options.lbsTypes ;
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
if isfield(options, 'decompT')
    decompT = options.decompT ;
end
if isfield(options, 'NLBSModes')
    NLBSModes = options.NLBSModes;
end
if isfield(options, 'NmodesToView')
    NmodesToView = options.NmodesToView ;
end
if isfield(options, 'nTimePoints2RmEnds')
    nTimePoints2RmEnds = options.nTimePoints2RmEnds ;
end
if isfield(options, 'drawArrowsInLBSPlane')
    drawArrowsInLBSPlane = options.drawArrowsInLBSPlane ;
end
if isfield(options, 'drawArrowsInLBS3D')
    drawArrowsInLBS3D = options.drawArrowsInLBS3D ;
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

decompTIDx = tubi.xp.tIdx(decompT) ;

%% Construct LBS for each kind of meshChoice
if nargout > 0
    results = {} ;
end
dmykk = 1 ; % dummy variable for filling in results
for meshChoiceID = 1:numel(meshChoice)
    meshStyle = meshChoice{meshChoiceID} ;
    outputDirRoot = fullfile(tubi.dir.LBSoverTime, meshStyle) ;
    
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
        allVNc = zeros(size(allVc,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
        for tidx = 1:numel(tubi.xp.fileMeta.timePoints)
            Vc = squeeze(allVc(:,:,tidx)) ;

            % Generate coordinate frame field on each face
            TRc = triangulation(Fc,Vc);
            allFNc(:,:,tidx) = TRc.faceNormal;
            
            % Generate vertex normals
            [~, F2Vc] = meshAveragingOperators(Fc, Vc);
            VNc = F2Vc * TRc.faceNormal;
            VNc = VNc ./ repmat(sqrt(sum(VNc.^2, 2)), 1, 3);
            allVNc(:,:,tidx) = VNc;
            
        end
        
        if any(contains(lbsTypes, 'vnScalar'))
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
        if any(contains(lbsTypes, 'H2vn'))
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
        if any(contains(lbsTypes, 'divv'))
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
        if any(contains(lbsTypes, 'gdot'))
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
            VNc = spcutMeshSmRSC.vn;
            
            % Use open periodic BC mesh
            % spcutMeshSmRS = tubi.getCurrentSPCutMeshSmRS; % The Lagrangian mesh
            % F = spcutMeshSmRS.f; % FACE ORDER HAS VERTEX NORMALS POINTING INWARD
            % V = spcutMeshSmRS.v; % ASK NOAH ABOUT THIS
            % Vn = spcutMeshSmRS.vn;

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
            % allVN(:,:,tidx) = VN;
            allVc(:,:,tidx) = Vc;
            allFNc(:,:,tidx) = TRc.faceNormal;
            allVNc(:,:,tidx) = VNc;
            
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
        if any(contains(lbsTypes, 'vnScalar'))
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
        if any(contains(lbsTypes, 'H2vn'))
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
        if any(contains(lbsTypes, 'divv'))
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
        if any(contains(lbsTypes, 'gdot'))
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
    
    % Construct Mesh Laplacian Eigenbasis ---------------------------------
    
    lbsV = allVc(:,:,decompTIDx);
    lbsF = Fc;
    lbsDEC = DiscreteExteriorCalculus(lbsF, lbsV);
    
    % The (unweighted, symmetric, cotangent) Laplace-Beltrami operator
    % Equivalent to 'cotmatrix.m' in gptoolbox
    LBO = lbsDEC.dd1 * lbsDEC.hd1 * lbsDEC.d0;
    
    % The Laplacian eigenbasis. Returns Nmodes eigenvalues that are closest
    % to 0 in decreasing order of value
    % LBVec = #V X #Nmodes to view matrix whose columns are the orthonormal
    % eigen vectors listed
    % LBVal = #Nmodes X #Nmodes diagonal matrix whose diagonal elements are
    % the corresponding eigenvalues.
    [LBVec, LBVal] = eigs(LBO, NLBSModes, 0);
    
    % Generate Mode Visualizations ----------------------------------------
    % NOTE: This section generates visualizations of the mesh
    % Laplacian eigenbasis used to generate the decompositions. It
    % is solely a property of the mesh and does NOT depend on the
    % signal that was decomposed.
    
    % Extract Lagrangian mesh
    V0c = lbsV; % allVc(:,:,t0idx);
    DEC = lbsDEC; % DiscreteExteriorCalculus(Fc, V0c);
    [V2Fc, F2Vc] = meshAveragingOperators(Fc, V0c);
    VNc = per_vertex_normals(V0c, Fc, 'Weighting', 'angle');
    
    % Chose the mode that will be visualized
    displayTypes = {'signal', 'displacement'} ;
    saveStrPrefix = 'Laplacian_Eigenbasis';
    outputDirLBS = fullfile(outputDirRoot, 'Laplacian_Eigenbasis') ;
    if ~exist(outputDirLBS, 'dir')
        mkdir(outputDirLBS)
    end
    for dispID = 1:length(displayTypes)
        
        plotType = displayTypes{dispID} ;
        for modeID = 1:NmodesToView
            
            if (modeID == 0)
                error('Broken handling');
            else
                vmode = LBVec(:, modeID);
            end
            
            % Extract visualization parameters
            [~, ~, ~, xyzlim] = tubi.getXYZLims();
            xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
            ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
            zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;
            
            % String for naming this mode figure
            saveStrMode = sprintf([saveStrPrefix, '_mode_%d'], modeID);
            
            % Plot displacement of vertices from reference mesh
            if strcmpi(plotType, 'displacement')
                
                % Scale LBS compnents to produce visible displacements
                vscale = 400;
                
                VMc = V0c + vscale * vmode .* VNc;
                saveStr = [saveStrMode, '_', plotType, '.png'];
                outputFigFn = fullfile(outputDirLBS, saveStr);
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
                        
                        titleString = ['Laplacian eigenbasis - ' sprintf('mode %d', modeID)];
                        labelString = ['Basis coefficients'];
                        
                        % if convert_to_period
                        %     labelString = ...
                        %         ['v_n' ...
                        %         ' [' tubi.spaceUnits '/T]'];
                        % else
                        %     labelString = ...
                        %         ['V_n' ...
                        %         ' [' tubi.spaceUnits '/' tubi.timeUnits ']'];
                        % end
                        
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
                    
                    % Resize Figure for Paper -----------------------------
                    set(fig, 'Units', 'centimeters');
                    set(fig, 'PaperPositionMode', 'auto');
                    set(fig, 'position', [0, 0, 16, 10])
                    set(fig.Children, 'FontName', 'Helvetica');
                    
                    % Save figure
                    disp(['Saving figure: ' outputFigFn])
                    export_fig(outputFigFn, '-png', '-r500');
                    close all
                end
                
            elseif strcmpi(plotType, 'signal')
                
                % Save image of scalar field alone
                saveStr = [saveStrMode '.png'];
                outputFigFn = fullfile(outputDirLBS, saveStr);
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
                        
                        titleString = ['Laplacian eigenbasis - ' sprintf('mode %d', modeID)];
                        labelString = ['Basis coefficients'];
                        
                        % if convert_to_period
                        %     labelString = ...
                        %         ['v_n' ...
                        %         ' [' tubi.spaceUnits '/T]'];
                        % else
                        %     labelString = ...
                        %         ['V_n' ...
                        %         ' [' tubi.spaceUnits '/' tubi.timeUnits ']'];
                        % end
                        
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

                        % cb = colorbar;
                        % cb.Label.String = labelString;
                        % cb.Label.Interpreter = 'latex';
                        
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
    
    % clear FN1 vn vt pcaType v_pca
    % clear V0 V2F F2V VN DEC modeID titleString
    % clear vmode VM lapVMode xyzlim xmin ymin zmin xmax ymax zmax
    % clear cmap vrange vcolors scaleFactor


    %% Visualize Velocity LBS Modes ===========================================
    close all; clc;

    % The full vertex-based 3D velocity field (#V x 3 x (#t-1))
    fprintf('First computing vertex-based velocity field over time... ');
    v3d = zeros(size(lbsV,1), 3, ntps-1) ;    
    for tidx = 1:numel(tubi.xp.fileMeta.timePoints)-1 
        V0c = squeeze(allVc(:,:,tidx)) ; % vertices at t0
        V1c = squeeze(allVc(:,:,tidx+1)) ; % vertices at t1

        % velocities on vertices (Lagrangian frame) in units 
        % of tubi.spaceUnits / [tubi.timeUnits * (timepoint1-timepoint0)]
        v3d(:, :, tidx) = V1c - V0c ;

    end

    %% Get subsampling for quiver plots
    % V0c = allVc(:,:,1);
    if plotArrowsOnModes
        disp('Load/compute sampling for vector field')
        fn = fullfile(tubi.dir.LBSoverTime, 'quiver_subsampling.mat') ;
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

    %% Load velocities and Perform LBS -----------------------------------------
    fprintf('Performing LBS... ');
    tic
    % The normal component of the velocity field (#V x 3 x #t)
    vn = sum(allVNc(:,:,1:(end-1)) .* v3d, 2) .* allVNc(:,:,1:(end-1));

    % The tangential component of each velocity field
    vt = v3d - vn;

    % Consider doing LBS on each kind of velocity field
    for lbsTypeID = 1:length(lbsTypes)
        lbsType = lbsTypes{lbsTypeID} ;
        fnResult = fullfile(outputDirRoot, sprintf('lbsResults_%s.mat', lbsType)) ;
        outputDirLBS = fullfile(outputDirRoot, lbsType) ;
        if ~exist(outputDirLBS, 'dir')
            mkdir(outputDirLBS)
        end
        
        % load from disk or recompute
        if exist(fnResult, 'file') && ~overwrite
            disp(['Found result on disk. Loading: ' fnResult])
            if nargout > 0
                load(fnResult, 'lbsModeError', 'coeff', ...
                'lbsType', 'convert_to_period', 'T', 't0', ...
                'NmodesToView', 'nTimePoints2RmEnds', 'meshChoice') ;
                results{dmykk} = struct('lbsModeError', lbsModeError, ...
                    'coeff', coeff, ...
                    'lbsType', lbsType, ...
                    'convert_to_period', convert_to_period, ...
                    'T', T, 't0', t0, ...
                    'NmodesToView', NmodesToView, ...
                    'nTimePoints2RmEnds', nTimePoints2RmEnds, ...
                    'meshChoice', meshChoice) ;
                dmykk = dmykk + 1 ; 
            end
        else
            disp('Did not find result on disk...')
            lbsType = lbsTypes{lbsTypeID} ;

            % Choose which component of the velocity field to analyze
            if strcmpi(lbsType, 'v3d')
                v_lbs = v3d;
                titleString0 = '3D velocity $$\bf v$$';
                saveStrPrefix = 'LBS_Full_Velocity';
            elseif strcmpi(lbsType, 'vnScalar')
                v_lbs = vnScalar;
                titleString0 = 'normal velocity $$v_n$$';
                saveStrPrefix = 'LBS_Normal_Velocity_Scalar';
            elseif strcmpi(lbsType, 'vnVector')
                v_lbs = vn ;
                titleString0 = 'normal velocity, $$v_n \hat{n}$$';
                saveStrPrefix = 'LBS_Normal_Velocity_Vector';
            elseif strcmpi(lbsType, 'vt')
                v_lbs = vt;
                titleString0 = ['tangential velocity ' ...
                    '$${\bf v}_\parallel$$'];
                saveStrPrefix = 'LBS_Tangent_Velocity';
            elseif strcmpi(lbsType, 'H2vn')
                v_lbs = H2Vn;
                titleString0 = 'out-of-plane deformation $$2H v_n$$';
                saveStrPrefix = 'LBS_H2vn';
            elseif strcmpi(lbsType, 'divv')
                v_lbs = divVt;
                titleString0 = 'in-plane divergence $$\nabla \cdot {\bf v}_\parallel$$';
                saveStrPrefix = 'LBS_divv';
            elseif strcmpi(lbsType, 'gdot')
                v_lbs = gdotS;
                titleString0 = 'area rate of change $$\mathrm{Tr}[g^{-1} \dot{g}] / 2$$';
                saveStrPrefix = 'LBS_gdot';
            else
                error('Invalid velocity type for LBS analysis');
            end

            if convert_to_period
                v_lbs = v_lbs * T;
            end
            
            
            if numel(size(v_lbs)) == 3
                
                % Clip ends to account for poor temporal smoothing
                if nTimePoints2RmEnds > 0
                    v_lbs(:,:,1:nTimePoints2RmEnds) = [];
                    v_lbs(:,:,end-nTimePoints2RmEnds+1:end) = [];
                end
                
                if convert_to_period, v_lbs = v_lbs * T; end
                
                % Decompose the signals onto the Laplacian eigenbasis
                coeff = zeros(NmodesToView, 3, ntps-1);
                for tidx = 1:(ntps-1)
                    coeff(:,:,tidx) = LBVec' * v_lbs(:,:,tidx);
                end
                
            elseif numel(size(v_lbs)) == 2
                
                % Scalar fields should be converted from #t x #V to #V x #t
                if ~isequal(size(v_lbs), [tubi.nU*(tubi.nV-1), ntps-1])
                    v_lbs = v_lbs.';
                    if ~isequal(size(v_lbs), [tubi.nU*(tubi.nV-1), ntps-1])
                        error(['Improperly sized scalar field. ' ...
                            'Scalar fields must be defined on ' ...
                            'mesh vertices']);
                    end
                end
                
                % Clip ends to account for poor temporal smoothing
                if nTimePoints2RmEnds > 0
                    v_lbs(:, 1:nTimePoints2RmEnds) = [];
                    v_lbs(:, end-nTimePoints2RmEnds+1:end) = [];
                end
                
                if convert_to_period, v_lbs = v_lbs * T; end
                
                % Decompose the signals onto the Laplacian eigenbasis
                coeff = LBVec' * v_lbs;
                
            else
                
                error(['Invalid input. Each timepoints signal must be ' ...
                    '3D or a scalar field and defined on mesh vertices']);
                
            end
                
            % % Cat into a large array
            % if numel(size(v_lbs)) == 3
            %     % Note we convert v_pca into a #timepoints x 3*#faces array
            %     v_lbs = [ v_lbs(:,:,1), v_lbs(:,:,2), v_lbs(:,:,3) ];
            % else
            %     disp('v_lbs is 2D, meaning it is a scalar field.')
            % end

            toc
            fprintf(' Done\n');

            %% Generate Mode Visualizations -------------------------------
            % NOTE: This section doesn't really make sense in the context
            % of the LBS. I'm keeping it here for now because I don't want
            % to have to re-write it if I figure out I'm wrong later
            
            % % Extract Lagrangian mesh
            % V0c = lbsV; % allVc(:,:,t0idx);
            % DEC = lbsDEC; % DiscreteExteriorCalculus(Fc, V0c);
            % [V2Fc, F2Vc] = meshAveragingOperators(Fc, V0c);
            % VNc = per_vertex_normals(V0c, Fc, 'Weighting', 'angle');
            % 
            % % Chose the mode that will be visualized
            % displayTypes = {'normal', 'displacement', 'rotational', 'dilatational'} ;
            % for dispID = 1:length(displayTypes)
            % 
            %     plotType = displayTypes{dispID} ;
            %     for modeID = 1:NmodesToView
            % 
            %         if (modeID == 0)
            %             error('Broken handling');
            %             vmode = 0 .* coeff(:,1);
            %         else
            %             vmode = coeff(:, modeID);
            %         end
            % 
            %         % Extract visualization parameters
            %         [~, ~, ~, xyzlim] = tubi.getXYZLims();
            %         xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
            %         ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
            %         zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;
            % 
            %         % String for naming this mode figure
            %         saveStrMode = sprintf([saveStrPrefix, '_mode_%d'], modeID);
            % 
            %         clearvars divV rotV harmV scalarP vectorP VM
            % 
            %         % Treat the visualization differently if Vector/Scalar
            %         isVectorField = numel(vmode) == size(Fc,1)*3 ;
            %         if isVectorField
            %             vmode = reshape(vmode, size(Fc,1), 3);
            %             vmode_F = vmode;
            %             vmode = F2Vc * vmode;
            %             % Take components of this 3D vector field
            %             [divV, rotV, harmV, scalarP, vectorP] = ...
            %                 DEC.helmholtzHodgeDecomposition(vmode_F);
            % 
            %             FNc = faceNormal(triangulation(Fc, V0c));
            %             normV = dot(harmV, FNc, 2) .* FNc;
            %             harmV = harmV - normV;
            % 
            %             % Scale PCA compnents to produce visible displacements
            %             vscale = 400;
            %             VMc = V0c + vscale * vmode;
            % 
            %             saveStr = [saveStrMode, '_', plotType, '.png'];
            %             outputFigFn = fullfile(outputDirLBS, saveStr);
            %             if ~exist(outputFigFn, 'file') || overwriteImages
            % 
            %                 % fig = figure('Visible', 'off',  'units', 'centimeters', ...
            %                 %     'position', [0,0,xwidth,ywidth]) ;
            %                 fig = figure('Visible', 'off',  'units', 'centimeters') ;
            %                 for viewAngleID = 1:4
            %                     subplot(2, 2, viewAngleID)
            % 
            %                     quiver_lw = 0.75;
            %                     if strcmpi(plotType, 'displacement')
            % 
            %                         cmap = brewermap(256, '*RdBu');
            %                         crange = max(sqrt(sum(vmode.^2, 2))) * [-1 1];
            %                         vcolors = mapValueToColor( ...
            %                             sqrt(sum(vmode.^2, 2)) .* sign(dot(vmode, VNc, 2)), ...
            %                             crange, cmap);
            % 
            %                         patch( 'Faces', Fc, 'Vertices', VMc, 'FaceVertexCData', vcolors, ...
            %                             'FaceColor', 'interp', 'EdgeColor', 'none', ...
            %                             'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
            %                             'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
            %                             'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
            % 
            %                         titleString = [titleString0 ' - ' sprintf('mode %d', modeID)];
            % 
            %                         if convert_to_period
            %                             labelString = ...
            %                                 ['sgn$$[\vec{\bf V} \cdot \hat{n}_0] ||\vec{\bf V}||$$' ...
            %                                 ' [' tubi.spaceUnits '/T]'];
            %                         else
            %                             labelString = ...
            %                                 ['sgn$$[\vec{\bf V} \cdot \hat{n}_0] ||\vec{\bf V}||$$' ...
            %                                 ' [' tubi.spaceUnits '/' tubi.timeUnits ']'];
            %                         end
            % 
            %                     elseif strcmpi(plotType, 'rotational')
            % 
            %                         vectorPV = F2Vc * vectorP;
            %                         cmap = brewermap(256, '*RdBu');
            %                         crange = max(abs(vectorPV)) * [-1 1];
            % 
            %                         plotRotV = normalizerow(rotV);
            % 
            %                         patch( 'Faces', Fc, 'Vertices', V0c, ...
            %                             'FaceVertexCData', vectorPV, ...
            %                             'FaceColor', 'interp', 'EdgeColor', 'none', ...
            %                             'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
            %                             'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
            %                             'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
            % 
            %                         hold on
            % 
            %                         if plotArrowsOnModes
            %                             quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
            %                                 COM([ssf; ssf; sst]), ...
            %                                 plotRotV([sst; ssf; ssf]), plotRotV([ssf; sst; ssf]), ...
            %                                 plotRotV([ssf; ssf; sst]), ...
            %                                 1, 'LineWidth', quiver_lw, 'Color', 'k' );
            % 
            %                             hold off
            %                         end
            % 
            %                         titleString = [titleString0 ' - ' ...
            %                             sprintf('rotational part, mode %d', modeID)];
            % 
            %                     elseif strcmpi(plotType, 'dilatational')
            % 
            %                         cmap = brewermap(256, '*RdBu');
            %                         crange = max(abs(scalarP)) * [-1 1];
            % 
            %                         plotDivV = normalizerow(divV);
            % 
            %                         patch( 'Faces', Fc, 'Vertices', V0c, ...
            %                             'FaceVertexCData', scalarP, ...
            %                             'FaceColor', 'interp', 'EdgeColor', 'none', ...
            %                             'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
            %                             'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
            %                             'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
            % 
            %                         if plotArrowsOnModes
            %                             hold on
            %                             quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
            %                                 COM([ssf; ssf; sst]), ...
            %                                 plotDivV([sst; ssf; ssf]), plotDivV([ssf; sst; ssf]), ...
            %                                 plotDivV([ssf; ssf; sst]), ...
            %                                 1, 'LineWidth', quiver_lw, 'Color', 'k' );
            %                             hold off
            %                         end
            % 
            %                         titleString = [titleString0 ' - ' ...
            %                             sprintf('dilatational part, mode %d', modeID)];
            % 
            %                     elseif strcmpi(plotType, 'harmonic')
            % 
            %                         cmap = brewermap(512, '*RdBu');
            %                         cmap = cmap((1:256).', :);
            % 
            %                         absHVP = sqrt(sum((F2Vc * harmV).^2, 2));
            %                         crange = 0.01 * max(absHVP) * [0 1];
            % 
            %                         plotHarmV = normalizerow(harmV);
            % 
            %                         patch( 'Faces', Fc, 'Vertices', V0c, ...
            %                             'FaceVertexCData', absHVP, ...
            %                             'FaceColor', 'interp', 'EdgeColor', 'none', ...
            %                             'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
            %                             'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
            %                             'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
            % 
            %                         if plotArrowsOnModes
            %                             hold on
            %                             quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
            %                                 COM([ssf; ssf; sst]), ...
            %                                 plotHarmV([sst; ssf; ssf]), plotHarmV([ssf; sst; ssf]), ...
            %                                 plotHarmV([ssf; ssf; sst]), ...
            %                                 1, 'LineWidth', quiver_lw, 'Color', 'k' );
            % 
            %                             hold off
            %                         end
            % 
            %                         titleString = [titleString0 ' - ' ...
            %                             sprintf('harmonic part, mode %d', modeID)];
            % 
            %                         labelString = '$$||V_{harm}||$$';
            % 
            %                     elseif strcmpi(plotType, 'normal')
            % 
            %                         cmap = brewermap(256, '*RdBu');
            %                         crange = max(sqrt(sum((F2Vc*normV).^2, 2))) * [-1 1];
            %                         vcolors = mapValueToColor( ...
            %                             sqrt(sum((F2Vc*normV).^2, 2)) .* sign(dot(F2Vc * normV, VNc, 2)), ...
            %                             crange, cmap);
            % 
            %                         plotNormV = normalizerow(normV);
            % 
            %                         patch( 'Faces', Fc, 'Vertices', V0c, ...
            %                             'FaceVertexCData', vcolors, ...
            %                             'FaceColor', 'interp', 'EdgeColor', 'none', ...
            %                             'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
            %                             'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
            %                             'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
            % 
            %                         if plotArrowsOnModes
            %                             hold on;
            %                             quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
            %                                 COM([ssf; ssf; sst]), ...
            %                                 plotNormV([sst; ssf; ssf]), plotNormV([ssf; sst; ssf]), ...
            %                                 plotNormV([ssf; ssf; sst]), ...
            %                                 1, 'LineWidth', quiver_lw, 'Color', 'k' );
            % 
            %                             hold off
            %                         end
            % 
            %                         titleString = [titleString0 ' - ' ...
            %                             sprintf('normal part, mode %d', modeID)];
            %                     else
            %                         error('Invalid plot type');
            %                     end
            % 
            %                     axis equal
            % 
            %                     scaleFactor = 1.1;
            %                     xlim(scaleFactor * [xmin, xmax]);
            %                     ylim(scaleFactor * [ymin, ymax]);
            %                     zlim(scaleFactor * [zmin, zmax]);
            % 
            %                     % camlight
            %                     xlabel(['AP Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            %                     ylabel(['Lateral Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            %                     zlabel(['DV Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            % 
            %                     colormap(cmap);
            %                     set(gca, 'Clim', crange);
            % 
            %                     % cb = colorbar;
            %                     % cb.Label.String = labelString;
            %                     % cb.Label.Interpreter = 'latex';
            % 
            %                     set(gcf, 'Color', [1 1 1]);
            %                     axis off
            % 
            %                     if viewAngleID == 1
            %                         if contains(titleString, '$')
            %                             sgtitle(titleString, 'Interpreter', 'latex');
            %                         else
            %                             sgtitle(titleString, 'FontWeight', 'normal');
            %                         end
            %                     end
            % 
            %                     set(gca, 'FontSize', 6);
            %                     set(gca, 'FontWeight', 'normal');
            % 
            %                     view([0,(viewAngleID-1)*90])
            %                 end
            % 
            %                 % Resize Figure for Paper -------------------------------------------------
            %                 set(fig, 'Units', 'centimeters');
            % 
            %                 % ratio = fig.Position(4) ./ fig.Position(3);
            %                 % fig.Position(3) = 4;
            %                 % fig.Position(4) = ratio * fig.Position(3);
            % 
            %                 set(fig, 'PaperPositionMode', 'auto');
            %                 set(fig, 'position', [0, 0, 16, 10])
            %                 set(fig.Children, 'FontName', 'Helvetica');
            % 
            %                 % Save figure
            %                 disp(['Saving figure: ' outputFigFn])
            %                 export_fig(outputFigFn, '-png', '-r500');
            %                 close all
            %             end
            %         else
            %             disp('PCA field is scalar field')
            % 
            %             % Plot displacement of vertices from reference mesh
            %             if strcmpi(lbsType, 'vnScalar') && strcmpi(plotType, 'displacement')
            %                 % Scale PCA compnents to produce visible displacements
            %                 vscale = 400;
            % 
            %                 VMc = V0c + vscale * vmode .* VNc;
            %                 saveStr = [saveStrMode, '_', plotType, '.png'];
            %                 outputFigFn = fullfile(outputDirLBS, saveStr);
            %                 if ~exist(outputFigFn, 'file') || overwriteImages
            % 
            %                     fig = figure('Visible', 'off',  'units', 'centimeters') ;
            %                     for viewAngleID = 1:4
            %                         subplot(2, 2, viewAngleID)
            % 
            %                         cmap = brewermap(256, '*RdBu');
            %                         crange = max(sqrt(sum(vmode.^2, 2))) * [-1 1];
            %                         vcolors = mapValueToColor( vmode, ...
            %                             crange, cmap);
            % 
            %                         patch( 'Faces', Fc, 'Vertices', VMc, 'FaceVertexCData', vcolors, ...
            %                             'FaceColor', 'interp', 'EdgeColor', 'none', ...
            %                             'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
            %                             'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
            %                             'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
            % 
            %                         titleString = [titleString0 ' - ' sprintf('mode %d', modeID)];
            % 
            %                         if convert_to_period
            %                             labelString = ...
            %                                 ['v_n' ...
            %                                 ' [' tubi.spaceUnits '/T]'];
            %                         else
            %                             labelString = ...
            %                                 ['V_n' ...
            %                                 ' [' tubi.spaceUnits '/' tubi.timeUnits ']'];
            %                         end
            %                         axis equal
            % 
            %                         scaleFactor = 1.1;
            %                         xlim(scaleFactor * [xmin, xmax]);
            %                         ylim(scaleFactor * [ymin, ymax]);
            %                         zlim(scaleFactor * [zmin, zmax]);
            % 
            %                         xlabel(['AP Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            %                         ylabel(['Lateral Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            %                         zlabel(['DV Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            % 
            %                         colormap(cmap);
            %                         set(gca, 'Clim', crange);
            % 
            %                         % cb = colorbar;
            %                         % cb.Label.String = labelString;
            %                         % cb.Label.Interpreter = 'latex';
            % 
            %                         set(gcf, 'Color', [1 1 1]);
            %                         axis off
            % 
            %                         if viewAngleID == 1
            %                             if contains(titleString, '$')
            %                                 sgtitle(titleString, 'Interpreter', 'latex');
            %                             else
            %                                 sgtitle(titleString, 'FontWeight', 'normal');
            %                             end
            %                         end
            % 
            %                         set(gca, 'FontSize', 6);
            %                         set(gca, 'FontWeight', 'normal');
            % 
            %                         view([0,(viewAngleID-1)*90])
            %                     end
            % 
            %                     % Resize Figure for Paper -----------------
            %                     set(fig, 'Units', 'centimeters');
            %                     set(fig, 'PaperPositionMode', 'auto');
            %                     set(fig, 'position', [0, 0, 16, 10])
            %                     set(fig.Children, 'FontName', 'Helvetica');
            % 
            %                     % Save figure
            %                     disp(['Saving figure: ' outputFigFn])
            %                     export_fig(outputFigFn, '-png', '-r500');
            %                     close all
            %                 end
            % 
            %             elseif strcmpi(plotType, 'normal')
            % 
            %                 % Save image of scalar field alone
            %                 saveStr = [saveStrMode '.png'];
            %                 outputFigFn = fullfile(outputDirLBS, saveStr);
            %                 if ~exist(outputFigFn, 'file') || overwriteImages
            % 
            %                     % fig = figure('Visible', 'off',  'units', 'centimeters', ...
            %                     %     'position', [0,0,xwidth,ywidth]) ;
            %                     fig = figure('Visible', 'off',  'units', 'centimeters') ;
            %                     for viewAngleID = 1:4
            %                         subplot(2, 2, viewAngleID)
            % 
            %                         cmap = brewermap(256, '*RdBu');
            %                         crange = max(sqrt(sum(vmode.^2, 2))) * [-1 1];
            %                         vcolors = mapValueToColor( vmode, ...
            %                             crange, cmap);
            % 
            %                         patch( 'Faces', Fc, 'Vertices', V0c, 'FaceVertexCData', vcolors, ...
            %                             'FaceColor', 'interp', 'EdgeColor', 'none', ...
            %                             'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
            %                             'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
            %                             'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );
            % 
            %                         titleString = [titleString0 ' - ' sprintf('mode %d', modeID)];
            % 
            %                         scaleFactor = 1.1;
            %                         xlim(scaleFactor * [xmin, xmax]);
            %                         ylim(scaleFactor * [ymin, ymax]);
            %                         zlim(scaleFactor * [zmin, zmax]);
            % 
            %                         xlabel(['AP Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            %                         ylabel(['Lateral Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            %                         zlabel(['DV Position [' tubi.spaceUnits ']'], 'interpreter', 'latex');
            % 
            %                         colormap(cmap);
            %                         set(gca, 'Clim', crange);
            %                         set(gcf, 'Color', [1 1 1]);
            %                         axis off
            % 
            %                         if viewAngleID == 1
            %                             if contains(titleString, '$')
            %                                 sgtitle(titleString, 'Interpreter', 'Latex');
            %                             else
            %                                 sgtitle(titleString, 'FontWeight', 'normal');
            %                             end
            %                         end
            %                         set(gca, 'FontSize', 6);
            %                         set(gca, 'FontWeight', 'normal');
            % 
            %                         view([0,(viewAngleID-1)*90])
            %                     end
            % 
            %                     % Resize Figure for Paper -----------------
            %                     set(fig, 'Units', 'centimeters');
            %                     set(fig, 'PaperPositionMode', 'auto');
            %                     set(fig, 'position', [0, 0, 16, 10])
            %                     set(fig.Children, 'FontName', 'Helvetica');
            % 
            %                     % Save figure
            %                     disp(['Saving figure: ' outputFigFn])
            %                     export_fig(outputFigFn, '-png', '-r500');
            %                     close all
            %                 end
            %             end
            %         end
            %     end
            % end
            % 
            % % clear FN1 vn vt pcaType v_pca 
            % % clear V0 V2F F2V VN DEC modeID titleString
            % % clear vmode VM lapVMode xyzlim xmin ymin zmin xmax ymax zmax
            % % clear cmap vrange vcolors scaleFactor

            %% Generate 2D LBS Mode Comparison Plot ===================================
            close all; clc;

            % Generate Visualization --------------------------------------------------
            % color by timestamp in proper units (tubi.timeInterval)
            tps = tubi.xp.fileMeta.timePoints(1:end-1) - tubi.t0;
            % Remove timepoints that we removed from the LBS
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
                outputFigFn = fullfile(outputDirLBS, saveStrPair);
                if ~exist(outputFigFn, 'file') || overwriteImages
                    
                    % The projection of the velocities along the component axes
                    % v1 = dot(v_lbs.', repmat(coeff(:, axID1), 1, size(v_lbs,1)), 1).';
                    % v2 = dot(v_lbs.', repmat(coeff(:, axID2), 1, size(v_lbs,1)), 1).';
                    
                    if size(v_lbs, 3) == 1
                        v1 = coeff(axID1, :).';
                        v2 = coeff(axID2, :).';
                    else
                        disp(['2D mode comparison plots not yet ' ...
                            'implemented for nonscalar signals']);
                    end

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

                    if drawArrowsInLBSPlane
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
                        xlabel(['Laplacian component ' num2str(axID1) ...
                            ' (' tubi.spaceUnits '/' tubi.timeUnits ')'], 'interpreter', 'latex');
                        ylabel(['Laplacian component ' num2str(axID2) ...
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

            %% Generate 3D LBS Mode Comparison Plot ===================================
            % Compares the 3 most important modes produced by PCA
            close all; clc;


            % Generate Visualization --------------------------------------------------
            outputFigFn = fullfile(outputDirLBS, [saveStrPrefix '_3D_trajectory.png']);
            if ~exist(outputFigFn, 'file') || overwriteImages 
                
                % The projection of the velocities along the component axes
                % v1 = dot(v_lbs.', repmat(coeff(:, 1), 1, size(v_lbs,1)), 1).';
                % v2 = dot(v_lbs.', repmat(coeff(:, 2), 1, size(v_lbs,1)), 1).';
                % v3 = dot(v_lbs.', repmat(coeff(:, 3), 1, size(v_lbs,1)), 1).';
                
                if size(v_lbs, 3) == 1
                    v1 = coeff(1, :).';
                    v2 = coeff(2, :).';
                    v3 = coeff(3, :).';
                else
                    disp(['3D mode comparison plots not yet ' ...
                        'implemented for nonscalar signals']);
                end

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

                if  drawArrowsInLBS3D 
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

                if drawArrowsInLBS3D 
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

            % Determine the Error Associated with Removing Each Mode Individually -----

            lbsModeError = zeros(size(LBVec,2), 1);
            for i = 1:numel(lbsModeError)
                if mod(i, 20) == 0
                    disp(['Checking for contribution from mode ' num2str(i)])
                end

                % The projection of the velocities along the current component axis
                % vi = dot(v_lbs.', repmat(coeff(:, i), 1, size(v_lbs,1)), 1) .* ...
                %     repmat(coeff(:, i), 1, size(v_lbs,1));
                vi = coeff(i) .* LBVec(:,i);

                curErr = sum(sum((vi.').^2, 2) ./ sum(v_lbs.^2, 2), 1) ./ size(v_lbs,1);

                lbsModeError(i) = curErr;

            end

            % Generate Visualization --------------------------------------------------

            plot_lw = 0.75;
            scatter_size = 2;

            fig = figure('Visible', 'off',  'units', 'centimeters') ;

            stem(1:numel(lbsModeError), lbsModeError.', ...
                'LineStyle', '-', ...
                'MarkerFaceColor', tubi.plotting.colors(1,:), ...
                'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
                'LineWidth', plot_lw, ...
                'MarkerSize', scatter_size);

            xlim([1 numel(lbsModeError)]);

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

            saveStr = ['Mode_Contribution_' lbsType '.pdf'];
            saveStr = fullfile(outputDirLBS, saveStr);
            export_fig(saveStr, '-pdf', '-r200');

            %% Save everything to disk
            disp(['Saving results to disk: ' fnResult])
            save(fnResult, 'lbsModeError', 'coeff', ...
                'lbsType', 'convert_to_period', 'T', 't0', ...
                'NmodesToView', 'nTimePoints2RmEnds', 'meshChoice') ;
            if nargout > 0
                results{dmykk} = struct('lbsModeError', lbsModeError, ...
                    'coeff', coeff, ...
                    'lbsType', lbsType, ...
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
