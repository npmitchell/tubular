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
% pairs of mode numbers to plot in 2d projections (modes are ordered by rank)
axPairs = [1,2;1,3;2,3] ;


%% Overwrite default options with supplied options
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = overwrite ;
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

if isfield(options, 'meshChoice')
    meshChoice = options.meshChoice ;
    isa(meshChoice, 'char')
    meshChoice = {meshChoice} ;
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
        allVx = reshape(permute(pathlines.vertices3d.vXrs, [2,3,1]), [tubi.nU*tubi.nV, 1, ntps]);
        allVy = reshape(permute(pathlines.vertices3d.vYrs, [2,3,1]), [tubi.nU*tubi.nV, 1, ntps]);
        allVz = reshape(permute(pathlines.vertices3d.vZrs, [2,3,1]), [tubi.nU*tubi.nV, 1, ntps]);
        allV = cat(2, allVx, allVy, allVz) ;

        % face list is static in Lagrangian frame  
        F = pathlines.refMesh.f ;

        % allF = repmat(F, [1, 1, ntps]) ;
        allFN = zeros(size(F,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
        for tidx = 1:numel(tubi.xp.fileMeta.timePoints)
            V = squeeze(allV(:,:,tidx)) ;

            % Generate coordinate frame field on each face
            TR = triangulation(F,V);
            fn = TR.faceNormal;
            allFN(:,:,tidx) = fn;
        end
        
    elseif strcmpi(meshStyle, 'sphi')

        %% Load Lagrangian Meshes =================================================
        close all; clc;

        allV = [];
        allFN = [];
        for tidx = 1:numel(tubi.xp.fileMeta.timePoints)
            progressbar(tidx, numel(tubi.xp.fileMeta.timePoints));

            % Load data
            t = tubi.xp.fileMeta.timePoints(tidx);
            tubi.setTime(t); % Set fit time to the current time
            spcutMeshSmRSC = tubi.getCurrentSPCutMeshSmRSC; % The Lagrangian mesh

            % Extract mesh
            F = spcutMeshSmRSC.f; % FACE ORDER HAS VERTEX NORMALS POINTING INWARD
            V = spcutMeshSmRSC.v; % ASK NOAH ABOUT THIS
            u = spcutMeshSmRSC.u;

            % Generate coordinate frame field on each face
            TR = triangulation(F,V);
            if isempty(allV)

                allV = zeros(size(V,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
                allFN = zeros(size(F,1), 3, numel(tubi.xp.fileMeta.timePoints)-1);
            end

            allV(:,:,tidx) = V;
            allFN(:,:,tidx) = TR.faceNormal;

        end
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

    clear spcutMeshSmRSC COM ssf sst TR
    clear tidx V u TR 

    %% Visualize Velocity PCA Modes ===========================================
    close all; clc;

    % The full face-based 3D velocity field (#t x #F x 3)
    fprintf('First computing face-based velocity field over time... ');
    v3d = zeros(ntps-1, size(F,1), 3) ;    
    for tidx = 1:numel(tubi.xp.fileMeta.timePoints) -1 
        V0 = squeeze(allV(:,:,tidx)) ; % vertices at t0
        V1 = squeeze(allV(:,:,tidx+1)) ; % vertices at t1

        % velocities on vertices (Lagrangian frame) in units 
        % of tubi.spaceUnits / [tubi.timeUnits * (timepoint1-timepoint0)]
        vel = V1 - V0 ;

        % Push the velocities to the face barycenters
        [V2F, F2V] = meshAveragingOperators(F, V0) ;
        velF = V2F * vel ;

        % Collimate velocities on barycenters
        v3d(tidx, :, :) = velF ;
    end

    % THIS IS ONLY APPROXIMATE SINCE (s,phi) is not fully Lagrangian
    % v3d = tubi.getVelocityAverage().vf;
    fprintf('Done\n');

    %% Load velocities and Perform PCA -----------------------------------------
    fprintf('Performing PCA... ');
    tic
    % The normal component of the velocity field
    FN1 = permute(allFN, [3 1 2]); 
    FN1(end, :, :) = [];
    vn = sum( FN1 .* v3d, 3 ) .* FN1;

    % The tangential component of each velocity field
    vt = v3d - vn;

    % Consider doing PCA on each kind of velocity field
    pcaTypes = {'v3d', 'vn', 'vt'};
    for pcaTypeID = 1:length(pcaTypes)
        pcaType = pcaTypes{pcaTypeID} ;
        fnResult = fullfile(outputDirRoot, sprintf('pcaResults_%s.mat', pcaType)) ;
        outputDirPCA = fullfile(outputDirRoot, pcaType) ;
        
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
                titleString = 'Full 3D Velocity $$\vec{\bf V}$$';
                saveStrPrefix = 'PCA_Full_Velocity';
            elseif strcmpi(pcaType, 'vn')
                v_pca = vn;
                titleString = 'Normal Velocity $$\vec{\bf V} \cdot \hat{n}$$';
                saveStrPrefix = 'PCA_Normal_Velocity';
            elseif strcmpi(pcaType, 'vt')
                v_pca = vt;
                titleString = ['Tangential Velocity ' ...
                    '$$\vec{\bf V}-(\vec{\bf V} \cdot \hat{n}) \hat{n}$$'];
                saveStrPrefix = 'PCA_Tangent_Velocity';
            else
                error('Invalid velocity type for PCA analysis');
            end

            if convert_to_period
                v_pca = v_pca * T;
            end

            % Clip ends to account for poor temporal smoothing
            v_pca = [ v_pca(:,:,1), v_pca(:,:,2), v_pca(:,:,3) ];
            if nTimePoints2RmEnds > 0
                v_pca(1:nTimePoints2RmEnds, :) = [];
                v_pca(end-nTimePoints2RmEnds+1:end, :) = [];
                v_pca = v_pca * T; 
            end

            % Perform PCA on the Cartesian components of the velocity
            coeff = pca(v_pca);

            toc
            fprintf(' Done\n');

            %% Generate Visualization --------------------------------------------------

            % Extract Lagrangian mesh
            V0 = allV(:,:,t0idx);
            DEC = DiscreteExteriorCalculus(F, V0);
            [V2F, F2V] = meshAveragingOperators(F, V0);
            VN = per_vertex_normals(V0, F, 'Weighting', 'angle');

            % Chose the mode that will be visualized
            displayTypes = {'normal', 'displacement', 'rotational', 'dilatational'} ;
            for dispID = 1:length(displayTypes)

                plotType = displayTypes{dispID} ;
                for modeID = 1:NmodesToView

                    if (modeID == 0)
                        error('Broken handling');
                        vmode = 0 .* coeff(:,1);
                        titleString = sprintf('Lagrangian Mesh T = %d', tubi.t0);
                    else
                        vmode = coeff(:, modeID);
                    end

                    vmode = reshape(vmode, size(F,1), 3);
                    vmode_F = vmode;
                    vmode = F2V * vmode;

                    [divV, rotV, harmV, scalarP, vectorP] = ...
                        DEC.helmholtzHodgeDecomposition(vmode_F);

                    FN = faceNormal(triangulation(F, V0));
                    normV = dot(harmV, FN, 2) .* FN;
                    harmV = harmV - normV;

                    % Scale PCA compnents to produce visible displacements
                    vscale = 1.5e3;
                    VM = V0 + vscale * vmode;

                    % Extract visualization parameters
                    [~, ~, ~, xyzlim] = tubi.getXYZLims();
                    xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
                    ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
                    zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;

                    % figure parameters
                    xwidth = 16 ; % cm
                    ywidth = 10 ; % cm

                    % String for naming this mode figure
                    saveStrMode = sprintf([saveStrPrefix, '_Mode_%d'], modeID);
                    
                    saveStr = [saveStrMode, '_', plotType, '.png'];
                    outputFigFn = fullfile(outputDirPCA, saveStr);
                    if ~exist(outputFigFn, 'file') || overwriteImages
                        
                        % fig = figure('Visible', 'on',  'units', 'centimeters', ...
                        %     'position', [0,0,xwidth,ywidth]) ;
                        fig = figure('Visible', 'on',  'units', 'centimeters') ;
                        for viewAngleID = 1:4
                            subplot(2, 2, viewAngleID)

                            quiver_lw = 0.75;
                            if strcmpi(plotType, 'displacement')

                                cmap = brewermap(256, '*RdBu');
                                crange = max(sqrt(sum(vmode.^2, 2))) * [-1 1];
                                vcolors = mapValueToColor( ...
                                    sqrt(sum(vmode.^2, 2)) .* sign(dot(vmode, VN, 2)), ...
                                    crange, cmap);

                                patch( 'Faces', F, 'Vertices', VM, 'FaceVertexCData', vcolors, ...
                                    'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                    'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                    'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                    'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                hold on

                                % ssf = false( size(V0,1), 1 );
                                % sst = false( size(V0,1), 1 );
                                % sst(1:2:end) = true;
                                % % sst = true( size(V0, 1), 1 );
                                %
                                % quiver3( V0([sst; ssf; ssf]), V0([ssf; sst; ssf]), ...
                                %     V0([ssf; ssf; sst]), ...
                                %     vmode([sst; ssf; ssf]), vmode([ssf; sst; ssf]), ...
                                %     vmode([ssf; ssf; sst]), ...
                                %     1, 'LineWidth', 2, 'Color', 'k' );

                                hold off

                                % titleString = [titleString, sprintf('\nMode %d', modeID)];
                                titleString = sprintf('Mode %d', modeID);

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

                                vectorPV = F2V * vectorP;
                                cmap = brewermap(256, '*RdBu');
                                crange = max(abs(vectorPV)) * [-1 1];

                                plotRotV = normalizerow(rotV);

                                patch( 'Faces', F, 'Vertices', V0, ...
                                    'FaceVertexCData', vectorPV, ...
                                    'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                    'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                    'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                    'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                hold on

                                COM = barycenter(V0, F);
                                ssf = false( size(F,1), 1 );
                                sst = false( size(F,1), 1 );
                                sst(1:7:end) = true;
                                % sst = true( size(F, 1), 1 );

                                quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                    COM([ssf; ssf; sst]), ...
                                    plotRotV([sst; ssf; ssf]), plotRotV([ssf; sst; ssf]), ...
                                    plotRotV([ssf; ssf; sst]), ...
                                    1, 'LineWidth', quiver_lw, 'Color', 'k' );

                                hold off

                                % titleString = [titleString, ...
                                %     sprintf('\n Rotational Part, Mode %d', modeID)];
                                titleString = sprintf('Rotational Part, Mode %d', modeID);

                                labelString = 'Vector Potential';

                            elseif strcmpi(plotType, 'dilatational')

                                cmap = brewermap(256, '*RdBu');
                                crange = max(abs(scalarP)) * [-1 1];

                                plotDivV = normalizerow(divV);

                                patch( 'Faces', F, 'Vertices', V0, ...
                                    'FaceVertexCData', scalarP, ...
                                    'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                    'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                    'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                    'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                hold on

                                COM = barycenter(V0, F);
                                ssf = false( size(F,1), 1 );
                                sst = false( size(F,1), 1 );
                                sst(1:7:end) = true;
                                % sst = true( size(F, 1), 1 );

                                quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                    COM([ssf; ssf; sst]), ...
                                    plotDivV([sst; ssf; ssf]), plotDivV([ssf; sst; ssf]), ...
                                    plotDivV([ssf; ssf; sst]), ...
                                    1, 'LineWidth', quiver_lw, 'Color', 'k' );

                                hold off

                                % titleString = [titleString, ...
                                %     sprintf('\n Dilatational Part, Mode %d', modeID)];
                                titleString = sprintf('Dilatational Part, Mode %d', modeID);

                                labelString = 'Scalar Potential';

                            elseif strcmpi(plotType, 'harmonic')

                                cmap = brewermap(512, '*RdBu');
                                cmap = cmap((1:256).', :);

                                absHVP = sqrt(sum((F2V * harmV).^2, 2));
                                crange = 0.01 * max(absHVP) * [0 1];

                                plotHarmV = normalizerow(harmV);

                                patch( 'Faces', F, 'Vertices', V0, ...
                                    'FaceVertexCData', absHVP, ...
                                    'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                    'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                    'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                    'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                hold on

                                COM = barycenter(V0, F);
                                ssf = false( size(F,1), 1 );
                                sst = false( size(F,1), 1 );
                                sst(1:7:end) = true;
                                % sst = true( size(F, 1), 1 );

                                quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                    COM([ssf; ssf; sst]), ...
                                    plotHarmV([sst; ssf; ssf]), plotHarmV([ssf; sst; ssf]), ...
                                    plotHarmV([ssf; ssf; sst]), ...
                                    1, 'LineWidth', quiver_lw, 'Color', 'k' );

                                hold off

                                % titleString = [titleString, ...
                                %     sprintf('\n Harmonic Part, Mode %d', modeID)];
                                titleString = sprintf('\n Harmonic Part, Mode %d', modeID);

                                labelString = '$$||V_{harm}||$$';

                            elseif strcmpi(plotType, 'normal')

                                cmap = brewermap(256, '*RdBu');
                                crange = max(sqrt(sum((F2V*normV).^2, 2))) * [-1 1];
                                vcolors = mapValueToColor( ...
                                    sqrt(sum((F2V*normV).^2, 2)) .* sign(dot(F2V * normV, VN, 2)), ...
                                    crange, cmap);

                                plotNormV = normalizerow(normV);

                                patch( 'Faces', F, 'Vertices', V0, ...
                                    'FaceVertexCData', vcolors, ...
                                    'FaceColor', 'interp', 'EdgeColor', 'none', ...
                                    'FaceLighting', 'gouraud', 'FaceAlpha', 1, ...
                                    'SpecularStrength', 0.2, 'DiffuseStrength', 0.5, ...
                                    'AmbientStrength', 0.3, 'SpecularExponent', 0.5 );

                                hold on

                                COM = barycenter(V0, F);
                                ssf = false( size(F,1), 1 );
                                sst = false( size(F,1), 1 );
                                sst(1:7:end) = true;
                                % sst = true( size(F, 1), 1 );

                                quiver3( COM([sst; ssf; ssf]), COM([ssf; sst; ssf]), ...
                                    COM([ssf; ssf; sst]), ...
                                    plotNormV([sst; ssf; ssf]), plotNormV([ssf; sst; ssf]), ...
                                    plotNormV([ssf; ssf; sst]), ...
                                    1, 'LineWidth', quiver_lw, 'Color', 'k' );

                                hold off

                                % titleString = [titleString, ...
                                %     sprintf('\n Harmonic Part, Mode %d', modeID)];
                                titleString = sprintf('\n Normal Part, Mode %d', modeID);

                                labelString = '$$||V_{N}||$$';

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

                            title(titleString, 'FontWeight', 'normal');
                            % title(titleString, 'interpreter', 'latex');
                            % set(gca, 'FontSize', 7);

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
                        set(fig.Children, 'FontName', 'Helvetica');

                        % Save figure
                        disp(['Saving figure: ' outputFigFn])
                        export_fig(outputFigFn, '-png', '-r500');
                        close all
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

                saveStrPair = sprintf([saveStr '_m%dm%d.pdf'], axID1, axID2);
                
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
                    scaleFactor = 1.25;
                    % xLim = 300 * [-1 1]; yLim = xLim;
                    % if convert_to_period
                    %     xLim = T * xLim;
                    %     yLim = T * yLim;
                    % end

                    fig = figure('Visible', 'on',  'units', 'centimeters') ;
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

                        % xlabel(['Velocity mode ' num2str(axID1) ...
                        %     ' (' tubi.spaceUnits '/T)'], 'interpreter', 'latex');
                        % ylabel(['Velocity mode ' num2str(axID2) ...
                        %     ' (' tubi.spaceUnits '/T)'], 'interpreter', 'latex');

                        xlabel(['Mode ' num2str(axID1) ' [' tubi.spaceUnits '/T]']);
                        ylabel(['Mode ' num2str(axID2) ' [' tubi.spaceUnits '/T]']);

                    else

                        xlabel(['Velocity principal component ' num2str(axID1) ...
                            ' (' tubi.spaceUnits '/' tubi.timeUnits ')'], 'interpreter', 'latex');
                        ylabel(['Velocity principal component ' num2str(axID2) ...
                            ' (' tubi.spaceUnits '/' tubi.timeUnits ')'], 'interpreter', 'latex');
                    end

                    % xticks([-2000 0 2000]);
                    % yticks([-2000 0 2000]);


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

                % title(titleString, 'interpreter', 'latex');

                set(gcf, 'Color', [1 1 1]);

                % Resize Figure for Paper -------------------------------------------------
                set(fig, 'Units', 'centimeters');

                % ratio = fig.Position(4) ./ fig.Position(3);
                % fig.Position(3) = 4;
                % fig.Position(4) = ratio * fig.Position(3);

                % set(gca, 'FontSize', 10);
                set(gca, 'FontWeight', 'normal');

                set(fig, 'PaperPositionMode', 'auto');
                set(fig.Children, 'FontName', 'Helvetica');

                export_fig(outputFigFn, '-pdf', '-r200');
            end
            end
            % clear v3d FN1 vn vt pcaType v_pca coeff
            % clear axID1 axID2 v1 v2 vvec scaleFactor xLim yLim zLim
            % clear cmap timeColors cb titleString vScale

            %% Generate 3D PCA Mode Comparison Plot ===================================
            % Compares the 3 most important modes produced by PCA
            close all; clc;


            % Generate Visualization --------------------------------------------------

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
            r = 5;
            if convert_to_period, r = T * r; end

            quiver_lw = 0.75;
            fig = figure('Visible', 'on',  'units', 'centimeters') ;
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

                xlabel(['Mode ' num2str(1)]);
                ylabel(['Mode ' num2str(2)]);
                zlabel(['Mode ' num2str(3)]);

            else
                xlabel(['Mode ' num2str(1) ...
                    ' [' tubi.spaceUnits '/' tubi.timeUnits ']'], 'interpreter', 'latex');
                ylabel(['Mode ' num2str(2) ...
                    ' [' tubi.spaceUnits '/' tubi.timeUnits ']'], 'interpreter', 'latex');
                zlabel(['Mode ' num2str(3) ...
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
            cb.Ticks = [0.2 1 2];
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

            % title(titleString);
            % title(titleString, 'interpreter', 'latex');

            set(gcf, 'Color', [1 1 1]);

            % Resize Figure for Paper -------------------------------------------------
            set(fig, 'Units', 'centimeters');

            ratio = fig.Position(4) ./ fig.Position(3);
            fig.Position(3) = 5.5;
            fig.Position(4) = ratio * fig.Position(3);

            set(gca, 'FontSize', 5);
            set(gca, 'FontWeight', 'normal');

            set(fig, 'PaperPositionMode', 'auto');
            set(fig.Children, 'FontName', 'Helvetica');


            outputFigFn = fullfile(outputDirPCA, saveStr);
            % export_fig(saveStr, '-pdf', '-r200');
            export_fig(outputFigFn, '-png', '-r500');

            % clear v3d FN1 vn vt pcaType v_pca coeff
            % clear axID1 axID2 v1 v2 v3 vvec scaleFactor xLim yLim zLim
            % clear cmap timeColors cb titleString
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

            fig = figure('Visible', 'on',  'units', 'centimeters') ;

            stem(1:numel(pcaModeError), pcaModeError.', ...
                'LineStyle', '-', ...
                'MarkerFaceColor', tubi.plotting.colors(1,:), ...
                'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
                'LineWidth', plot_lw, ...
                'MarkerSize', scatter_size);

            xlim([1 numel(pcaModeError)]);
            % ylim([0 0.5]);
            % yticks(0:0.1:0.5);

            xlabel('Mode number');
            ylabel('Mode contribution');

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

            % clear v3d FN1 vn vt pcaType v_pca
            % clear coeff score latent tsquared explained mu
            % clear pcaModeError i vi currErr plot_lw scatter_size ratio

            %% Save everything to disk
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
