function plotSPCutMeshSmSeriesUtility(QS, coordsys, options)
%plotSPCutMeshSmSeriesUtility(QS, coordsys, options)
%   Using (s,phi) pullback cutmeshes, smooth coord system in time via
%   simple triangular pulse averaging of positions in embedding space
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   Current QuapSlap object
% coordsys : str specifier ('spcutMeshSm', 'spcutMeshSmRS', 'spcutMeshSmRSC')
%   What kind of mesh to plot. 
%   Any string without 'RS' or 'RSC' directs to plotting 'spcutMeshSmRSC'
%   Any string with RS directs to plotting 'spcutMeshSm'
%   Any string with RSC directs to plotting 'spcutMeshSmRS'
% options : struct with fields
%   overwrite : bool, default=false
%       overwrite results on disk
%   preivew : bool, default=QS.plotting.preview
%       display intermediate results
%
% NPMitchell 2020

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
[~, ~, ~, xyzlim] = QS.getXYZLims() ;
timePoints = QS.xp.fileMeta.timePoints ;

%% Unpack options
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end

%% Create plots
if contains(coordsys, 'RS')
    imDir = fullfile(QS.dir.spcutMeshSmRS, 'images') ;
    disp(['Output to SmRS dir: ' imDir])
elseif contains(coordsys, 'RSC')
    imDir = fullfile(QS.dir.spcutMeshSmRSC, 'images') ;
    disp(['Output to SmRSC dir: ' imDir])
else
    imDir = fullfile(QS.dir.spcutMeshSm, 'images') ;
    disp(['Output to spcutMeshSm dir: ' imDir])
end

if ~exist(imDir, 'dir')
    mkdir(imDir)
end
pdir = ensureDir(fullfile(imDir, 'perspective')) ;
ddir = ensureDir(fullfile(imDir, 'dorsal')) ;
vdir = ensureDir(fullfile(imDir, 'ventral')) ;
ldir = ensureDir(fullfile(imDir, 'latL')) ;
rdir = ensureDir(fullfile(imDir, 'latR')) ;
for qq = 1:length(timePoints)
    tt = timePoints(qq) ;
    disp(['checking figures for time-smoothed meshes: t = ' num2str(tt)])
    % prep directories
    fig1fn = fullfile(imDir, [sprintf('%04d', tt ) '.png']) ;
    fp0fn = fullfile(pdir, [sprintf('%04d', tt ) '.png']) ;
    fp1fn = fullfile(ddir, ['dorsal_' sprintf('%04d', tt ) '.png']) ;
    fp2fn = fullfile(vdir, ['ventral_' sprintf('%04d', tt ) '.png']) ;
    fp3fn = fullfile(ldir, ['latL_' sprintf('%04d', tt ) '.png']) ;
    fp4fn = fullfile(rdir, ['latR_' sprintf('%04d', tt ) '.png']) ;
    e0 = ~exist(fig1fn, 'file') ;
    e1 = ~exist(fp0fn, 'file') ;
    e2 = ~exist(fp1fn, 'file') || ~exist(fp2fn, 'file') ;
    e3 = ~exist(fp3fn, 'file') || ~exist(fp4fn, 'file') ;
    if e0 || e1 || e2 || e3 || overwrite
        disp(['Consider smoothed sphi gridmesh for time ' num2str(tt)])
        % Load the spcutMesh for this timepoint
        
        if contains(coordsys, 'RS')
            load(sprintf(QS.fullFileBase.spcutMeshSmRS, tt), 'spcutMeshSmRS') ;
            vqq = spcutMeshSmRS.v ;
            nqq = spcutMeshSmRS.vn ;
            faces = spcutMeshSmRS.f ;
        elseif contains(coordsys, 'RSC')
            load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tt), 'spcutMeshSmRSC') ;
            vqq = spcutMeshSmRSC.v ;
            nqq = spcutMeshSmRSC.vn ;
            faces = spcutMeshSmRSC.f ;
        else
            load(sprintf(QS.fullFileBase.spcutMeshSm, tt), 'spcutMeshSm') 
            vqq = spcutMeshSm.v ;
            nqq = spcutMeshSm.vn ;
            faces = spcutMeshSm.f ;
        end
        
        % % Check against rotated & scaled version of the spcutMeshSm
        % disp('Loading full smoothed mesh timeseries stack')
        % [v3dsmM, nsmM] = QS.loadSPCutMeshSm() ;
        % spcutMeshSm_loaded = true ;
        % % Take this time slice
        % vqq = squeeze(v3dsmM(qq, :, :)) ;
        % nqq = squeeze(nsmM(qq, :, :)) ;
        % % rotate and scale
        % nqqrs = (rot * nqq')' ;
        % vqq = ((rot * vqq')' + trans) * resolution;
        % if flipy
        %     vqq(:, 2) = -vqq(:, 2) ;
        %     % Put the normal direction as inward
        %     nqqrs(:, 1) = -nqqrs(:, 1) ;
        %     nqqrs(:, 3) = -nqqrs(:, 3) ;
        % end
        % [rot, trans] = QS.getRotTrans() ;
        % resolution = QS.APDV.resolution ;
        % flipy = QS.flipy ;

        
    end
    
    % Plot embedding RS colored by normal vector in y direction
    if e0 || overwrite
        fig = figure('visible', 'off') ;
        % Color figure with normals taken via faceNormals in y dir 
        trisurf(faces, vqq(:, 1), vqq(:, 2), vqq(:, 3), ...
            nqq(:, 2), 'EdgeColor', 'none')
        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        title(['Smoothed mesh embedding, t=' num2str(tt)])
        disp(['Saving figure: ' fig1fn])
        saveas(fig, fig1fn)
        close all
    end
    
    % Plot embedding colored by phi
    if e1 || e2 || e3 || overwrite
        fig = figure('visible', 'off') ;
        % Plot embedding in 3d color coding (color=phi)
        pc = (1:nV) .* ones(nU, nV) / nV ;
        trisurf(faces, vqq(:, 1), vqq(:, 2), vqq(:, 3), pc(:), ...
            'EdgeColor', 'none')
        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        c = colorbar() ;
        c.Label.String = '\phi / 2\pi' ;
        title(['Smoothed mesh embedding, t=' num2str(tt)])
        disp(['Saving figure: ' fp0fn])
        saveas(fig, fp0fn)
        view(0, 90)
        saveas(fig, fp1fn)
        view(0, 270)
        saveas(fig, fp2fn)
        view(0, 0)
        saveas(fig, fp3fn)
        view(0, 180)
        saveas(fig, fp4fn)
        close all
    end
end
disp('done')
clearvars fig vqq nqqrs e0 e1 e2 e3 pdir ddir vdir rdir ldir

