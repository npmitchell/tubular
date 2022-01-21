function [acom_sm, pcom_sm, dcom] = computeAPDCOMs(QS, opts)
%[acom_sm, pcom_sm, dcom] = COMPUTEAPDCOMS(opts)
% Compute the anterior, posterior, and dorsal centers of mass from training
% for CENTERLINE computation. Note that these can be different than the
% APD points for ALIGNMENT computation. In particular, the posterior point
% might be a point which does NOT form an AP axis  with the anteriormost 
% point, as in the illustration of the midgut below:
% 
%         P for centerline
%        _x_         Dorsal
%       /  /     ___x_
%      /  /    /      \  Anterior
%     /  /____/        \x
%    |  xP for APDV     |
%     \________________/
%    (ventral here, unlabeled)
%
%
% Parameters
% ----------
% opts : struct with fields
%   - timePoints
%   - dorsal_thres : float between 0 and 1
%   - anteriorChannel : int
%   - posteriorChannel
%   - dorsalChannel
%   - overwrite : bool
%   - apdvoutdir : string
%   - meshDir : string
%   - preview_com : bool
%   - check_slices : bool
%   - axorder : length 3 int array
%   - smwindow : float or int (optional, default=30)
%       number of timepoints over which we smooth
%   - preview : bool (optional, default=false)
%       
%
% OUTPUTS
% -------
% apdv_coms_from_training.h5 (rawapdvname, QS.fileName.apdv)
%   Raw centers of mass for A, P, and D in subsampled pixels, in 
%   probability data space coordinate system
%   Saved to fullfile(meshDir, 'centerline/apdv_coms_from_training.h5')
% QS.fileName.dcom 
%   txt file with dorsal COM for APDV definition
% rawapdvmatname=fullfile(QS.dir.cntrline, 'apdv_coms_from_training.mat')
% 
%
% NPMitchell 2020

timePoints = QS.xp.fileMeta.timePoints ;
apdvoutdir = QS.dir.cntrline ;
meshDir = QS.dir.mesh ;
% axorder = QS.data.axisOrder ; % NOTE: axisorder for texture_axis_order
                                % invoked upon loading IV into QS.currentData.IV
ilastikOutputAxisOrder = QS.data.ilastikOutputAxisOrder ;

if isfield(opts, 'aProbFileName')
    aProbFileName = opts.aProbFileName ;
else
    aProbFileName = QS.fullFileBase.apCenterlineProb ;
end
if isfield(opts, 'pProbFileName')
    pProbFileName = opts.pProbFileName ;
else
    pProbFileName = QS.fullFileBase.apCenterlineProb ;
end
if isfield(opts, 'ilastikOutputAxisOrder')
    ilastikOutputAxisOrder = opts.ilastikOutputAxisOrder ;
end

% Default options
overwrite = false ; 
preview_com = false ;
check_slices = false ;

% Unpack opts
anteriorChannel = opts.anteriorChannel ;
posteriorChannel = opts.posteriorChannel ;
if isfield(opts, 'overwrite')
    overwrite = opts.overwrite ;
end
if isfield(opts, 'preview_com')
    preview_com = opts.preview_com ;
end
if isfield(opts, 'check_slices')
    check_slices = opts.check_slices ;
end

% Default valued options
smwindow = 30 ;
if isfield(opts, 'smwindow')
    smwindow = opts.smwindow ;
end

rawapdvname = QS.fileName.apdv ;
rawapdvmatname = fullfile(apdvoutdir, 'apdv_coms_from_training.mat') ;
preview = false ;
if isfield(opts, 'preview')
    preview = opts.preview ;
end


%% Iterate through each mesh to compute acom(t) and pcom(t). Prepare file.
acoms = zeros(length(timePoints), 3) ;
pcoms = zeros(length(timePoints), 3) ;
load_from_disk = false ;
if exist(QS.fileName.apdv, 'file') && ~overwrite
    load_from_disk = true ;
    try
        h5create(QS.fileName.apdv, '/acom_sm', size(acoms)) ;
        load_from_disk = false ;
    catch
        try
            acom_sm = h5read(QS.fileName.apdv, '/acom_sm') ;
            acoms = h5read(QS.fileName.apdv, '/acom') ;
            disp('acom_sm already exists')
        catch
            load_from_disk = false;
        end
        if load_from_disk
            if size(acoms, 1) ~= length(QS.xp.fileMeta.timePoints)
                disp(['#timepoints = ' num2str(length(QS.xp.fileMeta.timePoints)) ])
                disp(['#timepoints on disk = ', num2str(size(acoms, 1))])
                error(['Must first delete ' QS.fileName.apdv ' to overwrite with different number of timepoints'])
            end
        end
    end
    try
        h5create(QS.fileName.apdv, '/pcom_sm', size(pcoms)) ;
        load_from_disk = false ;
    catch
        try
            pcom_sm = h5read(QS.fileName.apdv, '/pcom_sm') ;
            pcoms = h5read(QS.fileName.apdv, '/pcom') ;
            disp('pcom_sm already exists')
        catch
            load_from_disk = false;
        end
        
        if load_from_disk
            if size(pcoms, 1) ~= length(QS.xp.fileMeta.timePoints) 
                error(['Must first delete ' QS.fileName.apdv ' to overwrite with different number of timepoints'])
            end
        end
    end
end
if ~load_from_disk
    disp('acom and/or pcom not already saved on disk. Compute them')
end

disp(['Load from disk? =>', num2str(load_from_disk)])

%% Compute smoothed acom and pcom if not loaded from disk -- RAW XYZ coords
if ~load_from_disk || overwrite
    % Compute raw acom and pcom if not loaded from disk -- RAW XYZ coords
    bad_size = false ;
    if exist(rawapdvmatname, 'file') && ~overwrite
        % load raw data from .mat
        load(rawapdvmatname, 'acoms', 'pcoms')
        bad_size = (size(acoms, 1) ~= length(QS.xp.fileMeta.timePoints)) ; 
    end
    
    if ~exist(rawapdvmatname, 'file') || overwrite || bad_size
        for tidx = 1:length(timePoints)
            tt = timePoints(tidx) ;
            %% Load the AP axis determination
            msg = ['Computing acom, pcom for ' num2str(tt) ] ;
            disp(msg)
            thres = 0.5 ;
            % load the probabilities for anterior posterior dorsal
            afn = sprintf(aProbFileName, tt);
            pfn = sprintf(pProbFileName, tt);
            disp(['Reading ', afn])
            adatM = h5read(afn, '/exported_data');
            if ~strcmp(afn, pfn)
                disp(['Reading ' pfn])
                pdatM = h5read(pfn, '/exported_data') ;
            else
                pdatM = adatM ;
            end

            % Load the training for anterior and posterior positions
            if strcmpi(ilastikOutputAxisOrder(1), 'c')
                adat = squeeze(adatM(anteriorChannel,:,:,:)) ;
                pdat = squeeze(pdatM(posteriorChannel,:,:,:)) ;
            elseif strcmpi(ilastikOutputAxisOrder(4), 'c')
                adat = squeeze(adatM(:,:,:,anteriorChannel)) ;
                pdat = squeeze(pdatM(:,:,:,posteriorChannel)) ;
            else
                error('Did not recognize ilastikAxisOrder. Code here')
            end
            
            if contains(lower(ilastikOutputAxisOrder), 'xyz')
                disp('no permutation necessary')
                % adat = permute(adat, [1,2,3]);
            elseif contains(lower(ilastikOutputAxisOrder), 'yxz')
                disp('permuting yxz')
                adat = permute(adat, [2,1,3]);
                pdat = permute(pdat, [2,1,3]);
            elseif contains(lower(ilastikOutputAxisOrder), 'zyx')
                disp('permuting zyx')
                adat = permute(adat, [3,2,1]);
                pdat = permute(pdat, [3,2,1]);
            elseif contains(lower(ilastikOutputAxisOrder), 'yzx')
                disp('permuting zyx')
                adat = permute(adat, [2,3,1]);
                pdat = permute(pdat, [2,3,1]);
            else
                error('unrecognized permutation in ilastikOutputAxisOrder')
            end
            
            % define axis order: 
            % if 1, 2, 3: axes will be yxz
            % if 1, 3, 2: axes will be yzx
            % if 2, 1, 3: axes will be xyz (ie first second third axes, ie --> 
            % so that bright spot at im(1,2,3) gives com=[1,2,3]
            
            % Note that axisOrder is applying upon invoking getCurrentData()
            % adat = permute(adat, axorder) ;
            % pdat = permute(pdat, axorder) ;
            
            options.check = preview_com ;
            disp('Extracting acom')
            options.color = 'red' ;
            
            acom = com_region(adat, thres, options) ;
            disp('Extracting pcom')
            options.color = 'blue' ;
            pcom = com_region(pdat, thres, options) ;
            clearvars options
            % [~, acom] = match_training_to_vertex(adat, thres, vertices, options) ;
            % [~, pcom] = match_training_to_vertex(pdat, thres, vertices, options) ;
            acoms(tidx, :) = acom ;
            pcoms(tidx, :) = pcom ;
            if preview
                disp('acom = ')
                acoms(tidx, :)
                disp('pcom = ')
                pcoms(tidx, :)
                
                clf
                mesh = read_ply_mod(sprintf(QS.fullFileBase.mesh, tt)) ;
                for ii = 1:3
                    subplot(1, 3, ii)
                    trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
                    hold on;
                    plot3(acom(1) * QS.ssfactor, acom(2) * QS.ssfactor, acom(3) * QS.ssfactor, 'o')
                    plot3(pcom(1) * QS.ssfactor, pcom(2) * QS.ssfactor, pcom(3) * QS.ssfactor, 'o')
                    axis equal
                    if ii == 1
                        view(0, 90)
                    elseif ii == 2
                        view(90, 0)
                    else
                        view(180, 0)
                    end
                end
                sgtitle(['t = ' num2str(tt)])
                pause(0.001)
            end
            
            % PLOT APD points on mesh
            if tidx == 1
                % load current mesh & plot the dorsal dot
                clf
                try
                    dcom = dlmread(QS.fileName.dcom) ;
                catch
                    error('Could not load dorsal COM: run QS.computeAPDVCoords() first')
                end
                for ii = 1:3
                    subplot(1, 3, ii)
                    mesh = read_ply_mod(sprintf(QS.fullFileBase.mesh, tt)) ;
                    trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
                    hold on;
                    plot3(acom(1) * QS.ssfactor, acom(2) * QS.ssfactor, acom(3) * QS.ssfactor, 'o')
                    plot3(pcom(1) * QS.ssfactor, pcom(2) * QS.ssfactor, pcom(3) * QS.ssfactor, 'o')
                    plot3(dcom(1) * QS.ssfactor, dcom(2) * QS.ssfactor, dcom(3) * QS.ssfactor, 'o')
                    axis equal
                    if ii == 1
                        view(0, 90)
                    elseif ii == 2
                        view(90, 0)
                    else
                        view(180, 0)
                    end
                end
                sgtitle('APD COMs for APD COMs for centerline')
                saveas(gcf, fullfile(QS.dir.mesh, 'apd_coms_centerline.png'))
            end
            
        end
        % Save raw data to .mat
        save(rawapdvmatname, 'acoms', 'pcoms')
        clearvars adat pdat
    end
    disp('done determining acoms, pcoms')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Smooth the acom and pcom data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(timePoints) > 2 
        if smwindow > 0
            disp('Smoothing acom and pcom...')
            acom_sm = 0 * acoms ;
            pcom_sm = 0 * acoms ;
            % fraction of data for smoothing window
            smfrac = smwindow / double(length(timePoints)) ;  
            acom_sm(:, 1) = smooth(timePoints, acoms(:, 1), smfrac, 'rloess');
            pcom_sm(:, 1) = smooth(timePoints, pcoms(:, 1), smfrac, 'rloess');
            acom_sm(:, 2) = smooth(timePoints, acoms(:, 2), smfrac, 'rloess');
            pcom_sm(:, 2) = smooth(timePoints, pcoms(:, 2), smfrac, 'rloess');
            acom_sm(:, 3) = smooth(timePoints, acoms(:, 3), smfrac, 'rloess');
            pcom_sm(:, 3) = smooth(timePoints, pcoms(:, 3), smfrac, 'rloess');
        else
            disp('No smoothing to acom and pcom...')
            acom_sm = acoms ;
            pcom_sm = pcoms ;
        end
    else
        acom_sm = acoms ;
        pcom_sm = pcoms ;
    end
    
    if preview
        plot(timePoints, acoms - mean(acoms,1), '.')
        hold on
        plot(timePoints, acom_sm - mean(acoms, 1), '-')
        sgtitle('Smoothed COMs for AP')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save smoothed anterior and posterior centers of mass ===============
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        h5create(rawapdvname, '/acom', size(acoms)) ;
    catch
        disp('acom already exists as h5 file. Overwriting.')
    end
    try
        h5create(rawapdvname, '/pcom', size(pcoms)) ;
    catch
        disp('pcom already exists as h5 file. Overwriting.')
    end
    try
        h5create(rawapdvname, '/acom_sm', size(acom_sm)) ;
    catch
        disp('acom_sm already exists as h5 file. Overwriting.')
    end
    try
        h5create(rawapdvname, '/pcom_sm', size(pcom_sm)) ;
    catch
        disp('pcom_sm already exists as h5 file. Overwriting.')
    end
    h5write(rawapdvname, '/acom', acoms) ;
    h5write(rawapdvname, '/pcom', pcoms) ;
    h5write(rawapdvname, '/acom_sm', acom_sm) ;
    h5write(rawapdvname, '/pcom_sm', pcom_sm) ;
else
    disp('Skipping, since already loaded acom_sm and pcom_sm')
    if preview
        acom_sm = h5read(rawapdvname, '/acom_sm');
        pcom_sm = h5read(rawapdvname, '/pcom_sm');
        plot3(acom_sm(:, 1), acom_sm(:, 2), acom_sm(:, 3))
        hold on;
        plot3(pcom_sm(:, 1), pcom_sm(:, 2), pcom_sm(:, 3))
        xlabel('x [subsampled pix]')
        ylabel('y [subsampled pix]')
        zlabel('z [subsampled pix]')
        axis equal
        pause(1)
    end
end

disp('done with AP COMs')

%% Display APDV COMS over time
try
    dcom = dlmread(QS.fileName.dcom) ;
    [xyzlim, ~, ~, ~] = QS.getXYZLims() ;
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        % Plot the APDV points
        clf
        plot3(acom_sm(tidx, 1), acom_sm(tidx, 2), acom_sm(tidx, 3), 'ro')
        hold on;
        plot3(acoms(tidx, 1), acoms(tidx, 2), acoms(tidx, 3), 'r.')
        plot3(pcom_sm(tidx, 1), pcom_sm(tidx, 2), pcom_sm(tidx, 3), 'b^')
        plot3(pcoms(tidx, 1), pcoms(tidx, 2), pcoms(tidx, 3), 'b.')
        plot3(dcom(1, 1), dcom(1, 2), dcom(1, 3), 'cs')
        axis equal
        title(['t = ', num2str(tp)]) 
        pause(0.01)
    end
catch
    disp('Could not display aligned meshes -- does dcom exist on file?')
end

disp('done')