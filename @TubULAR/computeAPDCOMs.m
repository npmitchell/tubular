function [acom_sm, pcom_sm, dcom] = computeAPDCOMs(tubi, opts)
%[acom_sm, pcom_sm, dcom] = COMPUTEAPDCOMS(opts)
% Compute the anterior, posterior, and dorsal centers of mass either from:
%   (1) Clicking on the points on the mesh from t=t0
%   (2) moments of inertia of the mesh surface, identifying the endpoints
%       near the long axis of the object at t=t0. (autoAP = true)
%   (3) iLastik training for CENTERLINE computation.  (use_iLastik=true)
%       Note that these are allowed to be different than the
%       APD points for ALIGNMENT computation. 
%       For example, the posterior point
%       might be a point which does NOT form an AP axis  with the 
%       anteriormost point, as in the illustration of the midgut below:
% 
%         P for centerline
%        _x_         Dorsal
%       /  /     ___x_
%      /  /    /      \    Anterior pt for both centerline and for defining APDV axes
%     /  /____/        \ x
%    |  x P for APDV    |
%     \________________/
%    (ventral here, unlabeled)
%
% The default behavior is to use iLastik training if found, unless 
% opts.use_iLastik is set to false. If use_iLastik is false or no iLastik 
% output is found, the default is to have the user click on the points.
% If opts.autoAP == true, then we automatically find A and P positions by
% simply point matching from line intersections onto the mesh.
%
%
% Parameters
% ----------
% opts : struct with fields
%   - use_iLastik : default=true if training h5s are present on disk
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
% apdv_coms_from_training.h5 (rawapdvname, tubi.fileName.apdv)
%   Raw centers of mass for A, P, and D in subsampled pixels, in 
%   probability data space coordinate system
%   Saved to fullfile(meshDir, 'centerline/apdv_coms_from_training.h5')
% tubi.fileName.dcom 
%   txt file with dorsal COM for APDV definition
% rawapdvmatname=fullfile(tubi.dir.cntrline, 'apdv_coms_from_training.mat')
% 
%
% NPMitchell 2020

timePoints = tubi.xp.fileMeta.timePoints ;
apdvoutdir = tubi.dir.cntrline ;
meshDir = tubi.dir.mesh ;
% axorder = tubi.data.axisOrder ; % NOTE: axisorder for texture_axis_order
                                % invoked upon loading IV into tubi.currentData.IV
ilastikOutputAxisOrder = tubi.data.ilastikOutputAxisOrder ;

if isfield(opts, 'aProbFileName')
    aProbFileName = opts.aProbFileName ;
else
    aProbFileName = tubi.fullFileBase.apCenterlineProb ;
end
if isfield(opts, 'pProbFileName')
    pProbFileName = opts.pProbFileName ;
else
    pProbFileName = tubi.fullFileBase.apCenterlineProb ;
end
if isfield(opts, 'ilastikOutputAxisOrder')
    ilastikOutputAxisOrder = opts.ilastikOutputAxisOrder ;
end
if isfield(opts, 'use_iLastik')
    use_iLastik = opts.use_iLastik ;
else
    use_iLastik = exist(aProbFileName, 'file') && ...
        exist(pProbFileName, 'file') ;
    if ~use_iLastik
        disp(['No ilastik training specifically for centerline computation', ...
        'was found, so define the centerline endpoints based on the mesh elongation axis'])
    end
end
if ~use_iLastik
    if isfield(opts, 'autoAP')
        autoAP = opts.autoAP ; % would you like to find automatic points for endcaps A and P?
    else
        autoAP = false ; % instead by default, click on the A and P points at t=t0
    end
end

% Default options
overwrite = false ; 
preview_com = false ;
check_slices = false ;

% Unpack opts
if isfield(opts, 'anteriorChannel')
    anteriorChannel = opts.anteriorChannel ;
else
    anteriorChannel = 1 ;
end
if isfield(opts, 'anteriorChannel')
    posteriorChannel = opts.posteriorChannel ;
else
    posteriorChannel = 2 ;
end
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

rawapdvname = tubi.fileName.apdv ;
rawapdvmatname = fullfile(apdvoutdir, 'apdv_coms_from_training.mat') ;
preview = false ;
if isfield(opts, 'preview')
    preview = opts.preview ;
end


%% Iterate through each mesh to compute acom(t) and pcom(t). Prepare file.
acoms = zeros(length(timePoints), 3) ;
pcoms = zeros(length(timePoints), 3) ;
load_from_disk = false ;
if exist(tubi.fileName.apdv, 'file') && ~overwrite
    load_from_disk = true ;
    try
        h5create(tubi.fileName.apdv, '/acom_sm', size(acoms)) ;
        load_from_disk = false ;
    catch
        try
            acom_sm = h5read(tubi.fileName.apdv, '/acom_sm') ;
            acoms = h5read(tubi.fileName.apdv, '/acom') ;
            disp('acom_sm already exists')
        catch
            load_from_disk = false;
        end
        if load_from_disk
            if size(acoms, 1) ~= length(tubi.xp.fileMeta.timePoints)
                disp(['#timepoints = ' num2str(length(tubi.xp.fileMeta.timePoints)) ])
                disp(['#timepoints on disk = ', num2str(size(acoms, 1))])
                disp(['Must first rename ' tubi.fileName.apdv ' to overwrite with different number of timepoints: moving file.'])
                movefile(tubi.fileName.apdv, [tubi.fileName.apdv '_backup'])
                load_from_disk = false ;
            end
        end
    end
    try
        h5create(tubi.fileName.apdv, '/pcom_sm', size(pcoms)) ;
        load_from_disk = false ;
    catch
        try
            pcom_sm = h5read(tubi.fileName.apdv, '/pcom_sm') ;
            pcoms = h5read(tubi.fileName.apdv, '/pcom') ;
            disp('pcom_sm already exists')
        catch
            load_from_disk = false;
        end
        
        if load_from_disk
            if size(pcoms, 1) ~= length(tubi.xp.fileMeta.timePoints) 
                disp(['Must first rename ' tubi.fileName.apdv ' to overwrite with different number of timepoints: moving file.'])
                movefile(tubi.fileName.apdv, [tubi.fileName.apdv '_backup'])
                load_from_disk = false ;
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
    if use_iLastik
        bad_size = false ;
        if exist(rawapdvmatname, 'file') && ~overwrite
            % load raw data from .mat
            load(rawapdvmatname, 'acoms', 'pcoms')
            bad_size = (size(acoms, 1) ~= length(tubi.xp.fileMeta.timePoints)) ; 
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
                    mesh = read_ply_mod(sprintf(tubi.fullFileBase.mesh, tt)) ;
                    for ii = 1:3
                        subplot(1, 3, ii)
                        trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
                        hold on;
                        plot3(acom(1) * tubi.ssfactor, acom(2) * tubi.ssfactor, acom(3) * tubi.ssfactor, 'o')
                        plot3(pcom(1) * tubi.ssfactor, pcom(2) * tubi.ssfactor, pcom(3) * tubi.ssfactor, 'o')
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
                        dcom = dlmread(tubi.fileName.dcom) ;
                    catch
                        error('Could not load dorsal COM: run tubi.computeAPDVCoords() first')
                    end
                    for ii = 1:3
                        subplot(1, 3, ii)
                        mesh = read_ply_mod(sprintf(tubi.fullFileBase.mesh, tt)) ;
                        trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
                        hold on;
                        plot3(acom(1) * tubi.ssfactor, acom(2) * tubi.ssfactor, acom(3) * tubi.ssfactor, 'o')
                        plot3(pcom(1) * tubi.ssfactor, pcom(2) * tubi.ssfactor, pcom(3) * tubi.ssfactor, 'o')
                        plot3(dcom(1) * tubi.ssfactor, dcom(2) * tubi.ssfactor, dcom(3) * tubi.ssfactor, 'o')
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
                    saveas(gcf, fullfile(tubi.dir.mesh, 'apd_coms_centerline.png'))
                end

            end
            % Save raw data to .mat
            save(rawapdvmatname, 'acoms', 'pcoms')
            clearvars adat pdat
        end
        disp('done determining acoms, pcoms')
    elseif autoAP
        % No ilastik files used to extract probabilities, instead use mesh
        % elongation axis at t0 and then pointmatch for t>t0 and t<t0.
        
        % First do t0
        ssfactor = tubi.ssfactor ;
        tubi.setTime(tubi.t0set()) ;
        meshfn = sprintf(tubi.fullFileBase.mesh, tubi.t0set()) ;
        disp(['Loading mesh ' meshfn])
        mesh = read_ply_mod(meshfn );
        cntrd = mean(mesh.v) ;
        moi = momentOfInertia3D(mesh.v) ; 
        % Ascertain the long axis of the surface at this timepoint
        [eigvect, eigvals] = eig(moi) ;
        [~, ind] = min(abs([eigvals(1,1), eigvals(2,2), eigvals(3,3)])) ;
        % which data axis does this correspond to most closely?
        dotprods = [dot(eigvect(:, ind), [1, 0, 0]), ...
            dot(eigvect(:, ind), [0, 1, 0]), ...
            dot(eigvect(:, ind), [0, 0, 1]) ] ;
        [~, xIndex] = max(abs(dotprods)) ;
        axx = circshift([1,2,3],-xIndex+1) ;
        assert(xIndex == axx(1))
        
        % Place anterior at the intersection of x=xmin and the ray
        % emanating from the centroid along eigenvector
        % Check that the moment of inertia eigvect is not in plane
        % and that they are not parallel. Otherwise use Descartes formula:
        % ax + by + cz + d = 0, where n = [a, b, c] is a vector normal to the plane
        
        % anterior point -- near along the elongated axis
        xval = min(mesh.v(:, xIndex)) ;
        ptInPlaneA = [0,0,0] ;
        ptInPlaneA(xIndex) = xval ;
        
        % posterior point -- far along the elongated axis
        xval = max(mesh.v(:, xIndex)) ;
        ptInPlaneP = [0,0,0] ;
        ptInPlaneP(xIndex) = xval ;
        
        normalToPlane = [0,0,0] ;
        normalToPlane(xIndex) = 1 ; 
        lineDirec = eigvect(:, ind) ;
        [acom, specialCaseA] = linePlaneIntersection(lineDirec, cntrd, ...
            normalToPlane, ptInPlaneA) ;
        [pcom, specialCaseP] = linePlaneIntersection(lineDirec, cntrd, ...
            normalToPlane, ptInPlaneP) ;
        try
            assert(~specialCaseA && ~specialCaseP)
        catch
            disp('Failed to compute automatic endcap points. Is the surface flat?')
        end
        
        % Assignment for t0
        t0 = tubi.t0set() ;
        tidx0 = tubi.xp.tIdx(t0) ;
        acoms(tidx0, :) = acom / ssfactor ;
        pcoms(tidx0, :) = pcom / ssfactor ;
                
        tidxGreater = find(timePoints > t0) ;
        tidxSmaller = find(timePoints < t0) ;
        prevA = acoms(tidx0, :) * ssfactor ;
        prevP = pcoms(tidx0, :) * ssfactor ;
        for tidx = tidxGreater
            tubi.setTime(t0) ;
            mesh = tubi.loadCurrentRawMesh() ;
            [~,idxA] = min(...
                (mesh.v(:,1) - prevA(1)).^2 + ...
                (mesh.v(:,2) - prevA(2)).^2 + ...
                (mesh.v(:,3) - prevA(3)).^2);
            [~,idxP] = min(...
                (mesh.v(:,1) - prevP(1)).^2 + ...
                (mesh.v(:,2) - prevP(2)).^2 + ...
                (mesh.v(:,3) - prevP(3)).^2);
            acom = mesh.v(idxA, :) ;
            pcom = mesh.v(idxP, :) ;
            
            % Check that this is working
            clf
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none') ;
            hold on;
            plot3(prevA(1), prevA(2), prevA(3), 'bo')
            plot3(prevP(1), prevP(2), prevP(3), 'ro')
            plot3(acoms(:, 1), acoms(:, 2), acoms(:,3), '.-')
            plot3(pcoms(:, 1), pcoms(:, 2), pcoms(:,3), '.-')
            
            acoms(tidx, :) = acom / ssfactor ;
            pcoms(tidx, :) = pcom / ssfactor ;
            prevA = acom ; 
            prevP = pcom ;
        end
        
        % Now point match backward in time
        prevA = acoms(tidx0, :) * ssfactor ;
        prevP = pcoms(tidx0, :) * ssfactor ;
        for tidx = fliplr(tidxSmaller)
            tubi.setTime(timePoints(tidx)) ;
            mesh = tubi.loadCurrentRawMesh() ;
            [~,idxA] = min(...
                (mesh.v(:,1) - prevA(1)).^2 + ...
                (mesh.v(:,2) - prevA(2)).^2 + ...
                (mesh.v(:,3) - prevA(3)).^2);
            [~,idxP] = min(...
                (mesh.v(:,1) - prevP(1)).^2 + ...
                (mesh.v(:,2) - prevP(2)).^2 + ...
                (mesh.v(:,3) - prevP(3)).^2);
            acom = mesh.v(idxA, :) ;
            pcom = mesh.v(idxP, :) ;
            acoms(tidx, :) = acom / ssfactor ;
            pcoms(tidx, :) = pcom / ssfactor ;
            prevA = acom ; 
            prevP = pcom ;
        end
    else
        % No ilastik files used to extract probabilities, instead just
        % click on points in 3D on mesh.
        % First do t0
        ssfactor = tubi.ssfactor ;
        tubi.setTime(tubi.t0set()) ;
        meshfn = sprintf(tubi.fullFileBase.mesh, tubi.t0set()) ;
        disp(['Loading mesh ' meshfn])
        mesh = read_ply_mod(meshfn );
        
        %% getpts3d Select points from a 3D scatter plot by clicking on plot
        close all
        h = figure(1);
        clf
        vrs = tubi.xyz2APDV(mesh.v) ;
        trisurf(triangulation(mesh.f, vrs), 'edgecolor', 'none', 'facealpha',0.5)
        axis equal
        xlabel('x'); ylabel('y'); zlabel('z')
        
        % View ANTERIOR endcap 
        msg = 'Rotate the mesh to view Anterior endcap, then press Enter/return';
        disp(msg)
        title(msg)
        key = 'none' ;
        while ~strcmpi(key, 'return')
            waitforbuttonpress;
            key=get(gcf,'CurrentKey');
        end
        datacursormode on
        
        % View ANTERIOR endcap 
        msg = "Select Anterior endcap point, then press 'a' (with Fig in foreground)";
        disp(msg)
        title(msg)
        key = 'none' ;
        while ~strcmpi(key, 'a')
            datacursormode on
            waitforbuttonpress;
            key=get(gcf,'CurrentKey');
        end
        dcm_obj = datacursormode(h);  
        f = getCursorInfo(dcm_obj);
        apt = f.Position ;
        acom = tubi.APDV2xyz(apt) ;
        hold on;
        plot3(apt(1), apt(2), apt(3), 'ro') ;
        legend({'surface', 'anterior'})
        
        %%%%%%%%  
        % Rotate to see posterior
        
        % View ANTERIOR endcap 
        msg = 'Rotate the mesh to view Posterior endcap, then press Enter/return';
        disp(msg)
        title(msg)
        key = 'none' ;
        while ~strcmpi(key, 'return')
            waitforbuttonpress;
            key=get(gcf,'CurrentKey');
        end
        datacursormode on
        msg = "Select the Posterior endcap point, then press 'p' (with Fig in foreground)";
        disp(msg)
        title(msg)
        key = 'none' ;
        while ~strcmpi(key, 'p')
            datacursormode on
            waitforbuttonpress;
            key=get(gcf,'CurrentKey');
        end
        dcm_obj = datacursormode(h);
        f = getCursorInfo(dcm_obj);
        ppt = f.Position ;
        pcom = tubi.APDV2xyz(ppt) ;
        hold on;
        plot3(ppt(1), ppt(2), ppt(3), 'rs') ;
        legend({'surface', 'anterior', 'posterior'})
        
        % Assignment for t0
        t0 = tubi.t0set() ;
        tidx0 = tubi.xp.tIdx(t0) ;
        acoms(tidx0, :) = acom / ssfactor ;
        pcoms(tidx0, :) = pcom / ssfactor ;
                
        tidxGreater = find(timePoints > t0) ;
        tidxSmaller = find(timePoints < t0) ;
        prevA = acoms(tidx0, :) * ssfactor ;
        prevP = pcoms(tidx0, :) * ssfactor ;
        for tidx = tidxGreater
            tubi.setTime(t0) ;
            mesh = tubi.loadCurrentRawMesh() ;
            [~,idxA] = min(...
                (mesh.v(:,1) - prevA(1)).^2 + ...
                (mesh.v(:,2) - prevA(2)).^2 + ...
                (mesh.v(:,3) - prevA(3)).^2);
            [~,idxP] = min(...
                (mesh.v(:,1) - prevP(1)).^2 + ...
                (mesh.v(:,2) - prevP(2)).^2 + ...
                (mesh.v(:,3) - prevP(3)).^2);
            acom = mesh.v(idxA, :) ;
            pcom = mesh.v(idxP, :) ;
            
            % Check that this is working
            clf
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none') ;
            hold on;
            plot3(prevA(1), prevA(2), prevA(3), 'bo')
            plot3(prevP(1), prevP(2), prevP(3), 'ro')
            plot3(acoms(:, 1), acoms(:, 2), acoms(:,3), '.-')
            plot3(pcoms(:, 1), pcoms(:, 2), pcoms(:,3), '.-')
            
            acoms(tidx, :) = acom / ssfactor ;
            pcoms(tidx, :) = pcom / ssfactor ;
            prevA = acom ; 
            prevP = pcom ;
        end
        
        % Now point match backward in time
        prevA = acoms(tidx0, :) * ssfactor ;
        prevP = pcoms(tidx0, :) * ssfactor ;
        for tidx = fliplr(tidxSmaller)
            tubi.setTime(timePoints(tidx)) ;
            mesh = tubi.loadCurrentRawMesh() ;
            [~,idxA] = min(...
                (mesh.v(:,1) - prevA(1)).^2 + ...
                (mesh.v(:,2) - prevA(2)).^2 + ...
                (mesh.v(:,3) - prevA(3)).^2);
            [~,idxP] = min(...
                (mesh.v(:,1) - prevP(1)).^2 + ...
                (mesh.v(:,2) - prevP(2)).^2 + ...
                (mesh.v(:,3) - prevP(3)).^2);
            acom = mesh.v(idxA, :) ;
            pcom = mesh.v(idxP, :) ;
            acoms(tidx, :) = acom / ssfactor ;
            pcoms(tidx, :) = pcom / ssfactor ;
            prevA = acom ; 
            prevP = pcom ;
        end
                
    end
    
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
            if smfrac > 1 
                smfrac = 1 ;
            end
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
    dcom = dlmread(tubi.fileName.dcom) ;
    % [xyzlim, ~, ~, ~] = tubi.getXYZLims() ;
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