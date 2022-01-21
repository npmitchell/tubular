function [acom,pcom,dcom, rot,trans] = computeAPDVCoords(QS, opts)
%[acom,pcom,dcom, rot,trans] = computeAPDVCoords(QS, opts)
% Compute the APDV coordinate system, defined WRT the data coordSys by a
% rotation, translation, and resolution. This is done by identifying a 
% dorsal point, so that z dim in APDV points from the AP axis to dorsal
% along the shortest linesegment emanating from the AP axis to the dorsal
% point. 
%
% Chirality: note that the probabilties field is assumed to be meshlike so
% the first two axis are swapped within com_region(). 
%
% Todo: back-save acom_for_rot and pcom_for_rot for datasets already
% processed, inferred from rot and trans
%
% Note that axisOrder is applying upon invoking getCurrentData()


%% Default options
anteriorMethod = 'InsideOrNearestVertex' ;  % default is to use COM if COM 
                            % is inside the mesh polyhedron, otherwise 
                            % project onto nearest mes vertex
normal_step = 1/QS.ssfactor ;
% Note that axisOrder is applying upon invoking getCurrentData()
% axorder = QS.data.axisOrder ;
ilastikOutputAxisOrder = QS.data.ilastikOutputAxisOrder ;
thres = 0.5 ;
ssfactor = QS.ssfactor ;

%% Unpack opts
if isfield(opts, 'aProbFileName')
    aProbFileName = opts.aProbFileName ;
else
    aProbFileName = QS.fullFileBase.apdProb ;
end
if isfield(opts, 'pProbFileName')
    pProbFileName = opts.pProbFileName ;
else
    pProbFileName = QS.fullFileBase.apdProb ;
end
if isfield(opts, 'dProbFileName')
    dProbFileName = opts.dProbFileName ;
else
    dProbFileName = QS.fullFileBase.apdProb ;
end
if isfield(opts, 'ilastikOutputAxisOrder')
    ilastikOutputAxisOrder = opts.ilastikOutputAxisOrder ;
end
if isfield(opts, 'normal_step')
    normal_step = opts.normal_step ;
end
if isfield(opts, 'thres')
    thres = opts.thres ;
end
% Which timepoint to use to define dorsal and AP axis?
if isfield(opts, 'tref')
    trefIDx = find(QS.xp.fileMeta.timePoints == opts.tref) ;
    tt = opts.tref ;
    if isempty(trefIDx)
        error(['could not find match trefIDx to opts.tref=' num2str(opts.tref)])
    end
else
    % by default use first timepoint
    trefIDx = 1 ;
    tt = QS.xp.fileMeta.timePoints(trefIDx) ;
end

% Default options
overwrite = false ; 
preview = false ;
check_slices = false ;

% Unpack opts
dorsal_thres = opts.dorsal_thres ;
anteriorChannel = opts.anteriorChannel ;
posteriorChannel = opts.posteriorChannel ;
dorsalChannel = opts.dorsalChannel ;
if isfield(opts, 'overwrite')
    overwrite = opts.overwrite ;
end
if isfield(opts, 'preview')
    preview = opts.preview ;
end
if isfield(opts, 'check_slices')
    check_slices = opts.check_slices ;
end

%% Prepare filenames
rotname = QS.fileName.rot ;
transname = QS.fileName.trans ;
dcomname = fullfile(QS.dir.mesh, 'dcom_for_rot.txt') ;  % QS.fileName.dcom
acomname = fullfile(QS.dir.mesh, 'acom_for_rot.txt') ;
pcomname = fullfile(QS.dir.mesh, 'pcom_for_rot.txt') ;

% Check if rotation and translation exist on disk
no_rot_on_disk = ~exist(rotname, 'file') ;
no_trans_on_disk = ~exist(transname, 'file') ;
redo_rot_calc = no_rot_on_disk || no_trans_on_disk || overwrite ;
if exist(rotname, 'file') 
    disp('rot exists on file')
else
    disp(['no rot txt file: ' rotname ])
end

if redo_rot_calc || overwrite
    %% Dorsal COM for first timepoint
    % load the probabilities for anterior posterior dorsal
    afn = sprintf(aProbFileName, tt);
    pfn = sprintf(pProbFileName, tt);
    dfn = sprintf(dProbFileName, tt);
    disp(['Reading ' afn])
    adatM = h5read(afn, '/exported_data');

    % Obtain posterior data probabilities file
    if ~strcmp(afn, pfn)
        pdatM = h5read(pfn, '/exported_data') ;
    else
        pdatM = adatM ;
    end

    % Obtain dorsal data probabilities file
    if strcmp(afn, dfn) 
        ddatM = adatM ;
    elseif strcmp(pfn, dfn)
        ddatM = pdatM ;
    else
        ddatM = h5read(dfn, '/exported_data') ;
    end

    % Obtain just the dorsal channel from probabilities
    if strcmpi(ilastikOutputAxisOrder, 'cxyz')
        ddat = squeeze(ddatM(dorsalChannel, :, :, :)) ;
    elseif strcmpi(ilastikOutputAxisOrder, 'xyzc')
        ddat = squeeze(ddatM(:, :, :, dorsalChannel)) ;
    elseif strcmpi(ilastikOutputAxisOrder, 'czyx')
        ddat = squeeze(ddatM(dorsalChannel, :, :, :)) ;
        ddat = permute(ddat, [3, 2, 1]) ;
    elseif strcmpi(ilastikOutputAxisOrder, 'cyxz')
        ddat = squeeze(ddatM(dorsalChannel, :, :, :)) ;
        ddat = permute(ddat, [2, 1, 3]) ;
    elseif strcmpi(ilastikOutputAxisOrder, 'zyxc')
        ddat = squeeze(ddatM(:, :, :, dorsalChannel)) ;
        ddat = permute(ddat, [3, 2, 1]) ;
    elseif strcmpi(ilastikOutputAxisOrder, 'yxzc')
        ddat = squeeze(ddatM(:, :, :, dorsalChannel)) ;
        ddat = permute(ddat, [2, 1, 3]) ;
    elseif strcmpi(ilastikOutputAxisOrder, 'czxy')
        ddat = squeeze(ddatM(dorsalChannel, :, :, :)) ;
        ddat = permute(ddat, [3, 1, 2]) ;
    else
        error('Did not recognize ilastikOutputAxisOrder')
    end
    
    % Convert from xyz to the mesh / QS class axis order, if not xyz
    % Note that axisOrder is applying upon invoking getCurrentData()
    % ddat = permute(ddat, axorder) ;

    options.check = preview ; 
    options.check_slices = check_slices ; 
    options.color = 'green' ;

    apcomsOK = exist(acomname, 'file') && exist(pcomname, 'file') ;

    % Load dcom if already on disk
    if exist(dcomname, 'file') && ~overwrite && apcomsOK
        disp('Loading dorsal COM from disk')
        dcom = dlmread(dcomname) ;
        startpt = dlmread(acomname) ;
        spt = startpt * ssfactor;
        pcom = dlmread(pcomname) ;
    else
        if exist(dcomname, 'file')
            disp('Overwriting existing dorsal COM on disk')
        else
            disp('Computing dorsal COM for the first time')
        end
        cont = input('APD COMS not on disk -- Compute them? [N/y]', 's') ;
        if ~contains(lower(cont), 'y')
            error('exiting to stop APD COM computation')
        end

        search4com = true ;
        % start with a threshold == dorsal_thres, iteratively lower if
        % necessary
        tmp_dorsal_thres = dorsal_thres ;
        while search4com 
            try
                msg = 'Finding com region of dorsal data of size ' ;
                msg = [msg '(' num2str(size(ddat, 1)) ', ' ] ;
                msg = [msg num2str(size(ddat, 2)) ', '] ;
                msg = [msg num2str(size(ddat, 3)) ')'] ;
                disp([msg 'thres=' num2str(tmp_dorsal_thres)])
                close all
                dcom = com_region(ddat, tmp_dorsal_thres, options) ;
                search4com = false ;
            catch
                msg = 'Could not find any dorsal signal within 1% probability' ;
                disp(msg)
                disp('Showing dorsal signal volume, leaf by leaf')
                clf; set(gcf, 'visible', 'on')
                for qq=1:size(ddat, 1)
                    imshow(squeeze(ddat(qq, :, :)))
                    title(['dorsal signal, x=' num2str(qq)])
                    pause(0.01)                            
                end
                error(msg)
            end
        end

        % load current mesh & plot the dorsal dot
        for ii = 1:3
            subplot(1, 3, ii)
            mesh = read_ply_mod(sprintf(QS.fullFileBase.mesh, tt)) ;
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            plot3(dcom(1) * ssfactor, dcom(2) * ssfactor, dcom(3) * ssfactor, 'o')
            axis equal
            if ii == 1
                view(0, 90)
            elseif ii == 2
                view(90, 0)
            else
                view(180, 0)
            end
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        sgtitle('Dorsal COM for APDV coordinates')
        saveas(gcf, [QS.fileName.dcom(1:end-3) 'png'])

        % SAVE DCOM
        dlmwrite(dcomname, dcom) ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Compute and save anterior and posterior points to use as definition of 
        % AP axis and for translation offset to put anterior at the origin
        if strcmpi(ilastikOutputAxisOrder(1), 'c')
            adat = squeeze(adatM(anteriorChannel,:,:,:)) ;
            pdat = squeeze(pdatM(posteriorChannel,:,:,:)) ;
        elseif strcmpi(ilastikOutputAxisOrder(4), 'c')
            adat = squeeze(adatM(:,:,:,anteriorChannel)) ;
            pdat = squeeze(pdatM(:,:,:,posteriorChannel)) ;
        else
            error('Did not recognize ilastikAxisOrder. Code here')
        end

        % define axis order: 
        % if 1, 2, 3: axes will be yxz
        % if 1, 3, 2: axes will be yzx
        % if 2, 1, 3: axes will be xyz (ie first second third axes, ie --> 
        % so that bright spot at im(1,2,3) gives com=[1,2,3]
        % Note that axisOrder is applying upon invoking getCurrentData()
        % adat = permute(adat, axorder) ;
        % pdat = permute(pdat, axorder) ;

        options.check = preview ;
        disp('Extracting acom')
        options.color = 'red' ;
        acom = com_region(adat, thres, options) ;
        disp('Extracting pcom')
        options.color = 'blue' ;
        pcom = com_region(pdat, thres, options) ;
        clearvars options
        % [~, acom] = match_training_to_vertex(adat, thres, vertices, options) ;
        % [~, pcom] = match_training_to_vertex(pdat, thres, vertices, options) ;
        
        xyzstring = erase(lower(ilastikOutputAxisOrder), 'c') ;
        xpos = strfind(xyzstring, 'x') ;
        ypos = strfind(xyzstring, 'y') ;
        zpos = strfind(xyzstring, 'z') ;
        if isempty(xpos) || isempty(ypos) || isempty(zpos)
            error(['did not recognize ilastikOutputAxisOrder=' ilastikOutputAxisOrder])
        end
        acomPermuted = [0, 0, 0] ;
        acomPermuted(xpos) = acom(1) ;
        acomPermuted(ypos) = acom(2) ;
        acomPermuted(zpos) = acom(3) ;
        acom = acomPermuted ;
        pcomPermuted = [0, 0, 0] ;
        pcomPermuted(xpos) = pcom(1) ;
        pcomPermuted(ypos) = pcom(2) ;
        pcomPermuted(zpos) = pcom(3) ;
        pcom = pcomPermuted ;
        % switch lower(xyzstring)
        %     case 'xyz'
        %         disp('no permuting necessary')
        %     case 'zyx'
        %         acom = [acom(3) acom(2) acom(1)];
        %         pcom = [pcom(3) pcom(2) pcom(1)];
        %     case 'yxz'
        %         acom = [acom(2) acom(1) acom(3)];
        %         pcom = [pcom(2) pcom(1) pcom(3)];
        %     case 'z'
        %         acom = [acom(2) acom(1) acom(3)];
        %         pcom = [pcom(2) pcom(1) pcom(3)];
        %     otherwise
        %         error('did not recognize ilastikOutputAxisOrder')   
        % end
        
        if preview
            disp('acom = ')
            acom
            disp('pcom = ')
            pcom
        end

        % PLOT APD points on mesh
        % load current mesh & plot the dorsal dot
        clf
        for ii = 1:3
            subplot(1, 3, ii)
            mesh = read_ply_mod(sprintf(QS.fullFileBase.mesh, tt)) ;
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            plot3(acom(1) * ssfactor, acom(2) * ssfactor, acom(3) * ssfactor, 'o')
            plot3(pcom(1) * ssfactor, pcom(2) * ssfactor, pcom(3) * ssfactor, 'o')
            plot3(dcom(1) * ssfactor, dcom(2) * ssfactor, dcom(3) * ssfactor, 'o')
            axis equal
            if ii == 1
                view(0, 90)
            elseif ii == 2
                view(90, 0)
            else
                view(180, 0)
            end
        end
        sgtitle('APD COMs for APDV coordinates')
        saveas(gcf, fullfile(QS.dir.mesh, 'apd_coms.png'))

        %% Define "start point" as anterior COM projected onto mesh
        if strcmpi(anteriorMethod, 'insideornearestvertex')
            disp('Defining start point via point matching for first TP')
            meshfn = sprintf(QS.fullFileBase.mesh, QS.xp.fileMeta.timePoints(1)) ;
            disp(['Loading mesh ' meshfn])
            mesh = read_ply_mod(meshfn );
            vtx_sub = mesh.v / ssfactor ;
            vn = mesh.vn ;
            fvsub = struct('faces', mesh.f, 'vertices', vtx_sub, 'normals', vn) ;
            % Check if acom is inside mesh. If so, use that as starting point.
            ainside = inpolyhedron(fvsub, acom(1), acom(2), acom(3)) ;

            if ainside
                disp('start point for centerline is inside mesh')
                startpt = acom' ;
            else
                % Point match for aind and pind
                disp(['Point matching mesh ' meshfn])
                adist2 = sum((vtx_sub - acom) .^ 2, 2);
                %find the smallest distance and use that as an index 
                aind = find(adist2 == min(adist2)) ;

                % move along the inward normal of the mesh from the matched vertex
                vtx = [vtx_sub(aind, 1), vtx_sub(aind, 2), vtx_sub(aind, 3)]' ;
                normal = fvsub.normals(aind, :) ;
                startpt = vtx(:) + normal(:);
                if ~inpolyhedron(fvsub, startpt(1), startpt(2), startpt(3)) 
                    % this didn't work, check point in reverse direction
                    startpt = vtx(:) - normal(:) * normal_step ;
                    if ~inpolyhedron(fvsub, startpt(1), startpt(2), startpt(3))
                        % Can't seem to jitter into the mesh, so use vertex
                        disp("Can't seem to jitter into the mesh, so using vertex for startpt")
                        startpt = vtx ;
                    end
                end
            end 
            % Note: Keep startpt in subsampled units
            acom = [startpt(1), startpt(2), startpt(3)] ;
            spt = startpt * ssfactor;
        elseif strcmpi(anteriorMethod, 'com')
            spt = acom * ssfactor ;
        else
            error(['Did not recognize anteriorMethod = ', anteriorMethod])
        end

        % Save to disk
        dlmwrite(acomname, acom) ;
        dlmwrite(pcomname, pcom) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute rotation to align AP axis with X axis
    disp('Obtaining dorsal direction from first TP')
    if overwrite && ~no_rot_on_disk
        disp('Overwriting rot calculation using dorsal pt')
    elseif redo_rot_calc
        disp('Computing rot calculation using dorsal pt for the first time')
    end

    if redo_rot_calc
        % load the probabilities for anterior posterior dorsal
        % Load dcom if already on disk
        disp(['Loading dorsal COM from disk: ' dcomname])
        dcom = dlmread(dcomname) ;

        % compute rotation -- Note we choose to use startpt instead of
        % acom here so that the dorsal point will lie in the y=0 plane.
        % Explanation: we will subtract off startpt * ssfactor as the
        % translation, so if this differs from acom in y dim, then
        % dorsal point gets shifted in y.
        origin = startpt ; % [startpt(1), startpt(2), startpt(3)] ;
        apaxis = pcom(:) - origin(:) ;
        aphat = reshape(apaxis(:) / norm(apaxis), [1, 3]) ;

        % SEE rotate3dToAlignAxis()
        % compute rotation matrix using this procedure: 
        % https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        xhat = [1, 0, 0] ;
        zhat = [0, 0, 1] ;
        ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
        RU = @(A,B) eye(3) + ssc(cross(A,B)) + ...
             ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2) ;
        % rotz aligns AP to xhat (x axis)
        rotx = RU(aphat(:), xhat(:)) ;

        % Rotate dorsal to the z axis
        % find component of dorsal vector from acom perpendicular to AP
        dcom = reshape(dcom, [1, 3]) ;
        origin = reshape(origin, [1, 3]) ;
        % Get perpendicular component of the dorsal vector wrt AP axis
        dvec = rotx * (dcom - origin)' - rotx * (dot(dcom - origin, aphat) * aphat)' ;
        dhat = dvec / norm(dvec) ;
        rotz = RU(dhat, zhat) ;
        rot = rotz * rotx  ;

        % % test arrow
        % aphx = rotx * aphat' ;
        % aphxz = rot * aphat' ;
        % vecs = [aphat; aphx'; aphxz'] ;
        % for qq = 1:length(vecs)
        %     plot3([0, vecs(qq, 1)], [0, vecs(qq, 2)], [0, vecs(qq, 3)], '.-') 
        %     hold on;
        % end
        % axis equal
        % t2 = dcom - origin ;
        % aphx = rotx * t2' ;
        % aphxz = rot * t2' ;
        % vecs = [t2; aphx'; aphxz'] ;
        % for qq = 1:length(vecs)
        %     plot3([0, vecs(qq, 1)], [0, vecs(qq, 2)], [0, vecs(qq, 3)], '.-') 
        %     hold on;
        % end
        % legend({'original', 'rotx', 'rot', 'dorsal', 'rotxd', 'rotd'})
        % axis equal
        % error('here')

        % Save the rotation matrix
        disp(['Saving rotation matrix to txt: ', rotname])
        dlmwrite(rotname, rot)
    else
        disp('Loading rot from disk...')
        rot = dlmread(rotname) ;
    end

    %% Compute the translation to put anterior to origin AFTER rot & scale
    if overwrite || ~exist(transname, 'file')
        % Save translation in units of mesh coordinates
        spt = reshape(spt, [1, 3]) ;
        trans = -(rot * spt')' ;
        disp(['Saving translation vector (post rotation) to txt: ', transname])
        dlmwrite(transname, trans)
    else
        trans = dlmread(transname, ',');
    end
else
    disp('Rot and trans already on disk')
    rot = dlmread(rotname) ;
    trans = dlmread(transname, ',');
    dcom = dlmread(dcomname) ;
    acom = dlmread(acomname) ;
    pcom = dlmread(pcomname) ;
end


%%%%%%%%%%%%%%%%%%%%%%
% % disp('Showing dorsal segmentation...')
% clf
% for slice=1:2:size(ddat, 2)
%     im = squeeze(ddat(:, slice, :)) ;
%     % im(im < dorsal_thres) = 0 ;
%     imshow(im)
%     xlabel('x')
%     ylabel('z')
%     hold on
%     plot(dcom(:, 1), dcom(:, 3), 'o')
%     title([num2str(slice) '/' num2str(size(apdat, 3))])
%     pause(0.001)
% end
%%%%%%%%%%%%%%%%%%%%%%
fig = figure ;
disp('Displaying mesh in figure ...')
mesh = read_ply_mod(sprintf(QS.fullFileBase.mesh, tt)) ;
for ii = 1:3
    subplot(1, 3, ii)
    trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
    hold on;
    if ~exist('acom', 'var')
        acom = dlmread(acomname) ;
    end
    if ~exist('pcom', 'var')
        pcom = dlmread(pcomname) ;
    end
    if ~exist('dcom', 'var')
        dcom = dlmread(dcomname) ;
    end
    
    plot3(acom(1) * ssfactor, acom(2) * ssfactor, acom(3) * ssfactor, 'o')
    plot3(pcom(1) * ssfactor, pcom(2) * ssfactor, pcom(3) * ssfactor, 'o')
    plot3(dcom(1) * ssfactor, dcom(2) * ssfactor, dcom(3) * ssfactor, 'o')
    axis equal
    xlabel('x [pix]')
    ylabel('y [pix]')
    zlabel('z [pix]')
    if ii == 1
        view(0, 90)
    elseif ii == 2
        view(90, 0)
    else
        view(180, 0)
    end
end
sgtitle('APD COMs for APDV coordinates')
saveas(gcf, fullfile(QS.dir.mesh, 'apd_coms_APDVCoords.png'))
axis equal
%%%%%%%%%%%%%%%%%%%%%%
if preview
    waitfor(fig)
end
close all