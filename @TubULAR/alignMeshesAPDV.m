function [rot, trans, xyzlim_raw, xyzlim, xyzlim_um, xyzlim_um_buff] = ...
    alignMeshesAPDV(QS, opts)
% ALIGNMESHESAPDV(opts) 
% Uses anterior, posterior, and dorsal training in ilastik h5 output to
% align meshes along APDV coordinate system with global rotation matrix 
% and translation vector. Extracted COMs from the 
% segmented training is loaded/saved in h5 file opts.rawapdvname (usually
% "apdv_coms_from_training.h5").
% Smoothed COMs from the segmented training --> opts.rawapdvname
% Smoothed rotated scaled COMs              --> opts.outapdvname
% Note that the global rotation matrix for the QuapSlap instance is defined
% previously in computeAPDCOMs()
% 
% 
% This is a function similar to the align_meshes_APDV.m script
%
% Parameters
% ----------
% opts : struct with fields
%   overwrite     : bool
%   overwrite_ims : bool
%   smwindow      : float or int
%       number of timepoints over which we smooth
%
% Returns
% -------
% xyzlim_raw : 
%   xyzlimits of raw meshes in units of full resolution pixels (ie not
%   downsampled)
% xyzlim : 
%   xyzlimits of rotated and translated meshes in units of full resolution 
%   pixels (ie not downsampled)
% xyzlim_um : 
%   xyz limits of rotated and translated meshes in microns
% xyzlim_um_buff : 
%   xyz limits of rotated and translated meshes in microns, with padding of
%   QS.normalShift * resolution in every dimension
%
% OUTPUTS
% -------
% xyzlim.txt 
%   xyzlimits of raw meshes in units of full resolution pixels (ie not
%   downsampled)
% xyzlim_APDV.txt 
%   xyzlimits of rotated and translated meshes in units of full resolution 
%   pixels (ie not downsampled)
% xyzlim_APDV_um.txt 
%   xyz limits of rotated and translated meshes in microns
% rotation_APDV.txt
%   rotation matrix to align mesh to APDV frame, saved to 
%   fullfile(meshDir, 'rotation_APDV.txt') ;
% translation_APDV.txt
%   translation vector to align mesh to APDV frame, saved to 
%   fullfile(meshDir, 'translation_APDV.txt') 
% xyzlim.txt 
%   raw bounding box in original frame (not rotated), in full res pixels.
%   Saved to fullfile(meshDir, 'xyzlim.txt')
% xyzlim_APDV.txt
%   bounding box in rotated frame, in full resolution pixels. Saved to 
%   fullfile(meshDir, 'xyzlim_APDV.txt')
% xyzlim_APDV_um.txt
%   bounding box in rotated frame, in microns. Saved to 
%   fullfile(meshDir, 'xyzlim_APDV_um.txt')
% apdv_coms_rs.h5 (outapdvname)
%   Centers of mass for A, P, and D in microns in rotated, scaled APDV
%   coord system. Note that this coord system is mirrored if flipy==true.
%   Also contains raw acom,pcom,dcom in subsampled pixels.
%   Saved to fullfile(meshDir, 'centerline/apdv_coms_rs.h5')
% startendpt.h5
%   Starting and ending points
%   Saved to fullfile(meshDir, 'centerline/startendpt_rs.h5') ;
% 
% NPMitchell 2020

% Unpack QS
timePoints = QS.xp.fileMeta.timePoints ;
meshDir = QS.dir.mesh ;
resolution = QS.APDV.resolution ;
[acom_sm, pcom_sm] = QS.getAPCOMSm() ;

% Default options
overwrite = false ;
overwrite_ims = false ;
preview = false ;
plot_buffer = 20 ;
ssfactor = QS.ssfactor ;
flipy = QS.flipy ; 

% Booleans & floats
if isfield(opts, 'overwrite')
    overwrite = opts.overwrite ;  % overwrite everything 
end
if isfield(opts, 'overwrite_ims')
    overwrite_ims = opts.overwrite_ims ;  % overwrite images, whether or not we overwrite everything else
end
if isfield(opts, 'preview')
    preview = opts.preview ;
end
if isfield(opts, 'plot_buffer')
    plot_buffer = opts.plot_buffer ;
end
normal_step = opts.normal_step ;
dcomname = fullfile(meshDir, 'dcom_for_rot.txt') ;

% Default valued options
timeinterval = 1 ;
timeunits = 'min' ;
if isfield(opts, 'timeinterval')
    timeinterval = opts.timeinterval ;
end
if isfield(opts, 'timeunits')
    timeunits = opts.timeunits ;
end

% Data file names
rotname = QS.fileName.rot ;
transname = QS.fileName.trans ;
xyzlimname_raw = QS.fileName.xyzlim_raw ;
xyzlimname_pix = QS.fileName.xyzlim_pix ;
xyzlimname_um = QS.fileName.xyzlim_um ;
xyzlimname_um_buff = QS.fileName.xyzlim_um_buff ;
% Name output directory for apdv info
apdvoutdir = QS.dir.cntrline ;
outapdvname = fullfile(apdvoutdir, 'apdv_coms_rs.h5') ;
outstartendptname = fullfile(apdvoutdir, 'startendpt.h5') ;
% Name the directory for outputting aligned_meshes
alignedMeshDir = QS.dir.alignedMesh ;
meshFileName = QS.fullFileBase.mesh ;
alignedMeshBase = QS.fullFileBase.alignedMesh ;
alignedMeshXYFigBaseName = [QS.fileBase.alignedMesh '_xy.png'] ;
alignedMeshXZFigBaseName = [QS.fileBase.alignedMesh '_xz.png'] ;
alignedMeshYZFigBaseName = [QS.fileBase.alignedMesh '_yz.png'] ;
% Note: used to define fn = QS.fileBase.name ;

% rotname
if isfield(opts, 'rotname')
    rotname = opts.rotname ;
end
if ~strcmp(rotname(end-3:end), '.txt') 
    rotname = [rotname '.txt'] ;
end

% transname
if isfield(opts, 'transname')
    transname = opts.transname ;
end
if ~strcmp(transname(end-3:end), '.txt') 
    transname = [transname '.txt'] ;
end

if isfield(opts, 'outapdvname')
    outapdvname = opts.outapdvname ;
end

% dcomname
if isfield(opts, 'dcomname')
    dcomname = opts.dcomname ;
end
if ~strcmp(dcomname(end-3:end), '.txt') 
    dcomname = [dcomname '.txt'] ;
end
if isfield(opts, 'rawapdvname')
    apdvoutdir = opts.apdvOutDir ;
end

% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm
colors = define_colors ;
blue = colors(1, :) ;
red = colors(2, :) ;
green = colors(5, :) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks, directories, and assertions
% Ensure that PLY files exist
for tidx = 1:length(timePoints)
    tt = timePoints(tidx) ;
    if ~exist(sprintf(meshFileName, tt), 'file')
        msg = ['Found no matching PLY files ', ...
                sprintf(meshFileName, tt), ' in ', meshDir ]; 
        error(msg)
    end
end

% Name the directory for outputting figures
figoutdir = fullfile(alignedMeshDir, 'images');
fig1outdir = fullfile(figoutdir, 'aligned_mesh_xy') ;
fig2outdir = fullfile(figoutdir, 'aligned_mesh_xz') ;
fig3outdir = fullfile(figoutdir, 'aligned_mesh_yz') ;

% Create the directories 
dirs2make = {apdvoutdir, alignedMeshDir, figoutdir, ...
    fig1outdir, fig2outdir, fig3outdir} ;
for kk = 1:length(dirs2make)
    thisdir = dirs2make{kk} ;
    if ~exist(thisdir, 'dir')
        mkdir(thisdir) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get axis limits from looking at all meshes =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
if exist(xyzlimname_raw, 'file')
    disp('loading xyzlimits from disk')
    xyzlim_raw = dlmread(xyzlimname_raw, ',', 1, 0);
    xmin = xyzlim_raw(1);
    ymin = xyzlim_raw(2);
    zmin = xyzlim_raw(3);
    xmax = xyzlim_raw(4);
    ymax = xyzlim_raw(5);
    zmax = xyzlim_raw(6);
else
    disp('Extracting xyzlimits for raw meshes...')
    disp([' ... since ' xyzlimname_raw ' does not exist'])
    for tidx = 1:length(timePoints)
        tt = timePoints(tidx) ;
        disp(['tt = ', num2str(tt)])
        % Get the timestamp string from the name of the mesh
        mesh = read_ply_mod(sprintf(meshFileName, tt)) ;

        minx = min(mesh.v) ;
        maxx = max(mesh.v) ;
        if tt == timePoints(1)
            xmin = minx(1) ;
            ymin = minx(2) ;
            zmin = minx(3) ;
            xmax = maxx(1) ;
            ymax = maxx(2) ;
            zmax = maxx(3) ;
        else
            xmin= min(xmin, minx(1)) ;
            ymin = min(ymin, minx(2)) ;
            zmin = min(zmin, minx(3)) ;
            xmax = max(xmax, maxx(1)) ;
            ymax = max(ymax, maxx(2)) ;
            zmax = max(zmax, maxx(3)) ;
        end 
    end

    % Save xyzlimits 
    % disp('Saving raw mesh xyzlimits for plotting')
    header = 'xyzlimits for original meshes in units of full resolution pixels' ; 
    write_txt_with_header(xyzlimname_raw, [xmin, xmax; ymin, ymax; zmin, zmax], header) ;
    % Now read it back in
    xyzlim_raw = dlmread(xyzlimname_raw, ',', 1, 0) ;
    % xyzlim_raw = [xmin, xmax; ymin, ymax; zmin, zmax] ;
end
disp('done')

%% With acoms and pcoms in hand, we compute dorsal and rot/trans ==========
overwrite_startendpts = false ;
for tidx = 1:length(timePoints)
    tic
    tt = timePoints(tidx) ;
    
    % Pick out the acom and pcom in SUBSAMPLED UNITS from smoothed sequence
    % NOTE: this is from the RAW data
    acom = acom_sm(tidx, :) ;
    pcom = pcom_sm(tidx, :) ; 
    
    %% Name the output centerline
    fig1outname = fullfile(fig1outdir, sprintf(alignedMeshXYFigBaseName, tt)) ;
    fig2outname = fullfile(fig2outdir, sprintf(alignedMeshXZFigBaseName, tt)) ;
    fig3outname = fullfile(fig3outdir, sprintf(alignedMeshYZFigBaseName, tt)) ; 
        
    %% Read the mesh  
    meshfn = sprintf(meshFileName, tt) ;
    disp(['Loading mesh ' meshfn])
    mesh = read_ply_mod(meshfn );
    vtx_sub = mesh.v / ssfactor ;
    vn = mesh.vn ;
    fvsub = struct('faces', mesh.f, 'vertices', vtx_sub, 'normals', vn) ;
    
    %% Does the output aligned mesh exist?
    alignedmeshfn = sprintf(alignedMeshBase, tt) ;
    meshfn_exist = exist(alignedmeshfn, 'file') ;    
    
    % Check normals
    % close all
    % plot3(vtx_sub(1:10:end, 1), vtx_sub(1:10:end, 2), vtx_sub(1:10:end, 3), '.')
    % hold on
    % plot3(vtx_sub(1:10:end, 1) + 10 * vn(1:10:end, 1),...
    %     vtx_sub(1:10:end, 2) + 10 * vn(1:10:end, 2), ...
    %     vtx_sub(1:10:end, 3) + 10 * vn(1:10:end, 3), 'o')
    
    % View the normals a different way
    % close all
    % plot3(vtx_sub(1:10:end, 1), vtx_sub(1:10:end, 2), vtx_sub(1:10:end, 3), '.')
    % for i=1:10:length(vtx_sub)
    %     hold on
    %     plot3([vtx_sub(i, 1), vtx_sub(i, 1) + 10*vn(i, 1)], ... 
    %     [vtx_sub(i, 2), vtx_sub(i, 2) + 10*vn(i, 2)], ...
    %     [vtx_sub(i, 3), vtx_sub(i, 3) + 10*vn(i, 3)], 'r-') 
    % end
    % axis equal
    
    % Must either downsample mesh, compute xyzgrid using ssfactor and
    % pass to options struct.
    % Here, downsampled mesh
    % mesh.vertex.x = xs ;
    % mesh.vertex.y = ys ;
    % mesh.vertex.z = zs ;
    
    %% Try loading the spt and ept, or recompute
    % If aligned mesh doesn't exist, redo point-matching
    if meshfn_exist
        try
            name = sprintf(QS.fileBase.name, tt) ;
            spt = h5read(outstartendptname, ['/' name '/spt']) ;
            ept = h5read(outstartendptname, ['/' name '/ept']) ;
            if any(spt) && any(ept)
                spt_ept_exist = true ;
            end
        catch
            spt_ept_exist = false;
        end
    else
        % Redo point-matching since aligned mesh isn't saved
        spt_ept_exist = false ;
    end
    
    if overwrite || ~spt_ept_exist 
        % Point match for aind and pind
        disp(['Point matching mesh ' meshfn])
        adist2 = sum((vtx_sub - acom) .^ 2, 2);
        %find the smallest distance and use that as an index 
        aind = find(adist2 == min(adist2)) ;
        % Next point match the posterior
        pdist2 = sum((vtx_sub - pcom) .^ 2, 2);
        % find the smallest distance and use that as an index
        pind = find(pdist2 == min(pdist2)) ;

        % Check it
        if preview
            disp('Previewing mesh in figure window')
            trimesh(mesh.f, vtx_sub(:, 1), vtx_sub(:, 2), vtx_sub(:, 3), vtx_sub(:, 1))
            hold on;
            plot3(vtx_sub(aind, 1), vtx_sub(aind, 2), vtx_sub(aind, 3), 'ko')
            plot3(vtx_sub(pind, 1), vtx_sub(pind, 2), vtx_sub(pind, 3), 'ro')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Define start point and endpoint 
        disp(['Defining start point and endpoint for tp = ' num2str(tt)])
        % Check if acom is inside mesh. If so, use that as starting point.
        ainside = inpolyhedron(fvsub, acom(1), acom(2), acom(3)) ;
        pinside = inpolyhedron(fvsub, pcom(1), pcom(2), pcom(3)) ;

        if ainside
            disp('start point for centerline is inside mesh')
            startpt = acom' ;
        else
            % move along the inward normal of the mesh from the matched vertex
            vtx = [vtx_sub(aind, 1), vtx_sub(aind, 2), vtx_sub(aind, 3)]' ;
            normal = fvsub.normals(aind, :) ;
            startpt = vtx + normal;
            if ~inpolyhedron(fvsub, startpt(1), startpt(2), startpt(3)) 
                % this didn't work, check point in reverse direction
                startpt = vtx - normal * normal_step ;
                if ~inpolyhedron(fvsub, startpt(1), startpt(2), startpt(3))
                    % Can't seem to jitter into the mesh, so use vertex
                    disp("Can't seem to jitter into the mesh, so using vertex for startpt")
                    startpt = vtx ;
                end
            end
        end 
        % Note: Keep startpt in subsampled units

        % Define end point
        if pinside
            disp('end point for centerline is inside mesh')
            endpt = pcom' ;
        else
            % move along the inward normal of the mesh from the matched vertex
            vtx = [vtx_sub(pind, 1), vtx_sub(pind, 2), vtx_sub(pind, 3)]' ;
            normal = fvsub.normals(pind, :) ;
            endpt = vtx + normal * normal_step;
            if ~inpolyhedron(fvsub, endpt(1), endpt(2), endpt(3)) 
                % this didn't work, check point in reverse direction
                endpt = vtx - normal * normal_step ;
                if ~inpolyhedron(fvsub, endpt(1), endpt(2), endpt(3))
                    % Can't seem to jitter into the mesh, so use vertex
                    disp("Can't seem to jitter into the mesh, so using vertex for endpt")
                    endpt = vtx ;
                end
            end
        end 
        % Note: Keep endpt in subsampled units

        % Check out the mesh
        if preview
            hold on
            trimesh(fvsub.faces, xs, ys, zs)
            % plot3(xs, ys, zs, 'ko')
            scatter3(startpt(1), startpt(2), startpt(3), 'ro')
            scatter3(endpt(1), endpt(2), endpt(3), 'ko')
            xlabel('x [subsampled pixels]')
            ylabel('y [subsampled pixels]')
            zlabel('z [subsampled pixels]')
            hold off
            axis equal
        end

        %% Rescale start point and end point to full resolution
        spt = [startpt(1), startpt(2), startpt(3)] * ssfactor;
        ept = [endpt(1), endpt(2), endpt(3)] * ssfactor;
        overwrite_startendpts = true ;
        
        clearvars startpt endpt ainside pinside normal vtx
    else
        disp('loading spt and ept')
    end
    
    %% Grab rot, trans
    disp('Loading rot from disk...')
    rot = dlmread(rotname) ;
    dcom = dlmread(dcomname) ;
    trans = dlmread(transname, ',');
    
    %% Rotate and translate (and mirror) acom, pcom, dcom
    try 
        apdcoms_rs_exist = true ;
        name = sprintf(QS.fileBase.name, tt) ;
        acom_rs = h5read(outapdvname, ['/' name '/acom_rs']) ;
        pcom_rs = h5read(outapdvname, ['/' name '/pcom_rs']) ;
        dcom_rs = h5read(outapdvname, ['/' name '/dcom_rs']) ;
        
        % Check that matches what is stored
        acom_rs_new = ((rot * (acom' * ssfactor))' + trans) * resolution ;
        pcom_rs_new = ((rot * (pcom' * ssfactor))' + trans) * resolution ;
        dcom_rs_new = ((rot * (dcom' * ssfactor))' + trans) * resolution ;
        if flipy
            acom_rs_new = [ acom_rs_new(1) -acom_rs_new(2) acom_rs_new(3) ] ;
            pcom_rs_new = [ pcom_rs_new(1) -pcom_rs_new(2) pcom_rs_new(3) ] ;
            dcom_rs_new = [ dcom_rs_new(1) -dcom_rs_new(2) dcom_rs_new(3) ] ;
        end
        assert(all(abs(acom_rs_new - acom_rs) < 1e-6))
        assert(all(abs(pcom_rs_new - pcom_rs) < 1e-6))
        assert(all(abs(dcom_rs_new - dcom_rs) < 1e-6))
    catch
        disp('Rotated & Scaled APD COMS do not exist or are different on disk')
        apdcoms_rs_exist = false ;
    end
    
    if overwrite || ~apdcoms_rs_exist
        acom_rs = ((rot * (acom' * ssfactor))' + trans) * resolution ;
        pcom_rs = ((rot * (pcom' * ssfactor))' + trans) * resolution ;
        dcom_rs = ((rot * (dcom' * ssfactor))' + trans) * resolution ;

        if flipy
            acom_rs = [ acom_rs(1) -acom_rs(2) acom_rs(3) ] ;
            pcom_rs = [ pcom_rs(1) -pcom_rs(2) pcom_rs(3) ] ;
            dcom_rs = [ dcom_rs(1) -dcom_rs(2) dcom_rs(3) ] ;
        end
    end
    
    %% Rotate and translate vertices and endpoints
    % NOTE: ars = ((rot * a')' + trans) * QS.APDV.resolution ;
    xyzrs = ((rot * (vtx_sub * ssfactor)')' + trans) * resolution;
    % xyzrs2 = (rot * (mesh.v)' + trans')' * resolution ;
    vn_rs = (rot * fvsub.normals')' ;
    sptr = (rot * spt')' + trans ; 
    eptr = (rot * ept')' + trans ;
    dpt = dcom' * ssfactor ;
    dptr = (rot * (dcom' * ssfactor))' + trans ; 
    
    % Scale to actual resolution
    sptrs = sptr * resolution ;
    eptrs = eptr * resolution ; 
    dptrs = dptr * resolution ;
    
    % Flip in Y if data is reflected across XZ
    if flipy
        % Note: since normals point inward along y when y is flipped, it
        % remains only to flip normals along X and Z in the second line.
        xyzrs = [xyzrs(:, 1), -xyzrs(:, 2), xyzrs(:, 3)] ;  % flip vertices
        % xyzrs2 = [xyzrs2(:, 1), -xyzrs2(:, 2), xyzrs2(:, 3)] ;  % flip vertices
        vn_rs = [-vn_rs(:, 1), vn_rs(:, 2), -vn_rs(:, 3)] ; % flip normals > normals point inward
        sptrs = [sptrs(1), -sptrs(2), sptrs(3)] ;           % flip startpt
        eptrs = [eptrs(1), -eptrs(2), eptrs(3)] ;           % flip endpt
        dptrs = [dptrs(1), -dptrs(2), dptrs(3)] ;           % flip dorsalpt
    else
        vn_rs = -vn_rs ;    % flip normals > normals point inward 
        % mesh.f = mesh.f(:, [1, 3, 2]) ;
    end
    
    %% Update our estimate for the true xyzlims
    
    %% Get a guess for the axis limits if this is first TP
    if tidx == 1 
        % Check if already saved, and load or recompute
        % fntmp = xyzlimname_um ;
        [~, ~, ~, xyzlim_um_buff] = QS.getXYZLims() ;
        % Expand xyzlimits for plots
        xminrs_plot = xyzlim_um_buff(1,1) - plot_buffer ;
        yminrs_plot = xyzlim_um_buff(2,1) - plot_buffer ;
        zminrs_plot = xyzlim_um_buff(3,1) - plot_buffer ;
        xmaxrs_plot = xyzlim_um_buff(1,2) + plot_buffer ;
        ymaxrs_plot = xyzlim_um_buff(2,2) + plot_buffer ;
        zmaxrs_plot = xyzlim_um_buff(3,2) + plot_buffer ;    
    end
    
    %% Check the rotation
    if tidx == 1
        close all
        fig = figure('Visible', 'off') ;
        tmp = trisurf(mesh.f, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
                    xyzrs(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.5) ;
        try
            [~,~,~] = apply_ambient_occlusion(tmp, 'SoftLighting', true) ; % 'ColorMap', viridis) ;
        catch
            disp('Could not apply ambient occlusion! GPToolbox likely not configured')
        end
        hold on;
        xyz = vtx_sub;
        
        % Aligned meshes have inward pointing normals, so flip them for
        % plotting ambient occlusion (irrespective of flipy, I believe)
        faces_to_plot = mesh.f(:, [2, 1, 3]) ;
        
        tmp2 = trisurf(faces_to_plot, xyz(:, 1), xyz(:,2), xyz(:, 3), ...
            xyz(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.5) ;
        clearvars faces_to_plot
        try
            [~,~,~] = apply_ambient_occlusion(tmp2, 'SoftLighting', true) ; % 'ColorMap', viridis) ;
        catch
            debugMsg(1, 'Could not apply ambient occlusion \n')
        end
        boxx = [xmin, xmin, xmin, xmin, xmax, xmax, xmax, xmax, xmin] ;
        boxy = [ymin, ymax, ymax, ymin, ymin, ymax, ymax, ymin, ymin] ;
        boxz = [zmin, zmin, zmax, zmax, zmax, zmax, zmin, zmin, zmin] ;
        box = [boxx', boxy', boxz'] ;
        box_sub = box / ssfactor ; 
        boxrs = resolution * ((rot * box')' + trans) ;
        plot3(box_sub(:, 1), box_sub(:, 2), box_sub(:, 3), 'k-')
        plot3(boxrs(:, 1), boxrs(:, 2), boxrs(:, 3), 'k-')
        for i=1:3
            plot3([boxrs(i, 1), box_sub(i, 1)], ...
                [boxrs(i, 2), box_sub(i, 2)], ...
                [boxrs(i, 3), box_sub(i, 3)], '--')
        end
          
        % plot the skeleton
        % for i=1:length(skelrs)
        %     plot3(skelrs(:,1), skelrs(:,2), skelrs(:,3),'-','Color',[0,0,0], 'LineWidth', 3);
        % end
        plot3(sptrs(1), sptrs(2), sptrs(3), 'ro')
        plot3(eptrs(1), eptrs(2), eptrs(3), 'bo')
        plot3(dptrs(1), dptrs(2), dptrs(3), 'go')
        plot3(acom(1), acom(2), acom(3), 'rx')
        plot3(pcom(1), pcom(2), pcom(3), 'bx')
        plot3(dcom(1), dcom(2), dcom(3), 'gx')

        xlabel('x [$\mu$m or pix]', 'Interpreter', 'Latex'); 
        ylabel('y [$\mu$m or pix]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m or pix]', 'Interpreter', 'Latex');
        title('Checking rotation')
        axis equal
        saveas(fig, fullfile(alignedMeshDir, 'rot_check.png'))
        close all
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save the rotated, translated, scaled to microns mesh ===============
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if overwrite || ~meshfn_exist    
        disp('Saving the aligned mesh...')
        disp([' --> ' alignedmeshfn])
        plywrite_with_normals(alignedmeshfn, mesh.f, xyzrs, vn_rs)
    else
        disp(['alignedMesh PLY exists on disk (' alignedmeshfn ')'])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot and save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Save plot of rotated and translated mesh
    figs_do_not_exist = ~exist(fig1outname, 'file') || ...
        ~exist(fig2outname, 'file') || ~exist(fig3outname, 'file');
    
    if overwrite_ims || figs_do_not_exist || ~meshfn_exist
        disp('Saving rotated & translated figure (xy)...')    
        close all
        fig = figure('Visible', 'Off') ;
        if flipy
            faces_to_plot = mesh.f(:, [2, 1, 3]) ;
        else
            faces_to_plot = mesh.f ;
        end
        
        th = trisurf(faces_to_plot, xyzrs(:, 1), xyzrs(:, 2), xyzrs(:, 3), ...
            'edgecolor', 'none', 'facecolor', 'w', 'FaceAlpha', 0.5) ;
        % 'FaceVertexCData',bsxfun(@times,(1-AO),C)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check xyzrs against xyz2APDV(xyz)
        % fig2=figure(2)
        % QS.setTime(tt)
        % rmesh = QS.getCurrentRawMesh() ;
        % vtx2 = QS.xyz2APDV(rmesh.v) ;
        % hold on;
        % th = trisurf(faces_to_plot, vtx2(:, 1), vtx2(:, 2), vtx2(:, 3), ...
        %     'edgecolor', 'none', 'facecolor', 'y', 'FaceAlpha', 0.5) ;
        % th = trisurf(faces_to_plot, xyzrs(:, 1), xyzrs(:, 2), xyzrs(:, 3), ...
        %     'edgecolor', 'none', 'facecolor', 'k', 'FaceAlpha', 0.1) ;
        % % th = trisurf(faces_to_plot, xyzrs2(:, 1), xyzrs2(:, 2), xyzrs2(:, 3), ...
        % %     'edgecolor', 'none', 'facecolor', 'g', 'FaceAlpha', 0.5) ;
        axis equal
        close(fig2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        try
            [~,~,~] = apply_ambient_occlusion(th, 'SoftLighting', true) ;
        catch
            disp('Could not apply ambient occlusion')
        end
        
        % Figure properties
        axis equal
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
        set(gcf, 'color', 'k')
        titlestr = ['Aligned mesh, $t=$' num2str(tt * timeinterval) ' ' timeunits] ;
        title(titlestr, 'interpreter', 'latex', 'color', 'w'); 
        grid off
        
        % check it 
        % set(gcf, 'visible', 'on')
        % waitfor(gcf)
        % error('here')

        hold on;
        plot3(sptrs(1), sptrs(2), sptrs(3), 'o', 'color', red)
        plot3(eptrs(1), eptrs(2), eptrs(3), 'o', 'color', blue)
        plot3(dptrs(1), dptrs(2), dptrs(3), 'o', 'color', green)
        plot3(acom_rs(1), acom_rs(2), acom_rs(3), 's', 'color', red)
        plot3(pcom_rs(1), pcom_rs(2), pcom_rs(3), '^', 'color', blue)
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex'); 
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex');
        
        % xy
        view(2)
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        disp(['Saving to ' fig1outname])
        % saveas(fig, fig1outname)
        export_fig(fig1outname, '-nocrop', '-r150')
        
        % yz
        disp('Saving rotated & translated figure (yz)...')    
        view(90, 0);
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);  % x_width=10cm y_width=15cm
        % saveas(fig, fig2outname)
        export_fig(fig2outname, '-nocrop', '-r150')
        
        % xz
        disp('Saving rotated & translated figure (xz)...')  
        view(0, 0)    
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);  % x_width=10cm y_width=15cm
        % saveas(fig, fig3outname)
        export_fig(fig3outname, '-nocrop', '-r150')
        close all
    end
    
    %% Preview and save coms
    % Check the normals 
    if preview 
        close all
        plot3(vtx_rs(1:10:end, 1), vtx_rs(1:10:end, 2), vtx_rs(1:10:end, 3), '.')
        for i=1:10:length(vtx_rs)
            hold on
            plot3([vtx_rs(i, 1), vtx_rs(i, 1) + 10*vn_rs(i, 1)], ... 
            [vtx_rs(i, 2), vtx_rs(i, 2) + 10*vn_rs(i, 2)], ...
            [vtx_rs(i, 3), vtx_rs(i, 3) + 10*vn_rs(i, 3)], 'r-') 
        end
        axis equal
    end    
    
    % Save acom, pcom and their aligned counterparts as attributes in an
    % hdf5 file            
    name = sprintf(QS.fileBase.name, tt) ;
    % Save if overwrite
    if overwrite || ~apdcoms_rs_exist
        try
            h5create(outapdvname, ['/' name '/acom'], size(acom)) ;
        catch
            disp('acom already exists as h5 file. Overwriting.')
        end
        try
            h5create(outapdvname, ['/' name '/pcom'], size(pcom)) ;
        catch
            disp('pcom already exists as h5 file. Overwriting.')
        end
        try 
            h5create(outapdvname, ['/' name '/dcom'], size(dcom)) ;
        catch
            disp('dcom already exists as h5 file. Overwriting.')
        end
        try
            h5create(outapdvname, ['/' name '/acom_rs'], size(acom_rs)) ;
        catch
            disp('acom_rs already exists as h5 file. Overwriting.')
        end
        try
            h5create(outapdvname, ['/' name '/pcom_rs'], size(pcom_rs)) ;
        catch
            disp('pcom_rs already exists as h5 file. Overwriting.')
        end
        try 
            h5create(outapdvname, ['/' name '/dcom_rs'], size(dcom_rs)) ;
        catch
            disp('dcom_rs already exists as h5 file. Overwriting.')
        end
        h5write(outapdvname, ['/' name '/acom'], acom) ;
        h5write(outapdvname, ['/' name '/pcom'], pcom) ;
        h5write(outapdvname, ['/' name '/dcom'], dcom) ;
        h5write(outapdvname, ['/' name '/acom_rs'], acom_rs) ;
        h5write(outapdvname, ['/' name '/pcom_rs'], pcom_rs) ;
        h5write(outapdvname, ['/' name '/dcom_rs'], dcom_rs) ;
        % h5disp(outapdvname, ['/' name]);
        disp('Saved h5: scom pcom dcom acom_rs pcom_rs dcom_rs')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save the startpt/endpt, both original and rescaled to um ===========
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save acom, pcom and their aligned counterparts as attributes in an
    % hdf5 file -- these are mesh-dependent since point-matched to vertices
    name = sprintf(QS.fileBase.name, tt) ;
    
    if overwrite || ~spt_ept_exist || overwrite_startendpts
        disp('Saving the startpt/endpt...')
        try
            h5create(outstartendptname, ['/' name '/spt'], size(spt)) ;
        catch
            disp('spt already exists as h5 file. Overwriting.')
        end
        try
            h5create(outstartendptname, ['/' name '/ept'], size(ept)) ;
        catch
            disp('ept already exists as h5 file. Overwriting.')
        end
        try 
            h5create(outstartendptname, ['/' name '/dpt'], size(dpt)) ;
        catch
            disp('dpt already exists as h5 file. Overwriting.')
        end
        try
            h5create(outstartendptname, ['/' name '/sptrs'], size(eptrs)) ;
        catch
            disp('sptrs already exists as h5 file. Overwriting.')
        end
        try
            h5create(outstartendptname, ['/' name '/eptrs'], size(eptrs)) ;
        catch
            disp('eptrs already exists as h5 file. Overwriting.')
        end
        try 
            h5create(outstartendptname, ['/' name '/dptrs'], size(dptrs)) ;
        catch
            disp('dptrs already exists as h5 file. Overwriting.')
        end

        h5write(outstartendptname, ['/' name '/spt'], spt) ;
        h5write(outstartendptname, ['/' name '/ept'], ept) ;
        h5write(outstartendptname, ['/' name '/dpt'], dpt) ;
        h5write(outstartendptname, ['/' name '/sptrs'], sptrs) ;
        h5write(outstartendptname, ['/' name '/eptrs'], eptrs) ;
        h5write(outstartendptname, ['/' name '/dptrs'], dptrs) ;
        disp('Saved h5: spt ept dpt sptrs eptrs dptrs')
    else
        disp('startpt/endpt already exist and not overwriting...')
    end
    toc
end

% % Save xyzlim_raw
% if overwrite || ~exist(xyzlimname_raw, 'file')
%     disp('Saving rot/trans mesh xyzlimits for plotting')
%     header = 'xyzlimits for raw meshes in units of full resolution pixels' ;
%     xyzlim = [xmin, xmax; ymin, ymax; zmin, zmax] ;
%     write_txt_with_header(xyzlimname_raw, xyzlim, header) ;
% else
%     xyzlim_raw = [xmin, xmax; ymin, ymax; zmin, zmax] ;
% end
% 
% % Save xyzlimits 
% if overwrite || ~exist(xyzlimname_pix, 'file')
%     disp('Saving rot/trans mesh xyzlimits for plotting')
%     header = 'xyzlimits for rotated translated meshes in units of full resolution pixels' ;
%     xyzlim = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] / resolution;
%     write_txt_with_header(xyzlimname_pix, xyzlim, header) ;
% else
%     xyzlim = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] / resolution;
% end
% 
% % Save xyzlimits in um
% if overwrite || ~exist(xyzlimname_um, 'file')
%     disp('Saving rot/trans mesh xyzlimits for plotting, in microns')
%     header = 'xyzlimits for rotated translated meshes in microns' ;
%     xyzlim_um = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] ;
%     write_txt_with_header(xyzlimname_um, xyzlim_um, header) ;
% end
% 
% % Save buffered xyzlimits in um
% if overwrite || ~exist(xyzlimname_um_buff, 'file')
%     disp('Saving rot/trans mesh xyzlimits for plotting, in microns')
%     header = 'xyzlimits for rotated translated meshes in microns, with padding (buffered)' ;
%     xyzlim_um = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] ;
%     xyzlim_um_buff = xyzlim_um + 2 * abs(QS.normalShift) * resolution * [-1, 1] ;
%     write_txt_with_header(xyzlimname_um_buff, xyzlim_um_buff, header) ;
% end

disp('done')