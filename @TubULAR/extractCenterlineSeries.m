function extractCenterlineSeries(QS, cntrlineOptions)
% EXTRACTCENTERLINESERIES(QS, cntrlineOptions)
%   Extract centerlines for all timepoints in the QuapSlap instance
%
% Parameters
% ----------
% QS : 
%   uses in particular
%       timePoints : 
%       meshDir
%       flipy : bool (optional, default is true)
%           APDV coordinate system is mirrored across XZ wrt data coordinate
%           system XYZ. 
% Options: struct with fields
%   - overwrite : bool (optional, default is false)
%       whether to overwrite previous results on disk
%   - exponent : float (optional, default is 1.0)
%       exponent of the distance transform to use as speed
%   - res : float (optional, default is 1.0)
%       resolution with which to sample data volume to build shortest path
%   - xx, yy, zz : Nx1 float or int arrays (optional)
%       coordinates of the data volume in which to find the shortest path
%   - xyzlim : 3x1 or 3x2 or 1x3 float (required if xx,yy,zz not given)
%       extrema (maxima or maxima & minima) of each dimension for mesh
%       extent
%   - meshAPDVFileName : str (optional)
%       If supplied, plot rotated and scaled centerline with the previously
%       saved APDV aligned meshes matching this filename pattern
%   - reorient_faces : bool (optional, default is true)
%       ensure proper face orientation on each mesh (slower but rigorous)
%
% Outputs
% -------
% outDirs : 2 x 1 cell array 
%   The 2 output directories: centerlineOutDir, outdir for images
%
% NPMitchell 2020

timePoints = QS.xp.fileMeta.timePoints ;
meshFileName = QS.fullFileBase.mesh ;
startendptH5FileName = QS.fileName.startendPt ;
fn = QS.fileBase.name ;
plotStyleMesh = 'surface' ;  % [surface or fast] whether to plot the mesh carefully as surface or as (fast) scatterplot
saveImages = true ;

% Required options transform dataspace [pixels] into lab space [um]
resolution = QS.APDV.resolution ;
[rot, trans] = QS.getRotTrans() ;

% Unpack optional options
[xyzlim, ~, ~, xyzlim_um] = QS.getXYZLims() ;  % xyz maxima for grid formation if xx,yy,zz not supplied
overwrite = false ;             % overwrite previous results
overwrite_ims = false ;         % overwrite previous images, whether or not centerline is new
weight = 0.1;                   % for speedup of centerline extraction. Larger is less precise
exponent = 1.0 ;                % how heavily to scale distance transform for speed through voxel
res = 1.0 ;                     % resolution of distance tranform grid
flipy = QS.flipy ;              % APDV coordinate system is mirrored across XZ wrt data coords           
reorient_faces = true ;         % whether to ensure proper face orientation on each mesh (slower but rigorous)
preview = false ;               % view intermediate results
xwidth = 16 ;                   % width of figure in cm
ywidth = 10 ;                   % height of figure in cm
skipErrors = true ;             % skip timepoints that return errors
epsilon = eps ;                 % small value to give as weight of outside region
useSavedAPDVMeshes = false ;    % load APDV meshes instead of transforming the data space meshes on the fly
meshAPDVFileName = QS.fullFileBase.alignedMesh ; 
permuteDataAxisOrder = 'xyz' ;

if isfield(cntrlineOptions, 'timePoints')
    timePoints = cntrlineOptions.timePoints ;
end
if isfield(cntrlineOptions, 'overwrite')
    overwrite = cntrlineOptions.overwrite ;
end
if isfield(cntrlineOptions, 'overwrite_ims')
    overwrite_ims = cntrlineOptions.overwrite_ims ;
end
if isfield(cntrlineOptions, 'weight')
    exponent = cntrlineOptions.weight ;
end
if isfield(cntrlineOptions, 'exponent')
    exponent = cntrlineOptions.exponent ;
end
if isfield(cntrlineOptions, 'res')
    res = cntrlineOptions.res ;
end
if isfield(cntrlineOptions, 'skipErrors')
    skipErrors = cntrlineOptions.skipErrors ;
end
if isfield(cntrlineOptions, 'meshAPDVFileName')
    useSavedAPDVMeshes = true ;
    if strcmp(cntrlineOptions.meshAPDVFileName(end-3:end), '.ply')
        meshAPDVFileName = cntrlineOptions.meshAPDVFileName ;
    else
        meshAPDVFileName = [cntrlineOptions.meshAPDVFileName '.ply'] ;
    end
end
if isfield(cntrlineOptions, 'xyzlim')
    xyzlim = cntrlineOptions.xyzlim ;
end
if isfield(cntrlineOptions, 'reorient_faces')
    reorient_faces = cntrlineOptions.reorient_faces ;
end
if isfield(cntrlineOptions, 'preview')
    preview = cntrlineOptions.preview ;
end
if isfield(cntrlineOptions, 'xyzlim_um')
    xyzlim_um = cntrlineOptions.xyzlim_um ;
end
if isfield(cntrlineOptions, 'epsilon')
    epsilon = cntrlineOptions.epsilon ;
end
if isfield(cntrlineOptions, 'plotStyleMesh')
    plotStyleMesh = cntrlineOptions.plotStyleMesh ;
end
if isfield(cntrlineOptions, 'saveImages')
    saveImages = cntrlineOptions.saveImages ;
end
if isfield(cntrlineOptions, 'permuteDataAxisOrder')
    permuteDataAxisOrder = cntrlineOptions.permuteDataAxisOrder ;
end


% Figure options
colors = define_colors ;
blue = colors(1, :) ;
red = colors(2, :) ;
green = colors(5, :) ;

%% Get xyz grid for distance transform 
% use extrema of mesh vertices to clip the volume a bit for speed
if isfield(cntrlineOptions, 'xx') && isfield(cntrlineOptions, 'yy') && isfield(cntrlineOptions, 'zz')
    xx = cntrlineOptions.xx ;
    yy = cntrlineOptions.yy ;
    zz = cntrlineOptions.zz ;
else
    % by construction, the data volume must live in positive definite space
    % Check if limits passed by shape: have we been given [xmax,ymax,zmax] 
    % or [xmin xmax; ymin ymax; zmin zmax]?
    if all(size(xyzlim) > 1)
        xmax = xyzlim(1, 2) ;
        ymax = xyzlim(2, 2) ;
        zmax = xyzlim(3, 2) ;
    else
        xmax = xyzlim(1) ;
        ymax = xyzlim(2) ;
        zmax = xyzlim(3) ;
    end
    
    % Check that xyzlim are valid
    if ~all([xmax ymax zmax] > 0)
        disp('xyz limits are invalid, computing from mesh series')
        for tt = timePoints
            meshfn = sprintf(meshFileName, tt) ;
            mesh = read_ply_mod(meshfn) ;
            vv = mesh.v ;
            xmax = max(xmax, max(vv(:, 1))) ;
            ymax = max(ymax, max(vv(:, 2))) ;
            zmax = max(zmax, max(vv(:, 3))) ;
        end
        disp('done computing xmax ymax zmax from meshes')
    end
    xx = 0:res:ceil(xmax + res) ;
    yy = 0:res:ceil(ymax + res) ;
    zz = 0:res:ceil(zmax + res) ;
    
    % PERMUTE WRT DATA if desired/needed
    if ~strcmp(permuteDataAxisOrder, 'xyz')
        x1 = xx ;
        y1 = yy ; 
        z1 = zz ;
        if strcmp(permuteDataAxisOrder, 'xzy')
            yy = z1;
            zz = y1 ;
        elseif strcmp(permuteDataAxisOrder, 'yxz')
            xx = y1 ;
            yy = x1 ;
        elseif strcmp(permuteDataAxisOrder, 'yzx')
            xx = y1 ;
            yy = z1 ;
            zz = x1 ;
        elseif strcmp(permuteDataAxisOrder, 'zxy')
            xx = z1 ;
            yy = x1 ;
            zz = y1 ;
        elseif strcmp(permuteDataAxisOrder, 'zyx')
            xx = z1 ;
            zz = x1 ;
        else
            error('did not recognize desired permutation')
        end
    end
end


% Unpack output directories && ensure they exist
outdir = QS.dir.cntrline ;
figoutdir = fullfile(QS.dir.cntrline, 'images') ;
fig1outdir = fullfile(figoutdir, 'centerline_xy') ;
fig2outdir = fullfile(figoutdir, 'centerline_xz') ;
fig3outdir = fullfile(figoutdir, 'centerline_yz') ;
dirs2make = {outdir, figoutdir, fig1outdir, fig2outdir, fig3outdir} ;
for qq = 1:length(dirs2make)
    dir2make = dirs2make{qq} ;
    if ~exist(dir2make, 'dir') 
        mkdir(dir2make)
    end
end

% Unpack parameters into strings
expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
resstr = strrep(num2str(res, '%0.1f'), '.', 'p') ;
extenstr = ['_exp' expstr '_res' resstr] ;

for tt = timePoints
    % input filename mesh ending in 'ply'
    name = sprintf(QS.fileBase.mesh, tt) ;
    meshfn = sprintf(QS.fullFileBase.mesh, tt) ;
    disp(['Seeking centerline for mesh: ' meshfn])
    assert(exist(meshfn, 'file') > 0)
    
    %% Name the output centerline
    outname = [fullfile(outdir, name) '_centerline' extenstr] ;
    skel_rs_outfn = [fullfile(outdir, name) '_centerline_scaled' extenstr ] ;
    fig1outname = [fullfile(fig1outdir, name) '_centerline' extenstr '_xy.png'] ;
    fig2outname = [fullfile(fig2outdir, name) '_centerline' extenstr '_xz.png'] ;
    fig3outname = [fullfile(fig3outdir, name) '_centerline' extenstr '_yz.png'] ;
    fig1anyres_fn = [fullfile(fig1outdir, name) '_centerline' '_exp' expstr '*_xy.png'] ;
    fig2anyres_fn = [fullfile(fig2outdir, name) '_centerline' '_exp' expstr '*_xz.png'] ;
    fig3anyres_fn = [fullfile(fig3outdir, name) '_centerline' '_exp' expstr '*_yz.png'] ;
    
    % Check if any centerline exists at a different scale
    other_exist = ~isempty(dir([fullfile(outdir, name) '_centerline*'])) ;
    
    %% Compute centerline if has not been saved 
    if overwrite || (~exist(outname, 'file') && ~other_exist)
        if (~exist(outname, 'file') && ~other_exist)
            disp(['no centerline for timepoint = ' num2str(tt)])
        else
            disp(['Overwriting centerline for t=' num2str(tt)])
        end
        % Load startpoint and endpoint
        tifname = sprintf(fn, tt) ;
        disp(['loading /' tifname '/spt,ept,dpt from ' startendptH5FileName])
        startpt = h5read(startendptH5FileName, ['/' tifname '/spt' ]) ;
        endpt = h5read(startendptH5FileName, ['/' tifname '/ept' ]) ;
        
        % Read the mesh
        mesh = read_ply_mod(meshfn) ;
        % create faces/vertices struct for inpolyhedron
        fv.faces = mesh.f ;
        fv.vertices = mesh.v ;
        
        % Declare in command output
        if exist(outname, 'file')
            disp(['OVERWRITING: ' outname '.txt'])
        else
            disp(['Computing for ' outname '.txt'])
        end
        
        tic 
        if reorient_faces
            disp('Reorienting mesh faces... ') ;
            fv.faces = reorient_facets( fv.vertices, fv.faces );
        end
        
        % check it
        % trisurf(fv.faces, fv.vertices(:, 1), fv.vertices(:, 2), fv.vertices(:, 3))
        % hold on;
        % plot3(xx, 0*xx, 0*xx, 's') ;
        % plot3(0*yy, yy, 0*yy, 'o') ;
        % plot3(0*zz, 0*zz, zz, '^') ;
        % set(gcf, 'visible', 'on')
        % error('df')
        
        
        disp(['Identifying pts in mesh with inpolyhedron: ' name]) ;
        % note that inpolyhedron outputs YXZ format
        insideM = inpolyhedron(fv, xx, yy, zz) ;
        disp('> Computed segmentation:')
        toc ; 
        
        % Preview insideM
        if preview
            for page = 1:size(insideM, 3)
                imagesc(squeeze(insideM(:,:,page)))
                axis equal
                pause(0.01)
            end
        end

        % Optionally dilate the solid segmentation
        if cntrlineOptions.dilation > 0
            disp(['dilating inside volume by ' num2str(cntrlineOptions.dilation) ' voxels'])
            for qq = 1:cntrlineOptions.dilation
                [xb,yb,zb] = ndgrid(-3:3);
                se = strel(sqrt(xb.^2 + yb.^2 + zb.^2) <=3);
                insideM = imdilate(insideM, se) ;
            end
        end
        
        % Check that segmentation is a single solid
        props = regionprops3(insideM, 'Volume') ;
        ndilate = 0 ;
        while length(props.Volume) > 1
            disp('dilating inside volume by two voxels')
            [xb,yb,zb] = ndgrid(-4:4);
            se = strel(sqrt(xb.^2 + yb.^2 + zb.^2) <=4);
            insideM = imdilate(insideM, se) ;
            props = regionprops3(insideM, 'Volume') ;
            ndilate = ndilate + 1;
            if ndilate > 10
                error('could not dilate array to connect components')
            end
        end
        
        % use the distanceTransform from Yuriy Mishchenko
        disp(['Computing DT for ' name]) ;
        tic
        outside = 1 - insideM ;
        Dbw = bwdistsc(outside) ;
        % DD = max(DD(:)) - DD ;
        DD = (Dbw + epsilon) ./ (max(Dbw(:)) + epsilon) ;
        % DD = 1 - DD ;
        DD = DD.^(exponent) ; 
        DD(logical(outside)) = epsilon ;
        disp('> Computed DT:')
        toc ; 
        

        if preview
            % Preview DD
            close all ;
            disp('Previewing the distance transform')
            for ll=1:3
                for kk=1:10:size(DD,1)
                    imagesc(squeeze(DD(kk,:,:)))
                    title(['DT: z=' num2str(kk)])
                    colorbar
                    pause(0.001)
                end
            end

            % A better way to plot it
            clf
            p = patch(isosurface(xx,yy,zz,insideM,0.5));
            % isonormals(x,y,z,v,p)
            p.FaceColor = 'red';
            p.EdgeColor = 'none';
            daspect([1 1 1])
            view(3); 
            axis tight
            camlight 
            % lighting gouraud
        end

        % Check points with subsampling
        % ssample = 10 ;
        % xp = X(:); yp=Y(:); zp=Z(:) ; dp=D(:);
        % scatter3(xp(1:ssample:end), yp(1:ssample:end), ...
        %          zp(1:ssample:end), 30, dp(1:ssample:end), 'filled') ;

        %% use Peyre's fast marcher
        disp(['Computing fast marching for DT of ' name]);
        
        tic
        % From example (DD is W, with low values being avoided)
        options.heuristic = weight * DD ;
        % Convert here to the gridspacing of xx,yy,zz
        startpt_transposed = [startpt(2), startpt(1), startpt(3)]' / res ;
        endpt_transposed = [endpt(2), endpt(1), endpt(3)]' / res ;
        
        [D2,S] = perform_fast_marching(DD, startpt_transposed, options);
        
        try
            assert(length(find(D2(:) < Inf)) > 1)
            assert(any(D2(D2(:)< 1e9)))
        catch
            % A better way to plot it
            clf
            p = patch(isosurface(insideM,0.5));
            hold on;
            scatter3(startpt(1)/res, startpt(2)/res, startpt(3)/res, 30, 'filled')
            scatter3(endpt(1)/res, endpt(2)/res, endpt(3)/res, 30, 'filled')
            % isonormals(x,y,z,v,p)
            p.FaceColor = 'red';
            p.EdgeColor = 'none';
            p.FaceAlpha = 0.2 ;
            daspect([1 1 1])
            view(3); 
            axis tight
            camlight 
            pause(1)
            if ~skipErrors
                error('Could not build path')
            end
        end
        
        % Note: S tells us if the pixel distance has been computed (<0)
        path = compute_geodesic(D2, endpt_transposed);
        % plot_fast_marching_3d(D2, S, path, startpt, endpt);

        % Show the intermediate result
        disp('> Found skel via geodesic fast marching')        
        toc
        if preview
            % Preview DD
            clf ;
            for kk=1:10:size(DD,1)
                imagesc(squeeze(DD(kk,:,:)))
                title(['DT for plane z=' num2str(kk)])
                pause(0.01)
            end

            % Preview D2
            for kk=1:10:size(D2,1)
                imagesc(squeeze(D2(kk,:,:)))
                title(['D2 for plane z=' num2str(kk)])
                caxis([0,  max(max(D2(D2(:) < Inf)), 1)])
                pause(0.01)
            end
            
            % A better way to plot it
            clf
            for ii = 1:4
                subplot(4, 1, ii)
                p = patch(isosurface(insideM,0.5));
                hold on;
                scatter3(startpt(1)/res, startpt(2)/res, startpt(3)/res, 30, 'filled')
                scatter3(endpt(1)/res, endpt(2)/res, endpt(3)/res, 30, 'filled')
                % isonormals(x,y,z,v,p)
                p.FaceColor = 'red';
                p.EdgeColor = 'none';
                p.FaceAlpha = 0.2;
                daspect([1 1 1])
                view(3); 
                axis tight
                plot3(path(2,:), path(1, :), path(3,:), '-')
                trisurf(triangulation(fv.faces, fv.vertices/res), 'edgecolor', 'none', 'facealpha', 0.5)
                view([0, ii*45])
            end
            pause(1)
        end

        % Convert skeleton's rows to columns and flip start/end
        skel_tmp = fliplr(path)' ;
        disp(['Found a path of length ' num2str(length(skel_tmp))])
        if length(skel_tmp) > 3
            % Transpose x<->y back to original and scale to mesh units
            skel = [ skel_tmp(:,2), skel_tmp(:,1), skel_tmp(:,3) ] * res ;

            % Save centerline as text file
            fid = fopen([outname '.txt'], 'wt');
            % make header
            fprintf(fid, 'centerline xyz in units of pixels (full resolution, same as mesh)');  
            fclose(fid);
            disp(['Saving centerline to txt: ', outname, '.txt'])
            dlmwrite([outname '.txt'], skel)

            %% Rotate and scale (and mirror if flipy==true)
            % get distance increment
            ds = vecnorm(diff(skel), 2, 2) ;
            % get pathlength at each skeleton point
            ss = [0; cumsum(ds)] ;
            sss = ss * resolution ;
            
            skelrs = QS.xyz2APDV(skel) ;
            % The above line does the following transformation(s):
            % skelrs = ((rot * skel')' + trans) * resolution ;
            % if flipy
            %     skelrs(:, 2) = -skelrs(:, 2) ;
            % end

            % Check it
            QS.setTime(tt)
            amesh0 = QS.getCurrentAlignedMesh() ;
            amesh = struct() ;
            amesh.f = fv.faces ;
            amesh.v = QS.xyz2APDV(fv.vertices) ;
            subplot(1, 3, 1)
            trisurf(triangulation(amesh.f, amesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            plot3(skelrs(:, 1), skelrs(:, 2), skelrs(:, 3), '.-')
            axis equal ; view(0, 90)
            subplot(1, 3, 2)
            trisurf(triangulation(amesh.f, amesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            plot3(skelrs(:, 1), skelrs(:, 2), skelrs(:, 3), '.-')
            axis equal ; view(0, 0)      
            subplot(1, 3, 3)
            trisurf(triangulation(amesh.f, amesh.v), 'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            plot3(skelrs(:, 1), skelrs(:, 2), skelrs(:, 3), '.-')
            axis equal ; view(0, 90)
            
            
            %% Save the rotated, translated, scaled curve
            disp(['Saving rotated & scaled skeleton to txt: ', skel_rs_outfn, '.txt'])
            skeloutfn = [skel_rs_outfn '.txt'] ;
            fid = fopen(skeloutfn, 'wt');
            % make header
            fprintf(fid, 'Aligned & scaled skeleton: sss [um], skelrs [um]');  
            fclose(fid);
            dlmwrite(skeloutfn, [sss, skelrs])
            
            result_exists = true ;
            result_changed = true ;
            resolution_matches = true ;
        else
            if skipErrors
                disp('WARNING: PATH FAILED. SKIPPING THIS TIMEPOINT')
                % Preview D2
                clf ;
                for kk=[1:2:size(D2,1), size(D2,1):-1:1] 
                    imagesc(squeeze(D2(kk,:,:)))
                    title(['D2 for plane z=' num2str(kk)])
                    pause(0.001)
                end

                result_exists = false ;
                result_changed = false ;
                resolution_matches = true ;
            else
                error('PATH FAILED.')
            end
        end
    else     
        if exist(outname, 'file')
            disp('centerline exists on disk')
            result_exists = true ;
            result_changed = false ;
            resolution_matches = true ;
        else
            disp(['t=' num2str(tt) ': Centerline exists on disk ', ...
                'with different resolution, skipping tp'])
            result_exists = true ;
            result_changed = false ;
            resolution_matches = false ;
        end
    end
    
    % Check for existing figures of this timepoint centerline, but of ANY
    % resolution
    fig1s = ~isempty(dir(sprintf(fig1anyres_fn, tt))) ;
    fig2s = ~isempty(dir(sprintf(fig2anyres_fn, tt))) ;
    fig3s = ~isempty(dir(sprintf(fig3anyres_fn, tt))) ;
    fig_saved = fig1s && fig2s && fig3s ;
    if result_exists && (result_changed || ~fig_saved || overwrite_ims) && saveImages
        %% Plot and save
        if resolution_matches
            disp(['Loading skelrs from disk: ' skeloutfn])
            skeloutfn = [skel_rs_outfn '.txt'] ;
        else
            other_skelrs = dir([fullfile(outdir, name) '_centerline_scaled*']) ;
            skeloutfn = fullfile(other_skelrs(1).folder, other_skelrs(1).name) ;
        end
        ssskelrs = importdata(skeloutfn) ;
        skelrs = ssskelrs(:, 2:4) ;
        % Load start, end, and dorsal points in APDV coord system
        tifname = sprintf(fn, tt) ;
        sptrs = h5read(startendptH5FileName, ['/' tifname '/sptrs' ]) ;
        eptrs = h5read(startendptH5FileName, ['/' tifname '/eptrs' ]) ;
        dptrs = h5read(startendptH5FileName, ['/' tifname '/dptrs' ]) ;
        
        
        % Save plot of rotated and translated mesh
        disp('Saving rotated & translated figure (xy)...')    
        close all
        % Load rotated mesh
        % if useSavedAPDVMeshes
        disp(['Loading ' sprintf(meshAPDVFileName, tt)])
        mesh = read_ply_mod(sprintf(meshAPDVFileName, tt)) ;
        xyzrs = mesh.v ;
        % If we are plotting the mesh, reverse the triangle ordering
        % for ambient occlusion to work properly, but the APDV mesh
        % should ALREADY have the correct y values (mirrored XZ)
        if flipy
            tri = mesh.f(:, [1 3 2]) ;
        else
            tri = mesh.f ;
        end
        % else
        %     xyzrs = ((rot * mesh.v')' + trans) * resolution ;
        %     if flipy
        %         xyzrs(:, 2) = -xyzrs(:, 2) ;
        %         tri = mesh.f(:, [1 3 2]) ;
        %     end
        % end

        % Plot the result
        fig = figure('Visible', 'Off') ;
        if strcmpi( plotStyleMesh, 'fast')
            subidx = 1:round(length(xyzrs(:, 1))/2000):length(xyzrs(:, 1)) ;
            
            tmp = scatter3(xyzrs(subidx, 1), xyzrs(subidx,2), xyzrs(subidx, 3), ...
                10, xyzrs(subidx, 3), 'filled', 'markerfacealpha', 0.3) ;
        else
            tmp = trisurf(tri, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
                'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        end
        % [~,~,~] = apply_ambient_occlusion(tmp, 'SoftLighting', true) ; 
        hold on;
        % plot the skeleton
        plot3(skelrs(:,1), skelrs(:,2), skelrs(:,3), ...
            '-','Color',[0,0,0], 'LineWidth', 3);

        % annotate figure with APDV
        plot3(sptrs(1), sptrs(2), sptrs(3), 's', 'color', red)
        plot3(eptrs(1), eptrs(2), eptrs(3), '^', 'color', blue)
        plot3(dptrs(1), dptrs(2), dptrs(3), 'o', 'color', green)
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex');
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex');
        title(['Centerline using $D^{' num2str(exponent) '}$: ' num2str(tt)], ...
            'Interpreter', 'Latex')
        axis equal
        % xy
        view(2)
        xlim(xyzlim_um(1, :)); 
        ylim(xyzlim_um(2, :)); 
        zlim(xyzlim_um(3, :)) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, sprintf(fig1outname, tt))
        % yz
        disp('Saving rotated & translated figure (yz)...')    
        view(90, 0);
        xlim(xyzlim_um(1, :)); 
        ylim(xyzlim_um(2, :)); 
        zlim(xyzlim_um(3, :)) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, sprintf(fig2outname, tt))
        % xz
        disp('Saving rotated & translated figure (xz)...')  
        view(0, 0)    
        xlim(xyzlim_um(1, :)); 
        ylim(xyzlim_um(2, :)); 
        zlim(xyzlim_um(3, :)) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, sprintf(fig3outname, tt))
        close all
    end
    
end
