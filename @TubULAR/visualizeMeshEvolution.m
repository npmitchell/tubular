function visualizeMeshEvolution(QS, options)
% Plot meshes as a morphsnakes demo
% Load all meshes in demo output directory, pass them through meshlab
% filters. Plot them in orthogonal sectioning view in 3D.
% Given plys from morphsnakes in ./growing_plys/, deposit smoothed meshes
% into ./growing_meshlab_plys/, then make figures in ./demofigs
% This method is based on figure_demo_morphsnakes_visualize.m script
%
% NPMitchell 2021

%% Unpack options
if nargin < 2
    options = struct() ;
end

% Default options 
overwrite = false ;         % bool, overwrite output files on disk
plot_growth = true ;       % bool, plot the growth of the mesh to fill the organ for a single (first?) timepoint
plot_evolution = true ;     % bool, plot the evolution of the mesh over time
faceon = false ;            % visualize the results laterally rather than perspective
preview = false ;           % preview intermediate results
plotXslice = false ;        % optionally create orthoview of Xslice
brighten = 255 ;            % factor by which to brighten scaled IV data
x0 = 0 ;                    % x value for orthogonal slice
y0 = 0 ;                   % y value for orthogonal slice
z0 = -5 ;                   % z value for orthogonal slice
viewAngles = [-0.75, -1.0, 0.7] ;  % viewing angles for perspective
ambientStrength = 0.5 ;     % intrinsic (isotropic) brightess for orthosections
ambientStrength_meshOrtho = 1.0 ; 
forcetrue = true ;
growtht0 = QS.xp.fileMeta.timePoints(1) ;
lwX = 2 ;

% texturepatch options
meshFileBase = QS.fullFileBase.mesh ;
normal_shift = QS.normalShift ;
flipy = QS.flipy ;
texture_axis_order = QS.data.axisOrder ;
t0 = QS.t0set() ;

% figure parameters
xwidth = 32 ; % cm
ywidth = 20 ; % cm

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'plot_growth')
    plot_growth = options.plot_growth ;
end
if isfield(options, 'plot_evolution')
    plot_evolution = options.plot_evolution ;
end
if isfield(options, 'growth_t0')
    growtht0 = options.growth_t0 ;
end
if isfield(options, 'faceon')
    faceon = options.faceon ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'plotXslice')
    plotXslice = options.plotXslice ;
end
if isfield(options, 'adjust_high')
    adjust_high = options.adjust_high ;
end
if isfield(options, 'x0')
    x0 = options.x0 ;
end
if isfield(options, 'y0')
    y0 = options.y0 ;
end
if isfield(options, 'z0')
    z0 = options.z0 ;
end
if isfield(options, 't0')
    t0 = options.t0 ;
end
if isfield(options, 'xwidth')
    xwidth = options.xwidth ;
end
if isfield(options, 'ywidth')
    ywidth = options.ywidth ;
end

%% Colors
colors = define_colors(7) ; 
blue   = colors(1, :) ;
red    = colors(2, :) ;
yellow = colors(3, :) ;
purple = colors(4, :) ;
green  = colors(5, :) ;
sky    = colors(6, :) ;
maroon = colors(7, :) ;
lscolor = sky ;

%% Directories
if plot_growth
    outdirGPLY = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', ...
        sprintf('growth_plys_%06dt0', growtht0)) ;
    outdirGF = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', ...
        sprintf('growth_faceon_%06dt0', growtht0)) ;
    outdirGP = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', ...
        sprintf('growth_perspective_%06dt0', growtht0)) ;
    dirs2do = {outdirGPLY, ...
        fullfile(outdirGP, 'mesh'), ...
        fullfile(outdirGF, 'mesh'), ...
        fullfile(outdirGP, 'mesh_only'), ...
        fullfile(outdirGF, 'mesh_only'), ...
        fullfile(outdirGP, 'texture'), ...
        fullfile(outdirGF, 'texture')} ;
    for qq =  1:length(dirs2do)
        dir2make = dirs2do{qq} ;
        if ~exist(dir2make, 'dir')
            mkdir(dir2make)
        end
    end
end
% for plotting evolution
% ----------------------
outdirF = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', 'faceon') ;
outdirP = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', 'perspective') ;
dirs2do = {fullfile(outdirP, 'mesh_only'), ...
    fullfile(outdirF, 'mesh_only')} ;
for qq =  1:length(dirs2do)
    dir2make = dirs2do{qq} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PART I: GROWING SEQUENCE ===============================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_growth
    %% First run something like:
    %     
    % # Construct mesh for each timepoint
    % datDir="/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/demo_morphsnakes_figs/growing_meshlab_plys/";
    % 
    % cd $datDir
    % tp=85;
    % for (( num=0; num<=50; num++ )); do        
    % 	
    % 	tpx=$(printf "%03d" $(( tp )));     
    % 	idx=$(printf "%03d" $(( num )));     
    % 	prev=$(printf "%03d" $(( num-1 )))     
    %         
    % 	mslsDir="${datDir}msls_output_growEvolveDemo_tp${tpx}/";
    % 
    % 	# for all iterations, use the selected timepoint's h5 output
    %   initls=${mslsDir}msls_grow000${prev}.h5
    %             
    % 	python /mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/run_morphsnakes.py -i Time_000${tpx}_c1_stab_Probabilities.h5 -init_ls $initls -o  $mslsDir -prenu 0 -presmooth 0 -ofn_ply mesh_grow000${idx}.ply -ofn_ls msls_grow000${idx}.h5 -l1 1 -l2 1 -nu 0.0 -postnu 3 -smooth 0.2 -postsmooth 5 -exit 0.00010 -channel 1 -dtype h5 -permute zyxc -ss 4 -include_boundary_faces -center_guess 175,75,100 -rad0 10 -n 4 ; 
    % done

    
    %% Options
    use_pointcloud = false ; % Use the pointcloud from the level set rather
                             % than the boundary mesh from marching cubes.
                             % Requires the level sets as h5 files
    dtype = '.h5' ;
    mlxprogram = QS.xp.detectOptions(1).mlxprogram ;
    
    %% Directories
    plyDir = fullfile(outdirGPLY, 'growing_plys') ;
    meshlabOutDir = fullfile(outdirGPLY, 'growing_meshlab_plys') ;
    plyfn = fullfile(plyDir, 'mesh_grow%06d.ply') ;  
    plySmFn = fullfile(meshlabOutDir, 'mesh_growSmooth_%06d.ply') ;
    ofn_ply = 'mesh_grow' ;
    ofn_smoothply = 'mesh_growSmooth_' ;
    plyglob = fullfile(plyDir, 'mesh_grow*.ply') ;
    
    
    %% find all ms_...ply files in mslsDir, and smooth them all
    files_to_smooth = dir(plyglob) ;
    lsfns_to_smooth = dir(fullfile(plyDir, ['msls_grow*' dtype])) ;

    idx2do = [1:length(files_to_smooth)]
    for i=idx2do 
        msls_mesh_outfn = files_to_smooth(i).name ;
        PCfile = fullfile( files_to_smooth(i).folder, msls_mesh_outfn );
        % Note that LS file is outputLs ;
        split_fn = strsplit(msls_mesh_outfn, ofn_ply) ;
        extension_outfn = split_fn{end} ;
        mesh_outfn = [ofn_smoothply, extension_outfn] ;
        outputMesh = fullfile(meshlabOutDir, mesh_outfn);

        disp(['outputMesh = ', outputMesh])
        %bad = so_bad
        if ~exist( outputMesh, 'file')
            if use_pointcloud
                % Use the pointcloud from the level set rather than the
                % boundary mesh from marching cubes
                %----------------------------------------------------------------------
                % Extract the implicit level set as a 3D binary array
                %----------------------------------------------------------------------

                % The file name of the current time point
                ls_outfn_ii = lsfns_to_smooth(i).name ;
                % The 3D binay array
                bwLS = h5read( ls_outfn_ii, '/implicit_levelset' );

                % Extract the (x,y,z)-locations of the level set boundary (in pixel
                % space)
                bwBdyIDx = bwperim( bwLS );

                clear bwBdy
                [ bwBdy(:,1), bwBdy(:,2), bwBdy(:,3) ] = ind2sub( size(bwLS), ...
                    find(bwBdyIDx) );

                %----------------------------------------------------------------------
                % Create output mesh
                %----------------------------------------------------------------------

                % Write the points to a .obj file as a point cloud for ouput to Meshlab
                clear OBJ
                OBJ.vertices = bwBdy;
                OBJ.objects(1).type='f';
                OBJ.objects(1).data.vertices=[];

                pointcloud_fn = [base_outfn_for_pointcloud '.obj'] ;
                disp(['Writing point cloud ' pointcloud_fn]);
                write_wobj(OBJ, pointcloud_fn );

                % Run the meshlab script
                system( ['meshlabserver -i ' pointCloudFileName, ...
                    ' -o ' outputMesh, ...
                    ' -s ' mlxprogram ' -om vn']);
            else
                % Use the marching cubes mesh surface to smooth
                command = ['meshlabserver -i ' PCfile ' -o ' outputMesh, ...
                    ' -s ' mlxprogram ' -om vn'];
                % Either copy the command to the clipboard
                clipboard('copy', command);
                % or else run it on the system
                disp(['running ' command])
                system(command)
            end
        else
            disp(['t=', num2str(i) ': smoothed mesh file found...'])
        end
    end
        
    %% APDV frame
    [rot, trans] = QS.getRotTrans ;
    resolution = QS.APDV.resolution ;
    
    first = true ;
    tp = growtht0 ;    
    QS.setTime(tp) ;
    nFiles = length(dir(plyglob)) ;
    
    for tt = -1:(nFiles-1)
        % Make init file
        if tt >= 0
            timestr = sprintf('%06d', tt) ;
            meshfn = sprintf(plySmFn, tt) ;
            disp(['Considering time ' timestr])

            fnGP = fullfile(outdirGP, 'mesh', [sprintf(QS.fileBase.name, tt) '.png']) ;
            fnGF = fullfile(outdirGF, 'mesh', [sprintf(QS.fileBase.name, tt) '.png']) ;
            fn_textureGP = fullfile(outdirGP, 'texture', [sprintf(QS.fileBase.name, tt) '.png']) ;
            fn_textureGF = fullfile(outdirGF, 'texture', [sprintf(QS.fileBase.name, tt) '.png']) ;
            use_mesh = true ;
        else                
            fnGP = fullfile(outdirGP, 'mesh', 'data_only.png') ;
            fnGF = fullfile(outdirGF, 'mesh', 'data_only.png') ;
            fn_textureGP = fullfile(outdirGP, 'texture', 'data_only.png') ;
            fn_textureGF = fullfile(outdirGF, 'texture', 'data_only.png') ;
            use_mesh = false ;
        end

        if ~exist(fnGP, 'file') || ~exist(fnGF, 'file') || overwrite 
            % Plot with h5
            ploth5 = false ;
            if ploth5
                rawh5fn = fullfile(projectDir, sprintf(h5fn, tp)) ;
                % Load the raw data from h5
                raw = h5read(rawh5fn, '/inputData') ;
            else
                % tiffn_i = fullfile(projectDir, sprintf(tiffn, str2double(timestr))) ;
                % Load the raw data from tiff
                % raw = readTiff4D(tiffn_i, 1, 'xyz') ;
                raw = QS.getCurrentData() ;
                IV = raw{1} ;
            end

            if use_mesh
                % Raw and APDV aligned Meshes
                gmesh = read_ply_mod(meshfn) ;
                amesh = gmesh ;
                amesh.v = QS.xyz2APDV(gmesh.v) ;

                % scaled mesh vertices
                xa = amesh.v(:, 1)  ;
                ya = amesh.v(:, 2)  ;
                za = amesh.v(:, 3)  ;
            end

            %% Plot with orthogonal sectioning (note: faceon gives face-on view(2))
            close all
            fig = figure('units', 'centimeters', ...
                'outerposition', [0 0 xwidth xwidth], 'visible', 'off') ;

            %% Interpolate the data and evaluate at rotated orthogonal planes     
            % PASSIVE ROTATION            
            % Obtain rotation matrix that undoes the APDV rotation 
            invRot = QS.invertRotation(rot) ;

            % Plane to image
            xMin = -50 ;
            xMax = 300 ;
            yMin = -90 ;
            yMax = 90 ;
            zMin = -90 ;
            zMax = 80 ;
            oversample_factor = 1.1 ;
            nX = round((xMax - xMin) / resolution * oversample_factor) ;
            nY = round((yMax - yMin) / resolution * oversample_factor) ; 
            nZ = round((zMax - zMin) / resolution * oversample_factor) ;
            x1 = x0 * ones(nZ, nY) ;
            y2 = -y0 * ones(nZ, nX) ;
            z3 = z0 * ones(nY, nX) ;
            xspace = linspace(xMin, xMax, nX) ;
            yspace = linspace(yMin, yMax, nY) ;
            zspace = linspace(zMin, zMax, nZ) ;
            [y1, z1] = meshgrid(yspace, zspace) ; 
            [x2, z2] = meshgrid(xspace, zspace) ; 
            [x3, y3] = meshgrid(xspace, yspace) ; 
            xyzX = (invRot * (([x1(:), y1(:), z1(:)]/resolution) - trans)')' ;
            xyzY = (invRot * (([x2(:), y2(:), z2(:)]/resolution) - trans)')' ;
            xyzZ = (invRot * (([x3(:), y3(:), z3(:)]/resolution) - trans)')' ;
            xpx = reshape(xyzX(:, 1), [nY, nZ]) ;
            xpy = reshape(xyzX(:, 2), [nY, nZ]) ; 
            xpz = reshape(xyzX(:, 3), [nY, nZ]) ;
            ypx = reshape(xyzY(:, 1), [nX, nZ]) ;
            ypy = reshape(xyzY(:, 2), [nX, nZ]) ; 
            ypz = reshape(xyzY(:, 3), [nX, nZ]) ;
            zpx = reshape(xyzZ(:, 1), [nX, nY]) ;
            zpy = reshape(xyzZ(:, 2), [nX, nY]) ; 
            zpz = reshape(xyzZ(:, 3), [nX, nY]) ;
            assert(all(size(x1) == size(y1)))
            assert(all(size(y1) == size(z1)))
            assert(all(size(x2) == size(y2)))
            assert(all(size(y2) == size(z2)))
            assert(all(size(x3) == size(y3)))
            assert(all(size(y3) == size(z3)))

            % preview with mesh
            disp('performing interpolation in Xplane')
            ix = interp3(double(IV), xyzX(:, 2), xyzX(:, 1), xyzX(:, 3)) ;
            ix = reshape(ix, size(x1)) ;

            disp('performing interpolation in Yplane')
            iy = interp3(double(IV), xyzY(:, 2), xyzY(:, 1), xyzY(:, 3)) ;
            iy = reshape(iy, size(x2)) ;

            disp('performing interpolation in Zplane')
            iz = interp3(double(IV), xyzZ(:, 2), xyzZ(:, 1), xyzZ(:, 3)) ;
            iz = reshape(iz, size(x3)) ;

            if preview
                % Check order of axes wrt IV
                figure ;
                subplot(2, 2, 1)
                imagesc(1:size(IV, 2), 1:size(IV, 3), squeeze(IV(round(mean(xyzX(:, 1))), :, :))')
                xlabel('size(IV, 2)')
                ylabel('size(IV, 3)'); axis equal; title('IX (approx, pre-rot)')
                subplot(2, 2, 2)
                imagesc(1:size(IV, 1), 1:size(IV, 3), squeeze(IV(:, round(mean(xyzY(:, 2))), :))')
                xlabel('size(IV, 1)')
                ylabel('size(IV, 3)'); axis equal; title('IY (approx, pre-rot)')
                subplot(2, 2, 3)
                imagesc(1:size(IV, 1), 1:size(IV, 2), squeeze(IV(:, :, round(mean(xyzZ(:, 3)))))')
                xlabel('size(IV, 1)')
                ylabel('size(IV, 2)'); axis equal; title('IZ (approx, pre-rot)')           

                % check interpolation
                figure; 
                subplot(2, 2, 1); imagesc(yspace, zspace, ix'); axis equal ;
                subplot(2, 2, 2); imagesc(xspace, zspace, iy'); axis equal ;
                subplot(2, 2, 3); imagesc(xspace, yspace, iz'); axis equal ;
                waitfor(gcf)
                close all
            end

            ix(isnan(ix)) = 0 ;
            iy(isnan(iy)) = 0 ;
            iz(isnan(iz)) = 0 ;
            ix = brighten * ix / 2^(16) ;
            iy = brighten * iy / 2^(16) ;
            iz = brighten * iz / 2^(16) ;

            %% Plot in APDV
            clf            
            hold on;
            surface(x2, -y2, z2, iy, 'FaceColor','texturemap', ...
                'EdgeColor','none','CDataMapping','direct', ...
                                'AmbientStrength', ambientStrength_meshOrtho)
            hold on;

            % Plot X slice (optional)
            if ~faceon && plotXslice                
                surface(x1, -y1, z1, ix, 'FaceColor','texturemap', ...
                    'EdgeColor','none','CDataMapping','direct',...
                                'AmbientStrength', ambientStrength_meshOrtho)
            end

            % Plot Z slice (if not looking laterally)
            if ~faceon
                surface(x3, -y3, z3,  iz, 'FaceColor','texturemap', ...
                   'EdgeColor','none','CDataMapping','direct',...
                                'AmbientStrength', ambientStrength_meshOrtho)
                hold on;
            end
            colormap bone

            if use_mesh
                h = trimesh(amesh.f, xa, ya, za, 'EdgeColor', 'none', ...
                    'facecolor', lscolor) ;
            end
            axis equal
            hold on;
            axis off

            lighting gouraud    % preferred method for lighting curved surfaces
            material dull    % set material to be dull, no specular highlights

            % Lighting and view 
            % view(-20, 20)
            view(viewAngles) 
            lgt = camlight('headlight') ;

            % Obtain xyzlimits
            if first
                xlims = xlim ;
                ylims = ylim ;
                zlims = zlim ;
                first = false ;
            end
            xlim(xlims) ;
            ylim(ylims) ;
            zlim(zlims) ;
            disp(['Saving ' fnGP])

            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

            % saveas(fig, fn) ;
            export_fig(fnGP, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fnP);

            % Face-On view
            view([0,0]) 
            disp(['saving face-on view: ' fnGF])
            export_fig(fnGF, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fnF);
            close all
        end
        
        %% Create matching Texturepatch Image
        if ~exist(fn_textureGP, 'file') ||  ~exist(fn_textureGF, 'file') || overwrite
             % Plot with h5 or data
            ploth5 = false ;
            if ploth5
                rawh5fn = fullfile(projectDir, sprintf(h5fn, tp)) ;
                % Load the raw data from h5
                raw = h5read(rawh5fn, '/inputData') ;
            else
                % tiffn_i = fullfile(projectDir, sprintf(tiffn, str2double(timestr))) ;
                % Load the raw data from tiff
                % raw = readTiff4D(tiffn_i, 1, 'xyz') ;
                raw = QS.getCurrentData() ;
                IV = raw{1} ;
            end

            if use_mesh
                % Raw growing mesh
                gmesh = read_ply_mod(meshfn) ;
                amesh = gmesh ;
                amesh.v = QS.xyz2APDV(gmesh.v) ;
                
                % scaled mesh vertices
                xr = gmesh.v(:, 1) ;
                yr = gmesh.v(:, 2) ;
                zr = gmesh.v(:, 3) ;
            end
            
            % Psize is the linear dimension of the grid drawn on each triangular face
            Options = struct() ;
            Options.PSize = 8;
            Options.EdgeColor = 'none';
            %Options.Rotation = rot ;
            %Options.Translation = trans ;
            %Options.Dilation = resolution ;
            % OUTWARD, INWARD layers
            Options.numLayers = [5, 2];  
            % Note: [2, 2] marches ~0.524 um in either dir if layer spacing is 1.0
            Options.layerSpacing = 1 ;
            
            fig = figure('Visible', 'Off') ;
            disp(['creating texture patch ' num2str(tp, '%06d')])

            % Allow for axis swapping
            if use_mesh
                TV = gmesh.v(:, texture_axis_order) ;
                texture_patch_3d( gmesh.f, amesh.v, gmesh.f, TV, ...
                     QS.getCurrentData(), Options );
                colormap bone
            end
             
            % format the figure
            disp('formatting figure...')
            axis equal
            [~,~,~,xyzlim] = QS.getXYZLims() ;
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            view(viewAngles)
            
            titlestr = ['$t = $' num2str(tp*QS.timeInterval-t0) ' ' QS.timeUnits] ;
            title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
            % xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex')
            % ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
            % zlabel('DV position [$\mu$m]', 'Interpreter', 'Latex')
            axis off
            
            % Set size of figure to be same as orthoview figure
            fig = gcf ;
            set(fig, 'Units', 'centimeters');
            set(fig, 'Position', [0 0 xwidth ywidth]);

            export_fig(fn_textureGP, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_textureP);
            
            % Face-On view
            view([0,0]) 
            disp(['saving face-on view: ' fn_textureGF])
            export_fig(fn_textureGF, '-nocrop', '-transparent')            
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_textureF);
            
            %% Add orthoviews to texturepatch image
            % hold on;
            % surface(x2, -y2, z2, iy, 'FaceColor','texturemap', ...
            %     'EdgeColor','none','CDataMapping','direct',...
            %     'AmbientStrength', ambientStrength)
            % hold on;
            % 
            % % Plot X slice (optional)
            % if ~faceon && plotXslice                
            %     surface(x1, -y1, z1, ix, 'FaceColor','texturemap', ...
            %         'EdgeColor','none','CDataMapping','direct',...
            %         'AmbientStrength', ambientStrength)
            % end
            % 
            % % Plot Z slice (if not looking laterally)
            % if ~faceon
            %     surface(x3, -y3, z3,  iz, 'FaceColor','texturemap', ...
            %        'EdgeColor','none','CDataMapping','direct',...
            %        'AmbientStrength', ambientStrength)
            %     hold on;
            % end
            % 
            % % Set size of figure to be same as orthoview figure
            % fig = gcf ; 
            % set(fig, 'Units', 'centimeters');
            % set(fig, 'Position', [0 0 xwidth ywidth]);
            % 
            % view(viewAngles)
            % export_fig(fn_texture2P, '-nocrop', '-transparent')
            % % patchIm = getframe(gca);
            % % imwrite( patchIm.cdata, fn_texture2P);
            % 
            % % Face-On view
            % view([0,0]) 
            % disp(['saving face-on view: ' fn_texture2F])
            % export_fig(fn_texture2F, '-nocrop', '-transparent')
            % % patchIm = getframe(gca);
            % % imwrite( patchIm.cdata, fn_texture2F);
            % close all
        end
    end 
end
    
    
    
%     
% 
%     %% Plot with orthogonal sectioning (note: faceon gives face-on view(2))
%     % brighten the raw data
%     rawclip = raw;
%     thres = 5000 ;
%     clip = 2 * median(raw(raw > thres)) ;
%     rawclip(raw > clip) = clip ;
% 
%     close all
%     fig = figure('units', 'centimeters', ...
%         'position', [0 0 xwidth ywidth], 'visible', 'off') ;
%     
%     % X=const slice
%     xraw = double(squeeze(rawclip(xslice, :, :))) ;
%     xraw = 64 * xraw / max(xraw(:)) ;
%     [xpy, xpz] = meshgrid(1:size(xraw, 1), 1:size(xraw, 2)) ;
%     xpx = xslice * ones(size(xpy)) ;
%     
%     % Y=const slice
%     if ~faceon && plotYslice
%         yraw = double(squeeze(rawclip(:, yslice, :))) ;
%         yraw = 64 * yraw / max(yraw(:)) ;
%         [ypx, ypz] = meshgrid(1:size(yraw, 1), 1:size(yraw, 2)) ;
%         ypy = yslice * ones(size(ypx)) ;
%     end
%     
%     % Z=const slice
%     if ~faceon
%         zraw = double(squeeze(rawclip(:, :, zslice))) ;
%         zraw = 64 * zraw / max(zraw(:)) ;
%         [zpx, zpy] = meshgrid(1:size(zraw, 1), 1:size(zraw, 2)) ;
%         zpz = zslice * ones(size(zpx)) ;
%     end
%     
%     % Plot X slice
%     surface(xpz, xpy, xpx, xraw', 'FaceColor','texturemap', ...
%         'EdgeColor','none','CDataMapping','direct')
%     hold on;
% 
%     % Plot Y slice
%     if ~faceon && plotYslice
%         surface(ypz, ypy, ypx, yraw', 'FaceColor','texturemap', ...
%             'EdgeColor','none','CDataMapping','direct')
%     end
%     
%     % Plot Z slice
%     if ~faceon
%         surface(zpz, zpy, zpx, zraw', 'FaceColor','texturemap', ...
%            'EdgeColor','none','CDataMapping','direct')
%     end
%     % axis equal  % debug: toggle this off?
%     colormap bone
%     % hold on;
% 
%     % Plot each in turn
%     fns = dir(fullfile(meshlabOutDir, [ofn_smoothply '*.ply'])) ;
%     for ii = 1:length(fns)
%         disp(['Reading ply: ' fns(ii).name])
%         mesh = read_ply_mod(fullfile(fns(ii).folder, fns(ii).name)) ;
%         name = strsplit(fns(ii).name, '.ply') ;
%         name = name{1} ;
%         xx = mesh.v(:, 1);
%         yy = mesh.v(:, 2);
%         zz = mesh.v(:, 3);
%         h = trimesh(mesh.f, xx, yy, zz, 'EdgeColor', 'none', 'facecolor', lscolor) ;
%         axis equal
%         axis off
% 
%         lighting gouraud    % preferred method for lighting curved surfaces
%         material dull    % set material to be dull, no specular highlights
% 
%         % Obtain xyzlimits
%         if ii  == 1
%             if faceon
%                 lgt = camlight ;
%             else
%                 view([2,2,2])
%                 lgt = camlight('headlight') ;
%             end
%             xlims = xlim ;
%             ylims = ylim ;
%             zlims = zlim ;
%         end
%         xlim(xlims) ;
%         ylim(ylims) ;
%         if ~faceon
%             zlim(zlims) ;
%         end
% 
%         % xlabel('x')
%         % ylabel('y')
%         % zlabel('z')
%         fnGP = fullfile(outdirGP, [name '.png']) ;
%         disp(['Saving ' fnGP])
%         export_fig(fnGP, '-nocrop', '-transparent')
% 
%         delete(h)
%     end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART IIa: Mesh only ========================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% navigate to the outputdir for figures
if plot_evolution
    projectDir = QS.dir.data ;
    meshDir = QS.dir.mesh ;
    [rot, trans] = QS.getRotTrans ;
    resolution = QS.APDV.resolution ;
    plyfn = QS.fullFileBase.mesh ;
    
    timePoints = QS.xp.fileMeta.timePoints ;
    
    first = true ;
    tidx2do = [123-30] ;
    tidx2do = [tidx2do, setdiff(1:30:length(timePoints), tidx2do)] ;
    tidx2do = [tidx2do, setdiff(1:10:length(timePoints), tidx2do)] ;
    tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
    for ii = tidx2do 
        tp = timePoints(ii) ;
        timestr = sprintf('%06d', tp) ;
        disp(['Considering time ' timestr])
        
        fn_meshP = fullfile(outdirP, 'mesh_only', [sprintf(QS.fileBase.name, tp) '.png']) ;
        fn_meshF = fullfile(outdirF, 'mesh_only', [sprintf(QS.fileBase.name, tp) '.png']) ;
        if ~exist(fn_meshP, 'file') || overwrite
            
            QS.setTime(tp) ;

            % Raw and APDV aligned Meshes
            amesh = QS.loadCurrentAlignedMesh() ;

            % scaled mesh vertices
            amesh.v = laplacian_smooth(amesh.v, amesh.f, 'cotan', [], 0.001,'implicit',amesh.v,1000) ;
            xa = amesh.v(:, 1)  ;
            ya = amesh.v(:, 2)  ;
            za = amesh.v(:, 3)  ;
            
            close all            
            fig = figure('units', 'centimeters', ...
                'outerposition', [0 0 xwidth xwidth], 'visible', 'off') ;

            h = trimesh(amesh.f, xa, ya, za, 'EdgeColor', 'none', ...
                'facecolor', lscolor) ;
            axis equal
            hold on;
            axis off

            lighting gouraud    % preferred method for lighting curved surfaces
            material dull    % set material to be dull, no specular highlights

            % Lighting and view 
            % view(-20, 20)
            view(viewAngles) 
            lgt = camlight('headlight') ;
            set(gcf,'color','w');

            % Obtain xyzlimits
            if first
                [~,~,~,xyzlim] = QS.getXYZLims() ;
                xlims = xyzlim(1, :) ;
                ylims = xyzlim(2, :) ;
                zlims = xyzlim(3, :) ;
                first = false ;
            end
            xlim(xlims) ;
            ylim(ylims) ;
            zlim(zlims) ;
            disp(['Saving ' fn_meshP])

            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

            % saveas(fig, fn) ;
            export_fig(fn_meshP, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_meshP);
            
            % Face-On view
            view([0,0]) 
            disp(['saving face-on view: ' fn_meshF])
            export_fig(fn_meshF, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_meshF);
            close all
        end    
    end
end
disp('done')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART II: EVOLVING SEQUENCE =============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% navigate to the outputdir for figures
if plot_evolution
    outdirF = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', 'faceon') ;
    outdirP = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', 'perspective') ;
    outdirX = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', 'outline') ;
    dirs2do = {fullfile(outdirP, 'mesh_only'), ...
        fullfile(outdirP, 'mesh'), ...
        fullfile(outdirP, 'texture'), ...
        fullfile(outdirP, 'texture_orthoviews'), ...
        fullfile(outdirF, 'mesh_only'), ...
        fullfile(outdirF, 'mesh'), ...
        fullfile(outdirX), ...
        fullfile(outdirF, 'texture'), ...
        fullfile(outdirF, 'texture_orthoviews')} ;
    for qq =  1:length(dirs2do)
        dir2make = dirs2do{qq} ;
        if ~exist(dir2make, 'dir')
            mkdir(dir2make)
        end
    end
    projectDir = QS.dir.data ;
    meshDir = QS.dir.mesh ;
    [rot, trans] = QS.getRotTrans ;
    resolution = QS.APDV.resolution ;
        
    % rotfn = fullfile(meshDir, 'rotation_APDV.txt') ;
    % transfn = fullfile(meshDir, 'translation_APDV.txt') ;
    % rot = dlmread(rotfn) ;
    % trans = dlmread(transfn) ;
    % resolution = 0.2619 ;
    
    plyfn = QS.fullFileBase.mesh ;
    % h5fn = 'Time_%06d_c1_stab.h5' ;
    % tiffn = 'Time_%06d_c1_stab.tif' ;
    % fns = dir(fullfile(meshDir, [plyfn '*.ply']));
    
    timePoints = QS.xp.fileMeta.timePoints ;
    
    first = true ;
    tidx2do = [QS.xp.tIdx(QS.t0set())-30:30:QS.xp.tIdx(length(timePoints))] ;
    tidx2do = fliplr(tidx2do) ;
    tidx2do = [tidx2do, setdiff(1:30:length(timePoints), tidx2do)] ;
    tidx2do = [tidx2do, setdiff(1:10:length(timePoints), tidx2do)] ;
    tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
    for ii = tidx2do
    
        tp = timePoints(ii) ;
        QS.setTime(tp) ;
        timestr = sprintf('%06d', tp) ;
        disp(['Considering time ' timestr])
        
        fnP = fullfile(outdirP, 'mesh', [sprintf(QS.fileBase.name, tp) '.png']) ;
        fnF = fullfile(outdirF, 'mesh', [sprintf(QS.fileBase.name, tp) '.png']) ;
        fnX = fullfile(outdirX, [sprintf(QS.fileBase.name, tp) '.png']) ;
        fn_textureP = fullfile(outdirP, 'texture', [sprintf(QS.fileBase.name, tp) '.png']) ;
        fn_textureF = fullfile(outdirF, 'texture', [sprintf(QS.fileBase.name, tp) '.png']) ;
        fn_texture2P = fullfile(outdirP, 'texture_orthoviews', [sprintf(QS.fileBase.name, tp) '.png']) ;
        fn_texture2F = fullfile(outdirF, 'texture_orthoviews', [sprintf(QS.fileBase.name, tp) '.png']) ;
        
        if ~exist(fnP, 'file') || ~exist(fn_textureP, 'file') || ...
                ~exist(fnX, 'file') || overwrite 
            % Plot with h5
            ploth5 = false ;
            if ploth5
                rawh5fn = fullfile(projectDir, sprintf(h5fn, str2num(timestr))) ;
                % Load the raw data from h5
                raw = h5read(rawh5fn, '/inputData') ;
            else
                % tiffn_i = fullfile(projectDir, sprintf(tiffn, str2double(timestr))) ;
                % Load the raw data from tiff
                % raw = readTiff4D(tiffn_i, 1, 'xyz') ;
                raw = QS.getCurrentData() ;
                IV = raw{1} ;

            end

            % Raw and APDV aligned Meshes
            amesh = QS.loadCurrentAlignedMesh() ;

            % scaled mesh vertices
            xa = amesh.v(:, 1)  ;
            ya = amesh.v(:, 2)  ;
            za = amesh.v(:, 3)  ;

            % shift scaled mesh vertices
            if flipy 
                axyzshift = amesh.v - (QS.normalShift * resolution) * amesh.vn ;
            else
                axyzshift = amesh.v + (QS.normalShift * resolution) * amesh.vn ;
            end
            
            %% Plot with orthogonal sectioning (note: faceon gives face-on view(2))
            close all
            fig = figure('units', 'centimeters', ...
                'outerposition', [0 0 xwidth xwidth], 'visible', 'off') ;

            %% Interpolate the data and evaluate at rotated orthogonal planes     
            % PASSIVE ROTATION            
            % Obtain rotation matrix that undoes the APDV rotation 
            invRot = QS.invertRotation(rot) ;

            % Plane to image
            xMin = -50 ;
            xMax = 300 ;
            yMin = -90 ;
            yMax = 90 ;
            zMin = -90 ;
            zMax = 80 ;
            oversample_factor = 1.1 ;
            nX = round((xMax - xMin) / resolution * oversample_factor) ;
            nY = round((yMax - yMin) / resolution * oversample_factor) ; 
            nZ = round((zMax - zMin) / resolution * oversample_factor) ;
            x1 = x0 * ones(nZ, nY) ;
            y2 = -y0 * ones(nZ, nX) ;
            z3 = z0 * ones(nY, nX) ;
            xspace = linspace(xMin, xMax, nX) ;
            yspace = linspace(yMin, yMax, nY) ;
            zspace = linspace(zMin, zMax, nZ) ;
            [y1, z1] = meshgrid(yspace, zspace) ; 
            [x2, z2] = meshgrid(xspace, zspace) ; 
            [x3, y3] = meshgrid(xspace, yspace) ; 
            xyzX = (invRot * (([x1(:), y1(:), z1(:)]/resolution) - trans)')' ;
            xyzY = (invRot * (([x2(:), y2(:), z2(:)]/resolution) - trans)')' ;
            xyzZ = (invRot * (([x3(:), y3(:), z3(:)]/resolution) - trans)')' ;
            xpx = reshape(xyzX(:, 1), [nY, nZ]) ;
            xpy = reshape(xyzX(:, 2), [nY, nZ]) ; 
            xpz = reshape(xyzX(:, 3), [nY, nZ]) ;
            ypx = reshape(xyzY(:, 1), [nX, nZ]) ;
            ypy = reshape(xyzY(:, 2), [nX, nZ]) ; 
            ypz = reshape(xyzY(:, 3), [nX, nZ]) ;
            zpx = reshape(xyzZ(:, 1), [nX, nY]) ;
            zpy = reshape(xyzZ(:, 2), [nX, nY]) ; 
            zpz = reshape(xyzZ(:, 3), [nX, nY]) ;
            assert(all(size(x1) == size(y1)))
            assert(all(size(y1) == size(z1)))
            assert(all(size(x2) == size(y2)))
            assert(all(size(y2) == size(z2)))
            assert(all(size(x3) == size(y3)))
            assert(all(size(y3) == size(z3)))

            % preview with mesh
            if preview
                clf
                h = trimesh(rmesh.f, xr, yr, zr, 'EdgeColor', 'none', 'facecolor', lscolor) ;
                axis equal
                hold on;
                scatter3(xpx(:), xpy(:), xpz(:), 10, 'markeredgecolor', 'r')
                scatter3(ypx(:), ypy(:), ypz(:), 10, 'markeredgecolor', 'g')
                scatter3(zpx(:), zpy(:), zpz(:), 10, 'markeredgecolor', 'b')
                waitfor(gcf)
            end
            disp('performing interpolation in Xplane')
            ix = interp3(double(IV), xyzX(:, 2), xyzX(:, 1), xyzX(:, 3)) ;
            ix = reshape(ix, size(x1)) ;

            disp('performing interpolation in Yplane')
            iy = interp3(double(IV), xyzY(:, 2), xyzY(:, 1), xyzY(:, 3)) ;
            iy = reshape(iy, size(x2)) ;

            disp('performing interpolation in Zplane')
            iz = interp3(double(IV), xyzZ(:, 2), xyzZ(:, 1), xyzZ(:, 3)) ;
            iz = reshape(iz, size(x3)) ;

            if preview
                % Check order of axes wrt IV
                figure ;
                subplot(2, 2, 1)
                imagesc(1:size(IV, 2), 1:size(IV, 3), squeeze(IV(round(mean(xyzX(:, 1))), :, :))')
                xlabel('size(IV, 2)')
                ylabel('size(IV, 3)'); axis equal; title('IX (approx, pre-rot)')
                subplot(2, 2, 2)
                imagesc(1:size(IV, 1), 1:size(IV, 3), squeeze(IV(:, round(mean(xyzY(:, 2))), :))')
                xlabel('size(IV, 1)')
                ylabel('size(IV, 3)'); axis equal; title('IY (approx, pre-rot)')
                subplot(2, 2, 3)
                imagesc(1:size(IV, 1), 1:size(IV, 2), squeeze(IV(:, :, round(mean(xyzZ(:, 3)))))')
                xlabel('size(IV, 1)')
                ylabel('size(IV, 2)'); axis equal; title('IZ (approx, pre-rot)')           


                % check interpolation
                figure; 
                subplot(2, 2, 1); imagesc(yspace, zspace, ix'); axis equal ;
                subplot(2, 2, 2); imagesc(xspace, zspace, iy'); axis equal ;
                subplot(2, 2, 3); imagesc(xspace, yspace, iz'); axis equal ;
                waitfor(gcf)
                close all
            end

            ix(isnan(ix)) = 0 ;
            iy(isnan(iy)) = 0 ;
            iz(isnan(iz)) = 0 ;
            ix = brighten * ix / 2^(16) ;
            iy = brighten * iy / 2^(16) ;
            iz = brighten * iz / 2^(16) ;
        end
        
        %% Plot in APDV
        if ~exist(fnP, 'file') || overwrite 
            clf            
            hold on;
            surface(x2, -y2, z2, iy, 'FaceColor','texturemap', ...
                'EdgeColor','none','CDataMapping','direct', ...
                                'AmbientStrength', ambientStrength_meshOrtho)
            hold on;
            
            % Plot X slice (optional)
            if ~faceon && plotXslice                
                surface(x1, -y1, z1, ix, 'FaceColor','texturemap', ...
                    'EdgeColor','none','CDataMapping','direct',...
                                'AmbientStrength', ambientStrength_meshOrtho)
            end
            
            % Plot Z slice (if not looking laterally)
            if ~faceon
                surface(x3, -y3, z3,  iz, 'FaceColor','texturemap', ...
                   'EdgeColor','none','CDataMapping','direct',...
                                'AmbientStrength', ambientStrength_meshOrtho)
                hold on;
            end
            colormap bone
            h = trimesh(amesh.f, xa, ya, za, 'EdgeColor', 'none', ...
                'facecolor', lscolor) ;
            axis equal
            hold on;
            axis off

            lighting gouraud    % preferred method for lighting curved surfaces
            material dull    % set material to be dull, no specular highlights

            % Lighting and view 
            % view(-20, 20)
            view(viewAngles) 
            lgt = camlight('headlight') ;

            % Obtain xyzlimits
            if first
                xlims = xlim ;
                ylims = ylim ;
                zlims = zlim ;
                first = false ;
            end
            xlim(xlims) ;
            ylim(ylims) ;
            zlim(zlims) ;
            disp(['Saving ' fnP])

            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

            % saveas(fig, fn) ;
            export_fig(fnP, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fnP);
            
            % Face-On view
            view([0,0]) 
            disp(['saving face-on view: ' fnF])
            export_fig(fnF, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fnF);
            close all
            
        end 
        
        
        %% Plot x-section only
        if ~exist(fnX, 'file') || overwrite 
            close all
            
            fig = figure('units', 'centimeters', ...
                'outerposition', [0 0 xwidth xwidth], 'visible', 'off') ;
           
            hold on;
            % surface(x2, -y2, z2, iy, 'FaceColor','texturemap', ...
            %     'EdgeColor','none','CDataMapping','direct', ...
            %     'AmbientStrength', ambientStrength_meshOrtho)
            
            imagesc(x2(1, :), z2(:, 1), iy) ;
            % imagesc(iy) ;
            colormap bone   
            hold on;
            
            % Get boundary pixels via inpolyhedron
            xyzgrid = [x2(:), -y2(:), z2(:)] ;
            
            
            amesh_raw = QS.loadCurrentAlignedMesh() ;
            amesh_raw.v = amesh_raw.v - ...
               0.5 * QS.normalShift * amesh_raw.vn * QS.APDV.resolution ;
            inp = inpolyhedron(amesh_raw.f(:, [2, 1, 3]), ...
                amesh_raw.v, xyzgrid) ;
            ba = false(size(iy)) ;
            ba(inp) = true ;
            xz_outlines = bwboundaries(ba) ;
            for pp = 1:length(xz_outlines)
                outl = xz_outlines{pp} ;
                plot(diff(x2(1, 1:2)) * outl(:, 2) + min(x2(1, :)), ...
                     diff(z2(1:2, 1)) * outl(:, 1) + min(z2(1, :)), ...
                     '-', 'color', lscolor, 'linewidth', lwX) ;
            end
            
            axis equal
            axis off

            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

            % Face-On view
            % view([0,0]) 
            disp(['saving face-on outline: ' fnX])
            export_fig(fnX, '-r300', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fnF);
            close all
            
        end
        
        %% Create matching Texturepatch Image
        if ~exist(fn_textureP, 'file') ||  ~exist(fn_textureF, 'file') || overwrite
            
            % raw shifted mesh vertices  
            rmesh = QS.loadCurrentRawMesh() ;          
            rmesh.v = rmesh.v + QS.normalShift * rmesh.vn ;
            xr = rmesh.v(:, 1) ;
            yr = rmesh.v(:, 2) ;
            zr = rmesh.v(:, 3) ;
            
            % Psize is the linear dimension of the grid drawn on each triangular face
            Options = struct() ;
            Options.PSize = 8;
            Options.EdgeColor = 'none';
            %Options.Rotation = rot ;
            %Options.Translation = trans ;
            %Options.Dilation = resolution ;
            % OUTWARD, INWARD layers
            Options.numLayers = [5, 2];  
            % Note: [2, 2] marches ~0.524 um in either dir if layer spacing is 1.0
            Options.layerSpacing = 1 ;
            
            fig = figure('Visible', 'Off') ;
            disp(['creating texture patch ' num2str(tp, '%06d')])

            % Allow for axis swapping
            TV = rmesh.v(:, texture_axis_order) ;
            texture_patch_3d( amesh.f, axyzshift, rmesh.f, TV, ...
                 QS.getCurrentData(), Options );
            colormap bone
             
            % format the figure
            disp('formatting figure...')
            axis equal
            [~,~,~,xyzlim] = QS.getXYZLims() ;
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            view(viewAngles)
            
            titlestr = ['$t = $' num2str(tp*QS.timeInterval-t0) ' ' QS.timeUnits] ;
            title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
            % xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex')
            % ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
            % zlabel('DV position [$\mu$m]', 'Interpreter', 'Latex')
            axis off
            
            % Set size of figure to be same as orthoview figure
            fig = gcf ;
            set(fig, 'Units', 'centimeters');
            set(fig, 'Position', [0 0 xwidth ywidth]);

            export_fig(fn_textureP, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_textureP);
            
            % Face-On view
            view([0,0]) 
            disp(['saving face-on view: ' fn_textureF])
            export_fig(fn_textureF, '-nocrop', '-transparent')            
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_textureF);
            
            %% Add orthoviews to texturepatch image
            hold on;
            surface(x2, -y2, z2, iy, 'FaceColor','texturemap', ...
                'EdgeColor','none','CDataMapping','direct',...
                'AmbientStrength', ambientStrength)
            hold on;
            
            % Plot X slice (optional)
            if ~faceon && plotXslice                
                surface(x1, -y1, z1, ix, 'FaceColor','texturemap', ...
                    'EdgeColor','none','CDataMapping','direct',...
                    'AmbientStrength', ambientStrength)
            end
            
            % Plot Z slice (if not looking laterally)
            if ~faceon
                surface(x3, -y3, z3,  iz, 'FaceColor','texturemap', ...
                   'EdgeColor','none','CDataMapping','direct',...
                   'AmbientStrength', ambientStrength)
                hold on;
            end
            
            % Set size of figure to be same as orthoview figure
            fig = gcf ; 
            set(fig, 'Units', 'centimeters');
            set(fig, 'Position', [0 0 xwidth ywidth]);

            view(viewAngles)
            export_fig(fn_texture2P, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_texture2P);
            
            % Face-On view
            view([0,0]) 
            disp(['saving face-on view: ' fn_texture2F])
            export_fig(fn_texture2F, '-nocrop', '-transparent')
            % patchIm = getframe(gca);
            % imwrite( patchIm.cdata, fn_texture2F);
            close all
        end
    end
    
end
disp('done')
