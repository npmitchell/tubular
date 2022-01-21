function plotAlignedMeshesPretty(QS, options)
% plotAlignedMeshesPretty(QS, options)
%   Plot aligned APDV meshes with pretty lighting for timeseries 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   overwrite : bool
%       overwrite previous images on disk
%
%   ApplyAmbientOcclusion: Determines if the texture
%       colors should be modified by the ambient occlusion of the 
%       underlying triangulation {'false'}
%
%   AmbientOcclusion:    #Vx1 list of ambient occlusion
%       values
%
%   AmbientOcculsionFactor: A scalar between [0,1]
%       determining how strongly the ambient occlusion modifies the 
%       texture mapping {1}
%
%   AmbientOcculsionSamples: The number of samples to use
%       when computing the ambient occlusion {1000}
%
%   Unoriented: Treat the surface as unoriented when
%       applying ambient occlusion {'false'}
%
%   AddLights: Will add lighting to the plot {'false'}
%
%
% Returns
% -------
% <none>
%
% Saves to disk
% -------------
% metadat.mat : mat file with variables 'metadat'
%   options for how plotting was performed
%
% NPMitchell 2021

%% Default options
overwrite = false ;
plot_dorsal = true ;
plot_ventral = true ;
plot_left = true ;
plot_right = true ;
plot_perspective = true ;
channel = [] ;  % by default, plot all channels
AOSamples = 1000 ;
normal_shift = QS.normalShift ;
blackFigure = false ;
smoothIter = 1000 ;
smoothing_lambda = 0.0 ;

%% Unpack Options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'channel')
    channel = options.channel ;
end
if isfield(options, 'plot_dorsal')
    plot_dorsal = options.plot_dorsal ;
end
if isfield(options, 'plot_ventral')
    plot_ventral = options.plot_ventral ;
end
if isfield(options, 'plot_dorsal')
    plot_left = options.plot_left ;
end
if isfield(options, 'plot_dorsal')
    plot_right = options.plot_right ;
end
if isfield(options, 'plot_perspective')
    plot_perspective = options.plot_perspective ;
end
if isfield(options, 'AOSamples')
    AOSamples = options.AOSamples ;
end
if isfield(options, 'normalShift')
    normal_shift = options.normalShift ;
elseif isfield(options, 'normal_shift')
    normal_shift = options.normal_shift ;
end
if isfield(options, 'blackFigure')
    blackFigure = options.blackFigure ;
end 
if isfield(options, 'smoothing_lambda')
    smoothing_lambda = options.smoothing_lambda ;
end

% Collate boolean plot indicators to decide which views to plot
plot_view = [plot_dorsal, plot_ventral, plot_left, ...
    plot_right, plot_perspective] ;

%% Unpack QS
meshFileBase = QS.fullFileBase.alignedMesh ;
figoutdir = fullfile(QS.dir.alignedMesh, 'images_pretty') ;
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir)
end
% texture_axis_order = QS.data.axisOrder ;

try
    t0 = QS.t0set() ;
catch
    t0 = QS.xp.fileMeta.timePoints(1) ;
end

%% Load metadat and TexturePatchOptions if not supplied
metafn = fullfile(figoutdir, 'metadat.mat') ;
resave_metadat = false ;

if nargin < 3
    try
        load(metafn, 'metadat')
    catch
        % Define & Save metadata
        [~, ~, ~, xyzbuff] = QS.getXYZLims() ;
        xyzbuff(:, 1) = xyzbuff(:, 1) - 20 ; 
        xyzbuff(:, 2) = xyzbuff(:, 2) + 20 ; 
        metadat.normal_shift = QS.normalShift ;             % normal push, in pixels, along normals defined in data XYZ space
        metadat.xyzlim = xyzbuff ;                          % xyzlimits
        metadat.reorient_faces = false ;                    % set to true if some normals may be inverted
        resave_metadat = true ;
    end
else
    metadat = options ;
end

% Save it
if resave_metadat
    save(metafn, 'metadat')
end

% pass metadat to options
options.normal_shift = metadat.normal_shift ;
options.xyzlim = metadat.xyzlim ;
texture_axis_order = options.texture_axis_order ;
options.reorient_faces = metadat.reorient_faces ;

%% Name output directories
figdDir = fullfile(figoutdir, 'dorsal') ;
figvDir = fullfile(figoutdir, 'ventral') ;
figlat1Dir = fullfile(figoutdir, 'lateral1') ;
figlat2Dir = fullfile(figoutdir, 'lateral2') ;
figPerspDir = fullfile(figoutdir, 'perspective') ;
dirs = {figoutdir, figdDir, figvDir, figlat1Dir, figlat2Dir, figPerspDir} ;
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i}) ;
    end
end

%% Define output filenames
fns = {fullfile(figdDir, 'patch_dorsal_%06d.png'), ...
    fullfile(figvDir, 'patch_ventral_%06d.png'), ...
    fullfile(figlat1Dir, 'patch_lateral1_%06d.png'), ...
    fullfile(figlat2Dir, 'patch_lateral2_%06d.png'), ...
    fullfile(figPerspDir, 'patch_persp_%06d.png') };

% Unpack metadat and save 
xyzlim = options.xyzlim ;
reorient_faces = options.reorient_faces ;
timeinterval = QS.timeInterval ;
timeunits = QS.timeUnits ;
timePoints = QS.xp.fileMeta.timePoints ;

% Unpack xyzlim
xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;

%% Now draw for all TPs
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

tidx_todoA = 1:20:length(timePoints) ;
tidx_todoB = setdiff(1:10:length(timePoints), tidx_todoA) ;
tidx_todo = [tidx_todoA, tidx_todoB] ;
tidx_todoC = setdiff(1:length(timePoints), tidx_todo) ;
tidx_todo = [tidx_todo, tidx_todoC] ;

for tidx = tidx_todo
    tp = timePoints(tidx) ;
    ondisk = true ;
    for ii = 1:length(fns)
        ondisk = ondisk && exist(sprintf(fns{ii}, tp), 'file') ;
    end
    
    if (overwrite || ~ondisk) && any(plot_view)
        tic 
        close all

        QS.setTime(tp)
        
        % % Get the data in 3d for this timepoint
        % tidx = QS.xp.tIdx(tp) ;
        % % if t ~= xp.currentTime
        % xp.loadTime(tp) ;
        % xp.rescaleStackToUnitAspect() ;
        % % end
        % 
        % disp('Applying image...')
        % IV = xp.stack.image.apply();
        % IV = QS.adjustIV(IV, adjustlow, adjusthigh) ;

        % Read in the mesh file -----------------------------------------------
        disp('Reading mesh...')
        % Specfiy the mesh file to load
        meshfn = sprintf( meshFileBase, tp );
        mesh = read_ply_mod( meshfn );

        % If we smooth before pushing along the normal
        if smoothing_lambda > 0 
            disp('smoothing mesh via laplacian filter')
            mesh.v = laplacian_smooth(...
                mesh.v, mesh.f, 'cotan', [], smoothing_lambda, 'implicit', mesh.v, smoothIter) ;
            
        end
        
        % Make sure vertex normals are normalized
        mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
        % Normally evolve vertices
        mesh.v = mesh.v + normal_shift .* mesh.vn;
        % Re-orient faces
        if reorient_faces
            disp('reorienting faces...')
            mesh.f = reorient_facets( mesh.v, mesh.f );
        end

        % First args are physical vertices, then texture faces (same number as 
        % physical faces, but connectivity can be different), then texture
        % vertices, which can also be different. The vertices are in row, column, 
        % page format, so x and y get flipped. IV is the texture volume.
        % Options.PSize 
        % mesh.f = f{tidx} ;
        % mesh.v = v{tidx} ;
        % mesh.vn = vn{tidx} ;

        fig = figure('Visible', 'Off') ;
        % set(gcf, 'Visible', 'Off') 
        disp(['creating pretty aligned mesh image ' num2str(tp, '%06d')])
        
        disp('Creating ambient occlusion')
        AO = ambient_occlusion(mesh.v, mesh.f, mesh.v, mesh.vn, AOSamples);
        disp('trisurfing')
        trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ... %'facecolor', 'w')
            'FaceVertexCData',  AO * [1, 1, 1])
            
        % format the figure
        disp('formatting figure...')
        axis equal
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        zlim([zmin, zmax])
        titlestr = ['$t = $' num2str(tp*timeinterval-t0) ' ' timeunits] ;
        
        if blackFigure
            title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
        else
            title(titlestr, 'Interpreter', 'Latex') 
        end
        xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('DV position [$\mu$m]', 'Interpreter', 'Latex')

        % Extract the center of the current axis bounding box
        disp('adding lights')
        cen = [ mean([xmin, xmax]) mean([ymin, ymax]) mean([zmin, zmax]) ];
        light( 'Position', [cen(1) cen(2) 10*(zmax-zmin)+cen(3)], ...
            'Style', 'local', 'Color', [1 1 1]/3 );
        light( 'Position', [cen(1) 10*(zmax-zmin)+cen(2) cen(3)], ...
            'Style', 'local', 'Color', [1 1 1]/3 );
        light( 'Position', [cen(1) 10*(zmin-zmax)+cen(2) cen(3)], ...
            'Style', 'local', 'Color', [1 1 1]/3 );


        % Rotate the camera angle using rotation and translation 
        % camorbit(theta, phi)
        % camup() --> ax.CameraUpVector = [sin(45) cos(45) 1]
        % camroll(dtheta)
        % This is an active rotation
        % rotate(hSurface,direction,25)
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

        % Make background black & Make tick labels white
        grid off
        
        if blackFigure
            set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
            set(gcf, 'InvertHardCopy', 'off');
            set(gcf, 'Color', 'k')
            set(gcf, 'color', 'k')
        end
        
        % Capture all four views
        disp(['saving figure...' num2str(tp, '%06d')])
        % Save each figure
        for ii = 1:length(fns)
            % Only plot this view if plot_view(ii) is true
            if plot_view(ii)
                if ii == 1
                    % dorsal
                    disp('saving dorsal image...')
                    view(0, 90)
                elseif ii == 2
                    % ventral
                    disp('saving ventral image...')
                    view(0, 270)
                elseif ii == 3
                    % Lateral views
                    disp('saving lateral image...')
                    view(0, 0)
                elseif ii == 4
                    % lateral view 2
                    disp('saving second lateral image...')
                    view(0, 180)
                elseif ii == 5
                    % perspective view
                    disp('saving perspective image...')
                    view(-20, 20)
                else
                    error(['Exhausted DorsalVentralLeftRight indices. ',...
                        'What is going on here?'])
                end

                % Use export_fig instead, from plotting/export_fig/
                % saveas(fig, fullfile(figvdir, fnv))
                export_fig(sprintf(fns{ii}, tp), '-nocrop', '-r200')
            end
        end
        close all
        toc
    end
end