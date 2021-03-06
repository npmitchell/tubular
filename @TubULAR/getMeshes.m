function getMeshes(tubi, overwrite)
% Obtain mesh surfaces of volumetric data (like ImSAnE's surface
% detection methods), here using integralDetector methods for
% activecontouring
%
% Parameters
% ----------
% tubi : 
% overwrite : bool
%   overwrite previous meshes on disk
%   preview : bool
%   tension : 
%   prepressure
%   maxIterRelaxMeshSpikes : int >=0 (default=1000)
%       maximum number of iterations to relax spikes in the mesh
%
%
% Returns
% -------
% none
%
% Saves to disk
% -------------
% raw meshes as ply files 
% level sets as .mat files, 
% smoothed meshes as ply files
%
%
% NPMitchell 2022

if nargin < 2
    overwrite = false ;
end

opts = tubi.xp.detectOptions ;

if isfield(opts, 'preview')
    preview = opts.preview ;
else
    preview = true ;
end

if isfield(opts, 'pressure')
    pressure = opts.pressure ;
else
    pressure = 0 ;
end
if isfield(opts, 'tension')
    tension = opts.tension ;
else
    tension = 0 ;
end
if isfield(opts, 'pre_pressure')
    pre_pressure = opts.pre_pressure ;
else
    pre_pressure = 0 ;
end
% if isfield(opts, 'pre_tension')
%     pre_tension = opts.pre_tension ;
% else
%     pre_tension = 0 ;
% end
if isfield(opts, 'post_pressure')
    post_pressure = opts.post_pressure ;
else
    post_pressure = 0 ;
end
% if isfield(opts, 'post_tension')
%     post_tension = opts.post_tension ;
% else
%     post_tension = 0 ;
% end
if isfield(opts, 'target_edgelength')
    tar_length = opts.target_edgelength ;
else
    tar_length = 6 ;
end
if isfield(opts, 'enforceSingleComponent')
    enforceSingleComponent = opts.enforceSingleComponent ;
else
    enforceSingleComponent = true ;
end
if isfield(opts, 'enforceQuality')
    enforceQuality = opts.enforceQuality ;
else
    enforceQuality = false ;
end
if isfield(opts, 'maxIterRelaxMeshSpikes')
    maxIterRelaxMeshSpikes = opts.maxIterRelaxMeshSpikes ;
else
    maxIterRelaxMeshSpikes = 100 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify the surface using the loaded probabilities here
% Convert the current image to a level set using morphological snakes
ssfactor = opts.ssfactor ;
niter = opts.niter ;
niter0 = opts.niter0 ;
mslsDir = opts.mslsDir ;
exit_thres = opts.exit_thres ;
ofn_ply = opts.ofn_ply ;
ofn_ls = opts.ofn_ls ;
channel = opts.channel ;
% tpstamp = fileMeta.timePoints(first_tp);
timepoint = opts.timepoint ;
ofn_smoothply = opts.ofn_smoothply ;
mlxprogram = opts.mlxprogram;
init_ls_fn = opts.init_ls_fn ;
% Run MS on a series of timepoints in a parent directory
try
    dataDir = opts.dataDir;
catch
    dataDir = opts.run_full_dataset ;
end
% Radius of initial guess if init_ls_fn does not exist or is
% not supplied
radius_guess = opts.radius_guess ;
center_guess = opts.center_guess ;
plot_mesh3d = opts.plot_mesh3d ;
dtype = opts.dtype ; 
mask = opts.mask ;
use_pointcloud = opts.mesh_from_pointcloud ;
ilastikaxisorder = opts.ilastikaxisorder ;
smooth_with_matlab = opts.smooth_with_matlab ;

% Create the output dir if it doesn't exist
if ~exist(mslsDir, 'dir')
    mkdir(mslsDir)
end

if ~contains(ofn_ply, '%') || ~contains(ofn_ply, 'd')
    ofn_ply = [ofn_ply '%06d.ply'] ;
end

if ~contains(ofn_smoothply, '%') || ~contains(ofn_smoothply, 'd')
    ofn_smoothply = [ofn_smoothply '%06d.ply'] ;
end

for tidx = 1:length(tubi.xp.fileMeta.timePoints)
    tp = tubi.xp.fileMeta.timePoints(tidx) ;

    %% Define the mesh we seek
    outputLSfn = fullfile(mslsDir, sprintf([ofn_ls '%06d.' dtype], tp)) ;
    outputMesh = fullfile(mslsDir, sprintf(ofn_ply, tp)) ;

    disp(['does outputMesh exist: ', num2str(exist(outputMesh, 'file'))])
    disp(outputMesh)
    % error('here')
    if ~exist(outputMesh, 'file') || overwrite
        %% load the exported data out of the ilastik prediction
        fileName = fullfile(dataDir, ...
            [sprintf(opts.fileName, tp), '_Probabilities.h5']) ;
        disp(['Reading h5 file: ' fileName])
        h5fileInfo = h5info(fileName);
        if strcmp(h5fileInfo.Datasets.Name,'exported_data')
            file = h5read(fileName,'/exported_data');
        elseif strcmp(h5fileInfo.Datasets.Name,'volume')
            file = h5read(fileName,'/volume/prediction');
        else
            error(['Please provide a regular prediction from ilastik, either in', ...
                'the format of version 1.1 or 0.5 (ie with exported_data as a dataset)']);
        end

        % ilastik internally swaps axes. 1:x, 2:y, 3:z, 4:class
        % strategy: put into xyzc format, then pop last index
        if strcmp(opts.ilastikaxisorder, 'xyzc')
            pred = file ;
        elseif strcmp(opts.ilastikaxisorder, 'yxzc')
            % to convert yxzc to xyzc, put x=2 y=1 z=3 c=4
            pred = permute(file,[2,1,3,4]);
        elseif strcmp(opts.ilastikaxisorder, 'zyxc')
            % to convert yxzc to xyzc, put x=3 y=2 z=1 c=4
            pred = permute(file,[3,2,1,4]);
        elseif strcmp(opts.ilastikaxisorder, 'yzcx')
            % to convert yxzc to xyzc, put x=4 y=1 z=2 c=3
            pred = permute(file,[4,1,2,3]);
        elseif strcmp(opts.ilastikaxisorder, 'cxyz')
            % to convert yxzc to xyzc, put x=2 y=3 z=4 c=1
            pred = permute(file,[2,3,4,1]);
        elseif strcmp(opts.ilastikaxisorder, 'cyxz')
            % to convert yxzc to xyzc, put x=2 y=3 z=4 c=1
            pred = permute(file,[3,2,4,1]);
        elseif strcmp(opts.ilastikaxisorder, 'czyx')
            % to convert yxzc to xyzc, put x=1>4 y=2>3 z=3>2 c=4>1
            pred = permute(file,[4,3,2,1]);
        elseif strcmp(opts.ilastikaxisorder, 'cyzx')
            % to convert cyzx to xyzc put x=1>4 y=2>2 z=3>3 c=4>1
            pred = permute(file,[4,2,3,1]);
        elseif strcmp(opts.ilastikaxisorder, 'cxzy')
            % to convert cxzy to xyzc put x=1>2 y=2>4 z=3>3 c=4>1
            pred = permute(file,[2,4,3,1]);
        else
            error('Have not coded for this axisorder. Do so here')
        end


        % Check if previous time point's level set exists to use as a seed
        % First look for supplied fn from detectOptions.
        % If not supplied (ie init_ls_fn is none or empty string, then
        % seek previous timepoint output from MS algorithm.

        disp(['init_ls_fn = ', init_ls_fn])
        disp(['ofn_ls = ', ofn_ls])
        if tidx > 1 || strcmp(init_ls_fn, 'none') || strcmp(init_ls_fn, '')
            % User has NOT supplied fn from detectOptions
            init_ls_fn = [ofn_ls, ...
                num2str(tp - 1, '%06d' ) '.' dtype] ;
        end

        disp([ 'initial level set fn = ', init_ls_fn])
        if exist(init_ls_fn, 'file') || exist([init_ls_fn '.mat'], 'file') 
            % It does exist, and the given name is the RELATIVE path.    
            % Use it as a seed (initial level set) 
            disp('running using initial level set')
            init_ls = load(init_ls_fn, 'BW') ;
            init_ls = init_ls.BW ;
            niter_ii = niter ;
        elseif exist(fullfile(mslsDir, init_ls_fn), 'file') || ...
                exist(fullfile(mslsDir,[init_ls_fn '.mat']), 'file') 
            % It does exist, and given name is the relative path
            % without the extension. 
            % Use it as a seed (initial level set)
            disp('running using initial level set')
            init_ls = load(fullfile(mslsDir, init_ls_fn), 'BW') ;
            init_ls = init_ls.BW ;
            niter_ii = niter ; 
        else
            % The guess for the initial levelset does NOT exist, so use
            % a sphere for the guess.
            disp(['Using default sphere of radius ' num2str(radius_guess) ...
                ' for init_ls -- no such file on disk: ' fullfile(mslsDir, [ init_ls_fn '.h5'])])
            init_ls = zeros(size(squeeze(pred(:, :, :, opts.foreGroundChannel)))) ;
            SE = strel("sphere", radius_guess) ;
            SE = SE.Neighborhood ;
            se0 = size(SE, 1) ;
            rad = ceil(se0*0.5) ;
            assert( all(size(SE) == se0)) ;
            if isempty(center_guess)
                centers = size(pred) * 0.5 ;
                centers = centers(1:3) ;
            else
                centers = str2num(center_guess) ;
            end
            dd = centers - rad ;
            xmin = max(1, dd(1)) ;
            ymin = max(1, dd(2)) ;
            zmin = max(1, dd(3)) ;
            % xmax = min(size(init_ls, 1), dd(1)+se0);
            % ymax = min(size(init_ls, 2), dd(2)+se0) ;
            % zmax = min(size(init_ls, 3), dd(3)+se0) ;

            % distance from edge of data volume to edge of strel:
            sz0 = size(init_ls) ;

            % check if SE is entirely contained within init_ls
            if all(dd > 0)
                % minima are within the boundary
                init_ls(xmin:xmin+se0-1, ymin:ymin+se0-1, ...
                    zmin:zmin+se0-1) = SE ;
                init_ls = init_ls(1:sz0(1), 1:sz0(2), 1:sz0(3)) ;
            elseif all(dd + se0 < 0)
                % the maxima are all contained within the volume
                init_ls(xmin:xmin+se0+dd(1)-1, ...
                    ymin:ymin+se0+dd(2)-1, ...
                    zmin:zmin+se0+dd(3)-1) = ...
                    SE(-dd(1):se0, -dd(2):se0, -dd(3):se0) ;
                init_ls = init_ls(1:sz0(1), 1:sz0(2), 1:sz0(3)) ;
            else
                % the maxima are all contained within the volume
                chopx = dd(1) < 0 ;
                chopy = dd(2) < 0 ;
                chopz = dd(3) < 0 ;
                init_ls(xmin:xmin+se0+(dd(1)*chopx)-1, ...
                    ymin:ymin+se0+dd(2)*chopy-1, ...
                    zmin:zmin+se0+dd(3)*chopz-1) = ...
                        SE(-dd(1)*chopx+1:se0, ...
                        -dd(2)*chopy+1:se0, ...
                        -dd(3)*chopz+1:se0) ;
                init_ls = init_ls(1:sz0(1), 1:sz0(2), 1:sz0(3)) ;
            end

            try
                assert(any(init_ls(:)))
            catch
                error('The initial guess is outside the data volume')
            end

            niter_ii = niter0 ;
        end

        % Flip axis order LR of output mesh to return to MATLAB
        % orientation? NO, not helpful to do "-permute_mesh 'zyx'"
        % since already did morphsnakesaxisorder = fliplr() earlier
        % command = [command ' -adjust_for_MATLAB_indexing'] ;


        %% Extract contour/isosurface of levelset

        % Pre-processing
        if pre_pressure < 0
            disp('eroding input LS by post_pressure...')
            SE = strel('sphere', abs(pre_pressure)) ;
            init_ls = imerode(init_ls, SE) ;
        elseif pre_pressure > 0                
            disp('dilating input LS by post_pressure...')
            SE = strel('sphere', abs(pre_pressure)) ;
            init_ls = imdilate(init_ls, SE) ;
        end

        data = pred(:, :, :, opts.foreGroundChannel) ;
        % data_clipped = data - 0.1 ;
        % data_clipped(data_clipped < 0) = 0. ;

        disp(['niter is ', num2str(niter_ii)]);
        BW = activecontour(data, init_ls, 2*niter_ii, 'Chan-Vese', ...
            'SmoothFactor', tension, 'ContractionBias', -pressure) ;

        % Post processing
        if post_pressure < 0
            disp('eroding result by post_pressure...')
            SE = strel('sphere', abs(post_pressure)) ;
            BW = imerode(BW, SE) ;
        elseif post_pressure > 0          
            disp('dilating result by post_pressure...')      
            SE = strel('sphere', abs(post_pressure)) ;
            BW = imdilate(BW, SE) ;
        end

        % preview current results
        if preview
            clf
            if centers(3) < size(BW, 3) && centers(3) > 0.5 
                bwPage = squeeze(BW(:, :,round(centers(3)))) ;
                datPage = squeeze(data(:,:,round(centers(3)))) ;
                rgb = cat(3, bwPage, datPage, datPage) ;
                imshow(rgb)
                sgtitle('level set found...')
            else
                zframe = round(size(BW, 3) * 0.5) ;
                bwPage = squeeze(BW(:, :,zframe)) ;
                datPage = squeeze(data(:,:,zframe)) ;
                rgb = cat(3, bwPage, datPage, datPage) ;
                imshow(rgb)
                sgtitle('level set found...')

            end
            % pause to draw the figure and show it in foreground
            pause(1e-5)
        end

        % Remove all but biggest component
        if enforceSingleComponent
            CC = bwconncomp(BW, 6);
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [~,idx] = max(numPixels);
            BW = false(size(BW));
            BW(CC.PixelIdxList{idx}) = true;
            
            % bwareaopen(BW, PP, 6) ; % choose connectivity to be six so that face junctions are required
        end

        % Pad the walls with zeros
        BW2 = zeros(size(BW) + 2) ;
        BW2(2:end-1, 2:end-1, 2:end-1) = BW ;

        % Convert BW to mesh
        switch meshConstructionMethod
            case 'marchingCubes'
                mesh = isosurface(BW2, 0.5) ;
                mesh.vertices = mesh.vertices - 1 ;
                mesh.vertices = mesh.vertices * ssfactor ;
                % Swap X<->Y axes since MATLAB did this in isosurface/contour
                mesh.vertices = mesh.vertices(:, [2, 1, 3]) ;
                % mesh.faces = mesh.faces(:, [2, 1, 3]) ;
            case 'advancingFront'
                error('Implement advancing front from example here')
                
            case 'PoissonSurface'
                error('Implement poisson surface from example here')
                
        end
        
        % Check it
        % trisurf(triangulation(mesh.faces, mesh.vertices), 'edgecolor', 'none')

        % Write init_ls for next timepoint to disk
        save(outputLSfn, 'BW')

        % Write mesh to disk
        plywrite(outputMesh, mesh.faces, mesh.vertices)
    else
        disp(['output PLY already exists: ', ofn_ply])
    end

    %% Clean up mesh file for this timepoint using MeshLab --------

    % Here use the boundary mesh from marching cubes to make a
    % smooth mesh
    rawMesh = outputMesh ; 
    outputMesh = fullfile(mslsDir, sprintf(ofn_smoothply, tp));

    disp(['outputMesh = ', outputMesh])
    %bad = so_bad
    if ~exist( outputMesh, 'file') || overwrite
        % Smooth with either meshlab or matlab
        if smooth_with_matlab < 0
            % USE MESHLAB, not matlab
            if use_pointcloud
                % Use the pointcloud from the level set rather than the
                % boundary mesh from marching cubes
                %----------------------------------------------------------------------
                % Extract the implicit level set as a 3D binary array
                %----------------------------------------------------------------------

                % The file name of the current time point
                % The 3D binay array
                bwLS = h5read( outfnLS, '/implicit_levelset' );

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
                command = ['meshlabserver -i ' pointCloudFileName, ...
                    ' -o ' outputMesh, ' -s ' mlxprogram ' -om vn'] ;
                disp(['running ' command])
                system( command );
            else
                % Use the marching cubes mesh surface to smooth
                command = ['meshlabserver -i ' infile ' -o ' outputMesh, ...
                    ' -s ' mlxprogram ' -om vn'];
                % Either copy the command to the clipboard
                clipboard('copy', command);
                % or else run it on the system
                disp(['running ' command])
                system(command)
            end
        elseif smooth_with_matlab == 0
            disp('No smoothing, with either matlab or meshlab')
            mesh = read_ply_mod(infile) ;
            disp('Compute normals...')
            mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
            plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
        elseif smooth_with_matlab > 0 
            % Smooth with MATLAB
            disp(['Smoothing with MATLAB using lambda = given value of ' num2str(smooth_with_matlab)])
            mesh = read_ply_mod(rawMesh) ;

            % Check that this behaves the way we want
            % mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
            % size(mesh.vn)
            % size(mesh.v)
            % assert(size(mesh.vn, 1) == size(mesh.v, 1))

            if enforceSingleComponent
                [ F, V, oldVertexIDx, C ] = ...
                    remove_isolated_mesh_components(mesh.f, mesh.v) ;
            else
                F = mesh.f ;
                V = mesh.v ;
            end

            num_iter = 5;                
            protect_constraints = false;
            [V, F] = remesh_smooth_iterate(V,F, smooth_with_matlab,... 
                tar_length, num_iter, protect_constraints, enforceQuality, maxIterRelaxMeshSpikes) ;

            % newV = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], smooth_with_matlab) ;
            disp('Compute normals...')
            mesh.v = V ;
            mesh.f = F ;
            mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
            disp(['Saving smoothed mesh to ' outputMesh])
            plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
        end
    else
        disp(['t=', num2str(timepoint) ': smoothed mesh file found...'])
    end
end