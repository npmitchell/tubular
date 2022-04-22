function getMeshes(tubi, detectOptions)
    % Obtain mesh surfaces of volumetric data (like ImSAnE's surface
    % detection methods), here using integralDetector methods for
    % activecontouring
    if isfield(detectOptions, 'pressure')
        pressure = detectOptions.pressure ;
    else
        pressure = 0 ;
    end
    if isfield(detectOptions, 'tension')
        tension = detectOptions.tension ;
    else
        tension = 0 ;
    end
    if isfield(detectOptions, 'pre_pressure')
        pre_pressure = detectOptions.pre_pressure ;
    else
        pre_pressure = 0 ;
    end
    if isfield(detectOptions, 'pre_tension')
        pre_tension = detectOptions.pre_tension ;
    else
        pre_tension = 0 ;
    end
    if isfield(detectOptions, 'post_pressure')
        post_pressure = detectOptions.post_pressure ;
    else
        post_pressure = 0 ;
    end
    if isfield(detectOptions, 'post_tension')
        post_tension = detectOptions.post_tension ;
    else
        post_tension = 0 ;
    end
    opts = detectOptions ;



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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Identify the surface using the loaded probabilities here
    % Convert the current image to a level set using morphological snakes
    ssfactor = opts.ssfactor ;
    niter = opts.niter ;
    niter0 = opts.niter0 ;
    mslsDir = opts.mslsDir ;
    pressure = opts.pressure ;
    tension = opts.tension ;
    exit_thres = opts.exit_thres ;
    ofn_ply = opts.ofn_ply ;
    ofn_ls = opts.ofn_ls ;
    channel = opts.channel ;
    ms_scriptDir = opts.ms_scriptDir ;
    % tpstamp = fileMeta.timePoints(first_tp);
    timepoint = opts.timepoint ;
    ofn_smoothply = opts.ofn_smoothply ;
    mlxprogram = opts.mlxprogram;
    init_ls_fn = opts.init_ls_fn ;
    % Run MS on a series of timepoints in a parent directory
    run_full_dataset = opts.run_full_dataset ;
    % Radius of initial guess if init_ls_fn does not exist or is
    % not supplied
    radius_guess = opts.radius_guess ;
    save = opts.save ;
    center_guess = opts.center_guess ;
    plot_mesh3d = opts.plot_mesh3d ;
    dtype = opts.dtype ; 
    mask = opts.mask ;
    use_pointcloud = opts.mesh_from_pointcloud ;
    dataset_prob_searchstr = opts.prob_searchstr ;
    ilastikaxisorder = opts.ilastikaxisorder ;
    smooth_with_matlab = opts.smooth_with_matlab ;
    % since python flips axes wrt MATLAB, flip them here
    morphsnakesaxisorder = fliplr(ilastikaxisorder) ;

    % Create the output dir if it doesn't exist
    if ~exist(mslsDir, 'dir')
        mkdir(mslsDir)
    end

    outputMesh = fullfile(mslsDir, msls_mesh_outfn) ;

    % Check if previous time point's level set exists to use as a seed
    % First look for supplied fn from detectOptions.
    % If not supplied (ie init_ls_fn is none or empty string, then
    % seek previous timepoint output from MS algorithm.

    disp(['init_ls_fn = ', init_ls_fn])
    disp(['ofn_ls = ', ofn_ls])
    if strcmp(init_ls_fn, 'none') || strcmp(init_ls_fn, '')
        % User has NOT supplied fn from detectOptions
        init_ls_fn = [ofn_ls, ...
            num2str(timepoint - 1, '%06d' ) '.' dtype] ;
    end

    disp([ 'initial level set fn = ', init_ls_fn])
    if exist(init_ls_fn, 'file')
        % It does exist, and the given name is the RELATIVE path.    
        % Use it as a seed (initial level set) 
        disp('running using initial level set')
        init_ls = load(init_ls_fn) ;
    elseif exist(fullfile(mslsDir, init_ls_fn), 'file')
        % It does exist, and given name is the relative path
        % without the extension. 
        % Use it as a seed (initial level set)
        disp('running using initial level set')
        init_ls = load(fullfile(mslsDir, init_ls_fn)) ;
    elseif exist(fullfile(mslsDir, [ init_ls_fn '.h5']), 'file')
        % It does exist, and given name is the FULL path 
        % without the extension.
        disp('running using initial level set')
        init_ls = load(fullfile(mslsDir, init_ls_fn)) ;                
    else
        % The guess for the initial levelset does NOT exist, so use
        % a sphere for the guess.
        disp(['Using default sphere for init_ls -- no such file on disk: ' fullfile(mslsDir, [ init_ls_fn '.h5'])])
        init_ls = [] ;
    end

    % Flip axis order LR of output mesh to return to MATLAB
    % orientation? NO, not helpful to do "-permute_mesh 'zyx'"
    % since already did morphsnakesaxisorder = fliplr() earlier
    % command = [command ' -adjust_for_MATLAB_indexing'] ;

    disp(['does outputMesh exist: ', num2str(exist(outputMesh, 'file'))])
    disp(outputMesh)
    % error('here')
    if ~exist(outputMesh, 'file')

        if use_dataset_command
            % User has elected to run as a dataset, so pass a directory
            % with _Probabilities.h5 files to run on.
            error('handle here')
            msls_mesh_outfn = ofn_ply;
            ls_outfn = ofn_ls;
        else
            % We are running MS on a single file. Give the filename of
            % the ilatik output run on filename.h5
            prob_infn = [opts.fileName, '_Probabilities.h5'] ;
            msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d'), '.ply'];
            ls_outfn = [ofn_ls, num2str(timepoint, '%06d'), '.', dtype];

            BW = activecontour(data, init_ls, niter,...
                'SmoothFactor', tension, 'ContractionBias', pressure) ;

            % Convert BW to mesh
            mesh = isosurface(BW, 0.5) ;
        end

    else
        disp(['output PLY already exists: ', msls_mesh_outfn])
    end

    %% Clean up mesh file for this timepoint using MeshLab --------

    % Here use the boundary mesh from marching cubes to make a
    % smooth mesh
    % Check if we need to smooth the full dataset of meshes
    if use_dataset_command
        disp('Using dataset command for morphsnakes...')
        % find all ms_...ply files in mslsDir, and smooth them all
        files_to_smooth = dir(fullfile(mslsDir, [ofn_ply '*.ply'])) ;
        lsfns_to_smooth = dir(fullfile(mslsDir, [ls_outfn '*' dtype])) ;

        if abs(length(lsfns_to_smooth) - length(files_to_smooth))>1 
            error('The number of output levelsets does not equal the number of output unsmoothed PLYs. These must match.')
        end
        for i=1:length(files_to_smooth)
            msls_mesh_outfn = files_to_smooth(i).name ;
            infile = fullfile( mslsDir, msls_mesh_outfn );
            % Note that LS file is outputLs ;
            split_fn = strsplit(msls_mesh_outfn, ofn_ply) ;
            extension_outfn = split_fn{2} ;
            base_outfn_for_pointcloud = ofn_ply ;
            mesh_outfn = [ofn_smoothply, extension_outfn] ;
            outputMesh = fullfile(mslsDir, mesh_outfn);

            disp(['outputMesh = ', outputMesh])
            %bad = so_bad
            if ~exist( outputMesh, 'file')
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
                    mesh = read_ply_mod(infile) ;

                    % Check that this behaves the way we want
                    % mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                    % size(mesh.vn)
                    % size(mesh.v)
                    % assert(size(mesh.vn, 1) == size(mesh.v, 1))

                    newV = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], smooth_with_matlab) ;
                    disp('Compute normals...')
                    mesh.v = newV ;
                    mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                    disp(['Saving smoothed mesh to ' outputMesh])
                    plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
                end
            else
                disp(['t=', num2str(timepoint) ': smoothed mesh file found...'])
                disp([' --> file to smooth was ' files_to_smooth(i).name])
            end

        end
    else
        disp('Using individual timepoint command for morphsnakes')
        mesh_outfn = [ofn_smoothply, num2str(timepoint, '%06d'), '.ply'];
        outputMesh = fullfile(mslsDir, mesh_outfn) ;
        if ~exist( outputMesh, 'file')
            msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d' ), '.ply'];
            infile = fullfile( mslsDir, msls_mesh_outfn );
            if smooth_with_matlab < 0
                % Build meshlab command to smooth meshes
                command = ['meshlabserver -i ' infile ' -o ' outputMesh, ...
                    ' -s ' mlxprogram ' -om vn'];
                % Either copy the command to the clipboard
                clipboard('copy', command);
                % or else run it on the system
                disp(['running ' command])
                system(command)
            elseif smooth_with_matlab == 0
                disp('No smoothing, with either matlab or meshlab')
                 mesh = read_ply_mod(infile) ;
                 disp('Compute normals...')
                 mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                 plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)      
            elseif smooth_with_matlab > 0 
                disp(['Smoothing with MATLAB using lambda = given value of ' num2str(smooth_with_matlab)])
                mesh = read_ply_mod(infile) ;

                % Check that this behaves the way we want
                % mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                % size(mesh.vn)
                % size(mesh.v)
                % assert(size(mesh.vn, 1) == size(mesh.v, 1))

                newV = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], smooth_with_matlab) ;
                disp('Compute normals...')
                mesh.v = newV ;
                mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                disp(['Saving smoothed mesh to ' outputMesh])
                plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
            end
        else
            disp(['t=', num2str(timepoint) ': smoothed mesh file found, loading...'])    
        end    
    end