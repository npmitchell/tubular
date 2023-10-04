function cleanCylMeshes(tubi, cleanCylOptions) 
% cleanCylMeshes(tubi, cleanCylOptions)
%
% Parameters
% ----------
% tubi : TubULAR class instance
% cleanCylOptions : optional struct with optional fields
%   overwrite : bool
%   save_ims : bool
%
% Returns
% -------
%
% Saves to disk
% -------------
% tubi.fileName.aBoundaryDorsalPtsClean ;
% tubi.fileName.pBoundaryDorsalPtsClean 
% [tubi.fileName.aBoundaryDorsalPtsClean(1:end-3) '.txt'] ;
% [tubi.fileName.pBoundaryDorsalPtsClean(1:end-3) '.txt'] ;
%
% NPMitchell 2020

% Unpack options
if nargin < 2
    cleanCylOptions = struct() ;
end
overwrite = false ;
save_ims = true ;
if isfield(cleanCylOptions, 'overwrite')
    overwrite = cleanCylOptions.overwrite ;
end
if isfield(cleanCylOptions, 'save_ims')
    save_ims = cleanCylOptions.save_ims ;
end

% Define filepaths. Some of these are less pretty to be
% compatible with Windows operating systems now.
cylMeshCDir =  tubi.dir.cylinderMeshClean ;
cylMeshCBase = tubi.fileBase.cylinderMeshClean ;

cylMeshDir = tubi.dir.cylinderMesh ;
cylMeshBase = tubi.fileBase.cylinderMesh ;

dpFile = tubi.fileName.apBoundaryDorsalPts ;
outadIDxfn = tubi.fileName.aBoundaryDorsalPtsClean ;
outpdIDxfn = tubi.fileName.pBoundaryDorsalPtsClean ;

dir2make = fullfile( tubi.dir.cylinderMeshClean, 'images') ;
if ~exist(dir2make, 'dir')
    mkdir(dir2make)
end

cylinderMeshCleanFigBase = fullfile( tubi.dir.cylinderMeshClean, 'images', ...
    [ tubi.fileBase.mesh 'cylindercut_clean.png' ]) ;

[~, ~, xyzlim_um] = tubi.getXYZLims() ;

% Load or compute clean cylindrical mesh for each timepoint
for t = tubi.xp.fileMeta.timePoints 
    disp(['t = ' num2str(t)])
    mesh3dfn = fullfile(tubi.dir.cylinderMeshClean,...
            sprintf( cylMeshCBase, t )) ;
    if ~exist(mesh3dfn, 'file') || overwrite
        disp('Overwriting/Computing clean cylinderMesh')

        % Load the cylinder mesh
        cylmeshfn = fullfile(cylMeshDir, sprintf( cylMeshBase, t )) ;
        cylmesh = read_ply_mod( cylmeshfn );
        
        % Clean the mesh, including averaging normals onto vertices using
        % angle-weighted scheme
        mesh = cleanCylMesh(cylmesh) ;   

        % The dataset name base for the AD and PD points
        ADBase = ['/' tubi.fileBase.mesh '/adorsal'];
        PDBase = ['/' tubi.fileBase.mesh '/pdorsal'];
        
        % Point-match the AD and PD points to nearest clean-boundary pt
        [adIDx, pdIDx] = aux_adjust_dIDx(mesh, cylmesh, t, dpFile, ...
            ADBase, PDBase, cylMeshCDir, cylMeshCBase, ...
            outadIDxfn, outpdIDxfn, tubi.xp.fileMeta.timePoints) ;

        %% Save the 3d cut mesh with new indices
        % This is saving the cylinder meshes with no ears. Also adIDx
        % is saved in h5file.
        plywrite_with_normals(mesh3dfn, mesh.f, mesh.v, mesh.vn)
        disp('Saving a/pdIDx to disk')
        % Save adIDx with new indices
        save_to_h5(outadIDxfn, ['/' sprintf('%06d', t) ], adIDx, ['adIDx for t=' num2str(t) ' already exists, overwriting'])
        % Save pdIDx with new indices
        save_to_h5(outpdIDxfn, ['/' sprintf('%06d', t) ], pdIDx, ['pdIDx for t=' num2str(t) ' already exists, overwriting'])
        disp('done with cylindermesh cleaning')

        % View results --------------------------------------------------------
        mesh3dfigfn = sprintfm( cylinderMeshCleanFigBase, t ) ;
        if save_ims
            close all
            fig = figure('visible', 'off') ;
            xyzrs = tubi.xyz2APDV(mesh.v) ;
            trisurf( triangulation( mesh.f, xyzrs ), xyzrs(:, 2), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.5);
            hold on
            scatter3( xyzrs(adIDx,1), xyzrs(adIDx,2), xyzrs(adIDx,3), ...
                'filled', 'r' );
            scatter3( xyzrs(pdIDx,1), xyzrs(pdIDx,2), xyzrs(pdIDx,3), ...
                'filled', 'c' );
            hold off
            axis equal
            xlabel('x [\mum]')
            ylabel('y [\mum]')
            zlabel('z [\mum]')
            xlim(xyzlim_um(1, :)); 
            ylim(xyzlim_um(2, :)); 
            zlim(xyzlim_um(3, :));
            title(['cleaned cylinder mesh, t = ' num2str(t)])
            saveas(fig, mesh3dfigfn)
            close all
        end
    else
        disp('Loading from disk')
        mesh = read_ply_mod(mesh3dfn) ;
        adIDx = h5read(outadIDxfn, ['/' sprintf('%06d', t)]) ;
        pdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', t)]) ;

        % View results --------------------------------------------------------
        mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
        if (~exist(mesh3dfigfn, 'file') || overwrite) && save_ims
            [~, ~, ~, xyzlim] = tubi.getXYZLims() ;

            %% Contents of aux_plot_cleanCylMesh.m below
            close all
            fig = figure('visible', 'off') ;
            xyzrs = tubi.xyz2APDV(mesh.v) ;
            trisurf( triangulation( mesh.f, xyzrs ), xyzrs(:, 2), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            hold on
            scatter3( xyzrs(adIDx,1), xyzrs(adIDx,2), xyzrs(adIDx,3), ...
                'filled', 'r' );
            scatter3( xyzrs(pdIDx,1), xyzrs(pdIDx,2), xyzrs(pdIDx,3), ...
                'filled', 'c' );
            hold off
            axis equal
            xlabel('x [\mum]')
            ylabel('y [\mum]')
            zlabel('z [\mum]')
            xlim(xyzlim(1, :)); 
            ylim(xyzlim(2, :)); 
            zlim(xyzlim(3, :));
            title(['cleaned cylinder mesh, t = ' num2str(t)])
            saveas(fig, mesh3dfigfn)
            close all
        end
    end
end

% Collate all ad and pd points and save as txt file as well
adIDxs = zeros(size(tubi.xp.fileMeta.timePoints)) ;
pdIDxs = zeros(size(tubi.xp.fileMeta.timePoints)) ;
for t = tubi.xp.fileMeta.timePoints
    disp(['t = ', num2str(t)])
    tidx = tubi.xp.tIdx(t) ;
    adIDx = h5read(tubi.fileName.aBoundaryDorsalPtsClean, ['/' sprintf('%06d', t)]) ;
    pdIDx = h5read(tubi.fileName.pBoundaryDorsalPtsClean, ['/' sprintf('%06d', t)]) ;
    adIDxs(tidx) = adIDx ;
    pdIDxs(tidx) = pdIDx ;
end
fnA = [tubi.fileName.aBoundaryDorsalPtsClean(1:end-3) '.txt'] ;
fnP = [tubi.fileName.pBoundaryDorsalPtsClean(1:end-3) '.txt'] ;
write_txt_with_header(fnA, adIDxs', 'vertex index for anterior dorsal')
write_txt_with_header(fnP, pdIDxs', 'vertex index for posterior dorsal')


