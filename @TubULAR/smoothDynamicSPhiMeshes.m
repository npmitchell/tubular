function [v3dsmM, nsmM] = smoothDynamicSPhiMeshes(QS, options)
%SMOOTHDYNAMICSPHIMESHES(QS, options)
%   Using (s,phi) pullback cutmeshes, smooth coord system in time via
%   simple triangular pulse averaging of positions in embedding space
%   while preserving the pullback mesh. 
%   Note that the computed geodesic distance along the mesh is preserved in
%   the pullback mesh extent spcutMeshSm.u(:, 1) varies from (0, GeoLength)
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   Current QuapSlap object
% options : struct with fields
%   width : int, default=4
%       half-width of tripulse smoothing filter on output meshes
%   overwrite : bool, default=false
%       overwrite results on disk
%   preivew : bool, default=QS.plotting.preview
%       display intermediate results
%
% NPMitchell 2020

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
spcutMeshSmBase = QS.fullFileBase.spcutMeshSm ;
spcutMeshSmRSBase = QS.fullFileBase.spcutMeshSmRS ;
spcutMeshSmRSCBase = QS.fullFileBase.spcutMeshSmRSC ;
timePoints = QS.xp.fileMeta.timePoints ;
spcutMeshBase = QS.fullFileBase.spcutMesh ;
[rot, trans] = QS.getRotTrans() ;
resolution = QS.APDV.resolution ;
flipy = QS.flipy ;

%% Unpack options
% Default options
width = 4 ;
overwrite = false ;

if nargin < 2
    options = struct() ;
end
if isfield(options, 'width')
    width = options.width ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end


%% Prep filter
disp('Building tripulse filter equivalent to tripuls(-0.5:0.1:0.5)')
% made the width variable 2020-12-09, used to be 0:0.2:1, so 9 tp included
% --> ie width used to equal 6.
tripulse = linspace(0, 1, width) ;
tripulse = [tripulse, fliplr(tripulse(1:end-1))] ;
tripulse = tripulse ./ sum(tripulse(:)) ;
tripulse = reshape(tripulse, [length(tripulse), 1]) ;

%% Check if all already exist
redo_meshsmooth = overwrite ;
qq = 1 ;
disp('Checking if smoothed meshes already exist on file')
while ~redo_meshsmooth && qq < length(timePoints)
    tt = timePoints(qq) ;
    % Check that all meshes are saved--if not, declare that we must compute
    smfn = sprintf(spcutMeshSmBase, tt) ;
    smrsfn = sprintf(spcutMeshSmRSBase, tt) ;
    smrscfn = sprintf(spcutMeshSmRSCBase, tt) ;
    if ~exist(smfn, 'file') || ~exist(smrsfn, 'file') || ~exist(smrscfn, 'file')
        redo_meshsmooth = true ;
    end
    qq = qq + 1 ;
end

if redo_meshsmooth
    disp('Mesh smoothing does not exist on file, computing')
    % Load all spcutMesh objects 
    vM = zeros(length(timePoints), nU*nV, 3);
    nsmM = zeros(length(timePoints), nU*nV, 3) ;
    for i = 1:length(timePoints)
        tt = timePoints(i) ;
        % Load the spcutMesh for this timepoint
        disp(['Loading spcutMesh from disk... [t = ' num2str(tt) ']'])
        load(sprintf(spcutMeshBase, tt), 'spcutMesh') ;
        vM(i, :, :) = spcutMesh.v ;
    end
    disp('built v3d matrix')
    
    % Filter in time axis
    % linfilt = 0.1 * ones(10, 1, 1) ;
    % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
    v3dsmM = imfilter(vM, tripulse, 'replicate') ;
    % vsmM = permute(vsmM, [2,1,3]) ;
    % nsmM = permute(nsmM, [2,1,3]) ;

    % Alternative is to use Gaussian filter but seems no padding option here:
    % smoothdata(vM, 1, 'gaussian', 10)
    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the smoothed meshes, then smoothed/rotated/scaled meshes 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for qq = 1:length(timePoints)
        tt = timePoints(qq) ;
        load(sprintf(spcutMeshBase, tt), 'spcutMesh') ;

        % First build normals from smoothed vertices as facenormals
        vqq = squeeze(v3dsmM(qq, :, :)) ;
        nqq = per_vertex_normals(vqq, spcutMesh.f, 'Weighting', 'angle') ;
        % pack it into a struct mesh
        mesh.v = vqq ;
        mesh.vn = nqq ;
        mesh.f = spcutMesh.f ;

        % Average normals with neighboring normals
        nqq = average_normals_with_neighboring_vertices(mesh, 0.5) ;
        nsmM(qq, :, :) = nqq ;

        % Next save the mesh
        smfn = sprintf(spcutMeshSmBase, tt) ;
        smrsfn = sprintf(spcutMeshSmRSBase, tt) ;
        smrscfn = sprintf(spcutMeshSmRSCBase, tt) ;
        if ~exist(smfn, 'file') || ~exist(smrsfn, 'file') || ~exist(smrscfn, 'file') || overwrite
            vqq = squeeze(v3dsmM(qq, :, :)) ;
            nqq = squeeze(nsmM(qq, :, :)) ;
            nsmM(qq, :, :) = nqq ./ vecnorm(nqq, 2, 2) ;

            % rotate and scale
            nqqrs = (rot * nqq')' ;
            vqqrs = ((rot * vqq')' + trans) * resolution;
            if flipy 
                vqqrs(:, 2) = -vqqrs(:, 2) ;
                nqqrs(:, 1) = -nqqrs(:, 1) ;
                nqqrs(:, 3) = -nqqrs(:, 3) ;
            else
                disp('Keeping normals as they are. I think this is right...')
                % error('Should we flip normals here? I think we want normals to be outward in APDV coords?')
            end

            % check normals by plotting the meshes colored by normals
            % trisurf(spcutMesh.f, vqqrs(:, 1), vqqrs(:, 2), vqqrs(:, 3), ...
            %     nqqrs(:, 1), 'edgecolor', 'none')
            % pause(5)
            % clf
            % trisurf(spcutMesh.f, vqqrs(:, 1), vqqrs(:, 2), vqqrs(:, 3), ...
            %     nqqrs(:, 2), 'edgecolor', 'none')
            % pause(5)
            % clf
            % trisurf(spcutMesh.f, vqqrs(:, 1), vqqrs(:, 2), vqqrs(:, 3), ...
            %     nqqrs(:, 3), 'edgecolor', 'none')
            % pause(5)
            % clf
            
            spcutMeshSm.f = spcutMesh.f ;
            spcutMeshSm.v = vqq ;
            spcutMeshSm.vn = nqq ;
            spcutMeshSm.u = spcutMesh.sphi ;
            spcutMeshSm.nU = spcutMesh.nU ;
            spcutMeshSm.nV = spcutMesh.nV ;
            spcutMeshSm.pathPairs = spcutMesh.pathPairs ;

            % Compute relaxed aspect ratio
            tmp = spcutMeshSm.u ;
            tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
            arspsm = minimizeIsoarealAffineEnergy( spcutMeshSm.f, spcutMeshSm.v, tmp );
            spcutMeshSm.ar = arspsm ;

            % Resave s,phi smoothed mesh and their 3D embedding
            disp(['tp=' num2str(tt) ': Saving ' sprintf(spcutMeshSmBase, tt)])
            save(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;

            % Also save rotated and scaled (RS) copy of the time-smoothed mesh
            spcutMeshSmRS = spcutMeshSm ;
            spcutMeshSmRS.v = vqqrs ;
            spcutMeshSmRS.vn = nqqrs ;
            % Resave s,phi and their 3D embedding
            save(sprintf(spcutMeshSmRSBase, tt), 'spcutMeshSmRS') ;
            clearvars vqq vqqrs

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CLOSED & GLUED SMOOTHED MESHES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To close the mesh, do the following:
            tmp = spcutMeshSmRS ;
            tmp.u = spcutMesh.sphi ;
            spcutMeshSmRSC = glueCylinderCutMeshSeam(tmp) ;
            save(sprintf(spcutMeshSmRSCBase, tt), 'spcutMeshSmRSC') ;

            % check it
            % fig = figure ;
            % % triplot(spcutMeshSmRSC.f, spcutMeshSmRSC.u(:, 1), spcutMeshSmRSC.u(:, 2))
            % trisurf(spcutMeshSmRSC.f, spcutMeshSmRSC.v(:, 1),...
            %     spcutMeshSmRSC.v(:, 2), spcutMeshSmRSC.v(:, 3))
            % To fully close the mesh use:
            % anewpt = mean(spcutMeshSmRS.v(1:nV:end, :), 1)
            % pnewpt = mean(spcutMeshSmRS.v(nU:nV:end, :), 1)
            % waitfor(fig)
        end
    end
    disp('done smoothing meshes in time')
else
    disp('Mesh smoothing already exists on file, skipping')
end

