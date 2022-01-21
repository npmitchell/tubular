function generateCurrentCutMesh(QS, cutMeshOptions)
% generateCurrentCutMesh(QS, cutMeshOptions)
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   The object for which we generate the currentTime's cutMesh
% cutMeshOptions : struct with fields
%   nsegs4path (int, optional, default=2)
%       How many segments of piecewise geodesics to draw to cut the
%       axisymmetric surface so that the cut does not change winding
%       number with respect to the previous timepoint's cut around the
%       centerline
%   maxJitter (float, optional, default=100)
%       maximum displacement in mesh units to randomly displace nodes of
%       the path in a brute-force search for topologically connected path
%   maxTwChange (float, optional, default=0.15)
%       maximum allowed change in the value of the twist of the cutPath
%       with respect to the centerline, as compared to the previous
%       timepoint's cutPath twist about the centerline. This is a soft
%       proxy for the topology change, but in practice it is superior to
%       the error-prone methods of measuring the (discrete-valued) winding
%       number explored thus far
%
% NPMitchell 2020

% Parameters for cutMesh creation
nsegs4path = 2 ;
maxJitter = 100 ;
maxTwChange = 0.7 ;
preview = false ;
definePDviaRicci_t0 = true ;
definePDviaRicci = false ;
ricciOptions = struct() ;
try
    t0 = QS.t0set() ;
catch
    t0 = QS.xp.fileMeta.timePoints(1) ;
end
nargin
if nargin > 1
    if isfield(cutMeshOptions, 'nsegs4path')
        nsegs4path = cutMeshOptions.nsegs4path ;
    end
    if isfield(cutMeshOptions, 'maxJitter')
        maxJitter = cutMeshOptions.maxJitter ;
    end
    if isfield(cutMeshOptions, 'maxTwChange')
        maxTwChange = cutMeshOptions.maxTwChange ;
    end
    if isfield(cutMeshOptions, 'preview')
        preview = cutMeshOptions.preview ;
    end
    if isfield(cutMeshOptions, 't0')
        t0 = cutMeshOptions.t0 ;
    end
    if isfield(cutMeshOptions, 'definePDviaRicci')
        definePDviaRicci = cutMeshOptions.definePDviaRicci ;
    end
    if isfield(cutMeshOptions, 'ricciOptions')
        ricciOptions = cutMeshOptions.ricciOptions ;
    end
end
ricciOptions.t0 = t0 ;

% Unpack parameters
tt = QS.currentTime ;
cutMeshfn = sprintf(QS.fullFileBase.cutMesh, tt) ;
QS.getCleanCntrlines() ;

% Used to save time if mesh is loaded, but this can lead to problems
% mesh = QS.currentMesh.cylinderMeshClean ;
% if isempty(mesh)
QS.loadCurrentCylinderMeshClean() ;
mesh = QS.currentMesh.cylinderMeshClean ;
% end

cylinderMeshCleanBase = QS.fullFileBase.cylinderMeshClean ;
outcutfn = QS.fullFileBase.cutPath ;
% centerlines from QS
QS.getCleanCntrlines ;
cleanCntrlines = QS.cleanCntrlines ;

% Grab ad and pd indices for cylinder mesh
adIDx = h5read(QS.fileName.aBoundaryDorsalPtsClean,...
    ['/' sprintf('%06d', tt)]) ;
pdIDx = h5read(QS.fileName.pBoundaryDorsalPtsClean,...
    ['/' sprintf('%06d', tt)]) ;

% try geodesic if first timepoint
if tt == t0
    cutPath_ok = false ;
    
    % Modify adIDx/pdIDx if desired (for removing possible twist)
    if definePDviaRicci_t0
        [rawRicciMesh, ~] = ...
            QS.generateRawRicciMeshTimePoint(tt, ricciOptions) ;

        % Assert that adIDx vertex is very near phi= 0
        assert(rawRicciMesh.rectangle.u(adIDx,2) == 0)
        
        % Find pIDx nearest to phi = 0
        bnd = freeBoundary(triangulation(mesh.f, mesh.v)) ;
        [bnds, hairpins] = separateClosedCurveComponents(bnd) ;
        try
            assert(~any(isempty(hairpins)))
        catch
            error('There should be no hairpins in cleaned cylinder meshes')
        end
        assert(length(bnds) == 2)
        
        % Find anterior dorsal vertex in one of the boundary rings
        if ismember(adIDx, bnds{1}(:, 1))
            pInds = bnds{2} ;
        elseif ismember(adIDx, bnds{2}(:, 1))
            pInds = bnds{1} ;
        else
            error('could not find anterior dorsal vertex index in either boundary list')
        end
        
        % Find min phi in posterior vertex indices
        [~, minId] = min(rawRicciMesh.rectangle.u(pInds, 2)) ;
        pdIDx = pInds(minId) ; 
    end
        
    dmyk = 0 ;
    bbWeight = 10 ;
    while ~cutPath_ok
        disp(['Attempting to cutMesh with bbWeight = ', num2str(bbWeight)])
        cutOptions.method = 'fastest' ;
        cutOptions.bbWeight = bbWeight ;
        disp(['Cutting mesh using method ' cutOptions.method])
        cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );
        cutP = cutMesh.pathPairs(:, 1) ;
        assert(adIDx == cutP(1)) ;
        assert(pdIDx == cutP(end)) ;

        bndy = freeBoundary(triangulation(mesh.f, mesh.v)) ;
        % Count number of cutP vertices that are on freeBoundary
        nPVertsOnBndy = length(intersect(bndy, cutP)) ;
        if nPVertsOnBndy == 2
            cutPath_ok = true ;
        else
            % increase bulk-boundary edge weight in graph G
            bbWeight = bbWeight * 5 ;
        end
        if dmyk > 10
            error('Could not make proper cutPath after 10 attempts')
        end
    end
    
    if preview 
        trisurf(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2), ...
            cutMesh.v(:, 3), 'edgecolor', 'none', 'facealpha', 1.)
        hold on;
        plot3(cutMesh.v(cutP, 1), cutMesh.v(cutP, 2), cutMesh.v(cutP, 3), '.-')
        title('Cut mesh with cutPath')
        disp('Showing cutMesh with cutPath. Close figure to continue')
        waitfor(gcf)
        % 
        % % Reduce the mesh and check the results this way
        % fv.faces = mesh.f ;
        % fv.vertices = mesh.v ;
        % fv = reducepatch(fv, 0.99) ;    
        % [~, adIDx2] = min(vecnorm(fv.vertices - mesh.v(adIDx, :), 2, 2)) ;
        % [~, pdIDx2] = min(vecnorm(fv.vertices - mesh.v(pdIDx, :), 2, 2)) ;
        % vn2 = per_vertex_normals(fv.vertices, fv.faces, 'Weighting', 'angle') ;
        % cutMesh2 = cylinderCutMesh(fv.faces, fv.vertices, vn2, adIDx2, pdIDx2, cutOptions );
        % bdyIDx = freeBoundary( triangulation( cutMesh2.f, cutMesh2.v ) ) ;
        % idx = bdyIDx(:, 1) ;
        % % plot the remeshing
        % trisurf(cutMesh2.f, cutMesh2.v(:, 1), cutMesh2.v(:, 2), ...
        %     cutMesh2.v(:, 3), 'edgecolor', 'none', 'facealpha', 1.)
        % hold on;
        % plot3(cutMesh2.v(idx, 1), cutMesh2.v(idx, 2), cutMesh2.v(idx, 3), '-')
        % title('Remeshing of cutMesh with resampling of cutPath')
        % disp('Showing remeshing of cutMesh with cutPath. Close figure to continue')
        % waitfor(gcf)
        
        % check normals to faces
        % nn = faceNormals(cutMesh.f, cutMesh.v)
        % cc = incenter(triangulation(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2), cutMesh.v(:, 3)))
        % quiver3(cc(:, 1), cc(:, 2), cc(:, 3), nn(:, 1), nn(:, 2), nn(:, 3), 1)
    end
else
    
    % Modify adIDx/pdIDx if desired (for removing possible twist)
    if definePDviaRicci
        [rawRicciMesh, ~] = ...
            QS.generateRawRicciMeshTimePoint(tt, ricciOptions) ;

        % Assert that adIDx vertex is very near phi= 0
        assert(rawRicciMesh.rectangle.u(adIDx,2) == 0)
        
        % Find pIDx nearest to phi = 0
        bnd = freeBoundary(triangulation(mesh.f, mesh.v)) ;
        [bnds, hairpins] = separateClosedCurveComponents(bnd) ;
        try
            assert(~any(isempty(hairpins)))
        catch
            error('There should be no hairpins in cleaned cylinder meshes')
        end
        assert(length(bnds) == 2)
        
        % Find anterior dorsal vertex in one of the boundary rings
        if ismember(adIDx, bnds{1}(:, 1))
            pInds = bnds{2} ;
        elseif ismember(adIdx, bnds{2}(:, 1))
            pInds = bnds{1} ;
        else
            error('could not find anterior dorsal vertex index in either boundary list')
        end
        
        % Find min phi in posterior vertex indices
        [~, minId] = min(rawRicciMesh.rectangle.u(pInds, 2)) ;
        pdIDx = pInds(minId) ; 
    end
       
    
    if tt > t0
        prevTP = tt - 1 ;
    elseif tt < t0 
        prevTP = tt + 1 ;
    end
    % If a previous Twist is not held in RAM, compute it
    % if ~exist('prevTw', 'var')
    % Load previous mesh and previous cutP
    prevcylmeshfn = sprintf( cylinderMeshCleanBase, prevTP) ;
    disp(['Loading previous cylinderMeshClean: ' prevcylmeshfn])
    prevmesh = read_ply_mod( prevcylmeshfn ); 
    
    disp(['Loading previous cutPath: ' sprintf(outcutfn, prevTP)])
    prevcutP = dlmread(sprintf(outcutfn, prevTP), ',', 1, 0) ;
    previousP = prevmesh.v(prevcutP, :) ;
    % Load previous centerline in raw units
    prevcline = cleanCntrlines{QS.xp.tIdx(prevTP)} ; % use previous CORRECTED centerline (non-anomalous)
    prevcline = prevcline(:, 2:4) ;
    % Compute Twist for this previous timepoint
    prevTw = twist(previousP, prevcline) ;

    % [edgelen, annulusv2d] = DiscreteRicciFlow.EuclideanRicciFlow(mesh.f, mesh.v, ...
    %     'BoundaryType', 'fixed', 'BoundaryShape', 'Circles', ...
    %     'MaxIter', 25, 'MaxCircIter', 21, ...
    %     'Tolerance', 1e-6, 'CircTolerance', 1e-4, 'PCGTolerance', 1e-4) ;
    % findAnnularPathZeroWindingNumber(mesh.f, annulusv2d, adIDx, pdIDx)

    % Which path to match this one to: choose previous timepoint
    % Load previous mesh and previous cutP
    prevcylmeshfn = sprintf( cylinderMeshCleanBase, prevTP) ;
    prevmesh = read_ply_mod( prevcylmeshfn ); 
    prevcutP = dlmread(sprintf(outcutfn, prevTP), ',', 1, 0) ;
    previousP = prevmesh.v(prevcutP, :) ;

    % Current centerline: chop off ss to make Nx3
    cntrline = cleanCntrlines{QS.xp.tIdx(tt)} ;
    cntrline = cntrline(:, 2:4) ;
    
    % Check topology 
    eulerChar = eulerCharacteristic(mesh) ;
    assert(eulerChar == 0)
    
    % Cut the mesh
    [cutMesh, adIDx, pdIDx, cutP, ~] = ...
        generateCutMeshFixedTwist(mesh, adIDx, pdIDx, ...
        cntrline,...  % supply the current corrected centerline
        nsegs4path, prevTw, previousP, ...
        'MaxTwChange', maxTwChange, 'MaxJitter', maxJitter, ...
        'PrevCntrline', prevcline, 'centerlineIsErratic', true) ;
end

% Store this path for the next one to be nearby
% Limit the number of segments to nsegs4path
% previousP = cutMesh.v(cutP, :) ;
% pstep = round(length(cutP) / nsegs4path ) ;
% previousP = previousP(1:pstep:end, :) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Done with generating initial 3D CutMesh with cutPath\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the cutPath to txt file
header = 'cutP (path of cut), indexing into vertices' ;
write_txt_with_header(sprintf(outcutfn, tt), cutP, header)  
clearvars header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate pullback to rectangular domain ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The surface parameterization algorithm (optionally) takes four vertex IDs
% as input to specify the corners of the square parameterization domain.
% Maddeningly, the order in which these points are specified to not seem to
% effect the output. For consistency, we perform a post-hoc correction so
% that the final output has the following geometric ordering
%
%   (AD1)-------(PD1)
%     |           |
%     |           |
%     |           |
%     |           |
%   (AD2)-------(PD2)
%
% Note that the pathPairs variable has the following columns:
%   [ ( AD1 -> PD1 ), ( AD2 -> PD2 ) ]
%--------------------------------------------------------------------------

% View results --------------------------------------------------------
% P = cutMesh.pathPairs(:,1);
% 
% trisurf( triangulation( mesh.f, mesh.v ) );
% 
% hold on
%
% line( mesh.v(P,1), mesh.v(P,2), mesh.v(P,3), ...
%     'Color', 'c', 'LineWidth',2);
% 
% scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
%     'filled', 'r' );
% scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
%     'filled', 'm' );
% 
% hold off
% 
% axis equal
% 
% clear P

%----------------------------------------------------------------------
% Generate Pullback to Annular Orbifold Domain
%----------------------------------------------------------------------
fprintf('Relaxing network via Affine transformation... ');
try
    cutMesh = flattenAnnulus( cutMesh );
catch
    disp('Bad boundary! Remeshing')
    cutMesh = flattenAnnulus( cutMesh );
    % Reduce the mesh and check the results this way
    fv.faces = mesh.f ;
    fv.vertices = mesh.v ;
    fv = reducepatch(fv, 0.9) ;    
    [~, adIDx2] = min(vecnorm(fv.vertices - mesh.v(adIDx, :), 2, 2)) ;
    [~, pdIDx2] = min(vecnorm(fv.vertices - mesh.v(pdIDx, :), 2, 2)) ;
    vn2 = per_vertex_normals(fv.vertices, fv.faces, 'Weighting', 'angle') ;
    cutMesh = cylinderCutMesh(fv.faces, fv.vertices, vn2, adIDx2, pdIDx2, cutMeshOptions );
    bdyIDx = freeBoundary( triangulation( cutMesh.f, cutMesh.v ) ) ;
    idx = bdyIDx(:, 1) ;
    
    % plot the remeshing
    if preview
        trisurf(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2), ...
            cutMesh.v(:, 3), 'edgecolor', 'none', 'facealpha', 1.)
        hold on;
        plot3(cutMesh.v(idx, 1), cutMesh.v(idx, 2), cutMesh.v(idx, 3), '.-')
        title('Remeshed the surface and took new cutPath. Will retry flattening')
        disp('Remeshed the surface and took new cutPath. Close figure to retry flattening')
        waitfor(gcf)
    end
    
    cutMesh = flattenAnnulus( cutMesh );
end

if preview
    figure('visible', 'on');
    trisurf(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2), ...
        cutMesh.v(:, 3), 'edgecolor', 'none')
    axis equal
    figure('visible', 'on');
    trisurf(cutMesh.f, cutMesh.u(:, 1), cutMesh.u(:, 2), ...
        0*cutMesh.u(:, 2), 'facecolor', 'none')
    waitfor(gcf)
end

% Find lateral scaling that minimizes spring network energy
ar = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
% Assign scaling based on options: either a0 or a_fixed
% if tidx == 1 && ~a_fixed
%     a_fixed = ar ;
% end      
% a = a_fixed ;

% Scale the x axis by a or ar, also flip u if flipy is true
uvtx = cutMesh.u ;
if QS.flipy
    cutMesh.u = [ QS.a_fixed .* uvtx(:,1), 1.0 - uvtx(:,2) ];
else
    cutMesh.u = [ QS.a_fixed .* uvtx(:,1), uvtx(:,2) ];
end
cutMesh.ar = ar ;
cutMesh.umax = QS.a_fixed ;

% preview the pullback mesh uv coords
% plot(cutMesh.u(:, 1), cutMesh.u(:, 2), '.')
% title(['generateCurrentCutMesh(): showing cutMesh.u for t= ' num2str(tt)])
% pause(2)
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Done flattening cutMesh. Now saving.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save cutMesh
disp(['Saving cutMesh to ' cutMeshfn]) 
save(cutMeshfn, 'cutMesh', 'adIDx', 'pdIDx', 'cutP')

% Displacing mesh along normal direction
disp(['Evolving mesh along normal shift for pullback ',...
    'images: shift=' num2str(QS.normalShift)])
cutMesh.v = cutMesh.v + cutMesh.vn * QS.normalShift ;

disp('Assigning current cutMesh to self')
QS.currentMesh.cutMesh = cutMesh ;
disp('done with plotting & saving cut')
