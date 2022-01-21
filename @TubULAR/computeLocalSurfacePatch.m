function [subm, newpts] = computeLocalSurfacePatch(QS, pts, options)
% Make local surface parameterization enclosing pts in pullback space
% (could be used for computing angle between two cells during T1s)
% 
% Note: additional points can be found in mapped space via
% newpts = barycentricMap2d(subm.f, subm.u, subm.V2D_scaled_detg, pts)
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : optional struct with fields
%   coordSys : str specifier
%   overwrite : bool
%   buffer : float or int or numeric
%   bufferX : float or int or numeric, trumps buffer for X dim
%   bufferY : float or int or numeric, trumps buffer for Y dim
%   preview : preview
%   scaleByMetric : bool
%       rescale parameterization to have det(g) approx 1 on either:
%           - faces queried by pts (if both imW and imH are not supplied)
%           - in window com(pts) +/- [imW*0.5, imH*0.5] 
%   scaleByMetricComponents : bool
%       rescale parameterization to approximate isothermal coordinates
%   imW : numeric
%       if scaleByMetric or scaleByMetricComponents, supplying imW then 
%       determines the width in X dimension of V2D parameterization in which 
%       to query metric for fixing. If none supplied, uses pointLocation of pts
%       in the submesh to probe which faces are queried.
%       If supplied, probeMetricViaFaceLookup = false ;
%       (in units of quasi-embedding space --> ie registered/flattened
%       embedding submesh units, should be close to units of QS.spaceUnits if
%       submesh is not too strongly curved in embedding space).
%   imH : numeric
%       if scaleByMetric or scaleByMetricComponents, supplying imH then 
%       determines the width in Y dimension of V2D parameterization in which 
%       to query metric for fixing. If none supplied, uses pointLocation of pts
%       in the submesh to probe which faces are queried.
%       If supplied, probeMetricViaFaceLookup = false ;
%       (in units of quasi-embedding space --> ie registered/flattened
%       embedding submesh units, should be close to units of QS.spaceUnits if
%       submesh is not too strongly curved in embedding space).
%
% 
% Returns
% -------
% subm : struct with fields
%   f : (#faces x 3 int) face Connectivity List
%   v : (#vertices x 3 float) vertices in embedding space (must be 3D)
%   u : (#vertices x 2 float) original pullback coordinates of submesh
%   Vcut : (#vertices x 3 float) vertices in embedding space (must be 3D)
%       This is typically the same as v, in which case it is not returned
%   Fcut : (#faces x 3 int) face Connectivity List of cutMesh submesh
%       This is typically the same as f, in which case it is not returned
%   Ucut : (#vertices x 2 float) original pullback coordinates of submesh
%       This is typically the same as u, in which case it is not returned
%   V2D : #parameterized
%       as-rigid-as-possible re-parameterization of the embedding in a
%       local patch of pullback space 
%   V2D_scaled_detg : (returned if scaleByMetric==true)
%   V2D_scaled_g11g22 : (returned if scaleByMetricComponents==true)
%   connectivityCase : int indicator
%       1. submesh is within single Cover (most common, disk-like)
%            --> no problems here, all triangles are oriented
%            correctly
%       2. submesh has one edge on another cover, with some long faces,
%           but is still a topological disk, so we push those vertices 
%           reaching across the branch cut to the other side so that 
%           triangles are oriented correctly
%       3. submesh has faces on multiple covers (topological disk)
%       4. submesh spans the whole single Cover (topological annulus)
%
% See also
% --------
% visualizeSegmentationPatch(QS, options)
% visualizeDemoTracks(QS, options)
% 
% NPMitchell 2021

% Unpack options
bufferX = 0 ;
coordSys = 'spsme';
scaleByMetric = false ;
scaleByMetricComponents = true ;
bufferX = 0.1 ;
bufferY = 0.1 ;
preview = false ;

% Unpack options
if nargin < 3
    options = struct() ;
end

if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'buffer') 
    bufferX = options.buffer ; 
    bufferY = options.buffer ; 
end
if isfield(options, 'bufferX') 
    bufferX = options.bufferX ; 
end
if isfield(options, 'bufferY') 
    bufferY = options.bufferY ; 
end
if isfield(options, 'preview') 
    preview = options.preview ; 
end
if isfield(options, 'scaleByMetric')
    scaleByMetric = options.scaleByMetric ;
elseif isfield(options, 'scaledByMetric')
    % allow typo
    scaleByMetric = options.scaledByMetric ;
end
if isfield(options, 'scaleByMetricComponents')
    scaleByMetricComponents = options.scaleByMetricComponents ;
end

if scaleByMetric || scaleByMetricComponents
    % Decide on patch in which to probe metric
    if ~isfield(options, 'imW') && ~isfield(options, 'imH')
        % Use points to query max and min bounds for patch using
        % pointLocation
        probeMetricViaFaceLookup = true ;
    else
        if isfield(options, 'imW')
            imW = options.imW ;
            if ~isfield(options, 'imH')
                imH = imW ;
            end
        end
        if isfield(options, 'imH')
            imH = options.imH ;
            if ~isfield(options, 'imW')
                imW = imH ;
            end
        end
        probeMetricViaFaceLookup = false ;
    end
end

% Make sure we know what time it is
try
    assert(~isempty(QS.currentTime))
catch
    error('First set current time: QS.setTime(timepoint)')
end

% Get XY limits for regions
if strcmpi(coordSys, 'spsme')
    % Get the image and meshes
    disp(['Loading images and meshes for timepoint: ' num2str(QS.currentTime)])
    cutMesh = QS.getCurrentSPCutMeshSm() ;
    glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
    cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
    glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1)) ;
    bc = barycenter(cutMesh.u, cutMesh.f) ;
    im = imread(sprintf(QS.fullFileBase.im_sp_sme, QS.currentTime)) ;
    shiftY = size(im, 1) * 0.5 ;
    doubleCovered = true ;
else
    error('handle desired coordSys here')
end

% Get XY limits for fixed frame
Xbuffer = bufferX * size(im, 1) ;
Ybuffer = bufferY * size(im, 1) ;

% Incude all faces in [xmin,xmax], [ymin,ymax]
% Crop the patch by a bounding box
xmin = min(pts(:, 1)) - Xbuffer ;
xmax = max(pts(:, 1)) + Xbuffer ;
miny = min(pts(:, 2)) ;
maxy = max(pts(:, 2)) ;
bcXY = QS.uv2XY(im, bc, doubleCovered, 1., 1.) ;
faces2rm_initial = bcXY(:, 1) < xmin | bcXY(:, 1) > xmax |...
    (abs(bcXY(:, 2) - miny) > Ybuffer & ...
    abs(bcXY(:, 2) - maxy) > Ybuffer & ...
    abs(bcXY(:, 2) - miny - shiftY) > Ybuffer & ...
    abs(bcXY(:, 2) - maxy - shiftY) > Ybuffer & ...
    abs(bcXY(:, 2) - miny + shiftY) > Ybuffer & ...
    abs(bcXY(:, 2) - maxy + shiftY) > Ybuffer);

% Figure out which faces to remove, including ears
[ submV, submF, facesKept ] = removeMeshFaces(glueMesh.v, glueMesh.f,...
    faces2rm_initial ) ;
[~, ~, ~, facesKept2] = removeMeshEars(submV, submF) ;
% [newFace, newFaceIdx, newVertex, newVertexIdx] = clip_mesh_mod(submF, submV) ;
others2rm = setdiff(1:size(submF, 1), facesKept2) ;
others2rm = facesKept(others2rm) ;
faces2rm = unique([find(faces2rm_initial); others2rm]) ;
assert(length(faces2rm) == length(find(faces2rm_initial)) + length(others2rm))

% Now create ear-free submesh
[ submTV, submTF ] = removeMeshFaces(cutMesh.v, cutMesh.f,...
    faces2rm ) ;
[ submV, submF ] = removeMeshFaces(glueMesh.v, glueMesh.f,...
    faces2rm ) ;
[ submU, ~ ] = removeMeshFaces(glueMesh.u, glueMesh.f,...
    faces2rm ) ;
[ submUcut, submFcut ] = removeMeshFaces(cutMesh.u, cutMesh.f,...
    faces2rm ) ;

% There are a few cases for how this looks:
% 1. submesh is within single Cover (most common, disk-like)
%       --> no problems in sight, all triangles are oriented
%       correctly
% 2. submesh has one edge on another cover, with some long faces,
% but is still a topological disk
%   
% 3. submesh has faces on multiple covers (topological disk)
% 
% 4. submesh spans the whole single Cover (topological annulus)
bcdiff = barycenter(submUcut, submFcut) -barycenter(submU, submF) ;
assert(~any(abs(bcdiff(:, 1)) > 1e-13))
subm = struct('f', submF, 'v', submV) ; % for checking eulerCharacteristic
if ~any(abs(bcdiff(:, 2) > 1e-13 ))
    % No periodic boundary passed
    connectivityCase = 1;
elseif eulerCharacteristic(subm) == 1
    % must be either case 2 or 3
    if all(size(submUcut) == size(submU))
        % Same number of vertices, so we haven't actually crossed
        % the branch cut, just are hugging it
        connectivityCase = 2 ;

        if median(submUcut(:, 2)) > 0.5 
            submU(submU(:, 2) == 0, 2) = 1 ;
        elseif median(submUcut(:, 2)) < 0.5
            submU(submU(:, 2) == 1, 2) = 0 ;
        else
            error('handle how to cut the submesh here')
        end
    else
        connectivityCase = 3 ;
    end
elseif eulerCharacteristic(subm) == 0
    connectivityCase = 4; 
end

% Package into struct
subm = struct('f', submF, 'v', submV, 'u', submU, ...
    'connectivityCase', connectivityCase) ;
if connectivityCase > 1
    subm.Vcut = submTV ;
    subm.Fcut = submTF ;
    subm.Ucut = submUcut ;
end


% Check the submesh
if preview
    trisurf(triangulation(subm.f, subm.v))
end

% SURFACE PARAMETERIZATION
param = struct();
param.borderType = 2; % FIXED = 1, FREE = 2
param.fixedType = 1; % ARC_LENGTH = 1, UNIFORM = 2
param.fixedShape = 1; % CIRCLE = 1, SQUARE = 2
% FIXED: BARYCENTRIC = 1, AUTHALIC = 2, CONFORMAL = 3, MEAN = 4
% FREE: LSCM = 1, ARAP = 2
param.paramMethod =  2;
% param.corners = [];
% param.fixedPoints = [];
disp('Generating surface parameterization... ');
[ V2D, outFaces ] = surface_parameterization( subm.f, subm.v, param );

% Check for spuriously added points -- sometimes a point at the
% origin is added at the end it seems(?). Not sure why
if size(V2D, 1) ~= size(subm.v, 1)
    error('Points added/removed. Handle this case here. This is abnormal')
    if size(V2D, 1) == size(subm.v, 1) + 1
        assert(all(V2D( size(subm.v, 1) + 1, :) == 0))
        [ outFaces, V2D] = ...
            remove_vertex_from_mesh(outFaces, V2D,  size(subm.v, 1) + 1)
    end
end

if param.borderType == 1
    if param.fixedShape == 1
        V2D = 2 .* ( V2D - 0.5 ); % Map disk to unit disk
    end
end

% Register to orginal mesh in space
% Note: there could be some issues with the periodic BC felt by
% submU's influence here!! todo
[phi, u, v] = rigidParameterizationAlignment( ...
    submF, submV, submU, submF, submV, V2D ) ;
rotm = [cos(phi), sin(phi); -sin(phi), cos(phi)] ;
V2Dr = (V2D + [u, v]) * rotm  ;
subm.V2D = V2Dr ;



%% Account for possible total dilation/contraction near the
% tracked cells -- this can be significant for large patches
if scaleByMetric || scaleByMetricComponents
    % Compute metric in a small window near the origin
    [gg, ~] = constructFundamentalForms(submF, submV, V2Dr) ;
    bcV2D = barycenter(V2Dr, submF) ;
    
    % Either query the metrics whose faces are sampled by pts or are within
    % imW, imH of Ucom (com of supplied pts in registered pullback space)
    if probeMetricViaFaceLookup
        ptsU = QS.XY2uv(im, pts, doubleCovered, 1., 1.) ;
        cellu = barycentricMap2d(submF, submU, V2Dr, ptsU) ;
        cellFaces = unique(pointLocation(...
            triangulation(submF, V2Dr), cellu)) ;
        gF = zeros(length(cellFaces), 2, 2) ;
        for cfid = 1:length(cellFaces)
            gF(cfid, :, :) = gg{cellFaces(cfid)} ;
        end
        inBox = cellFaces ;
    else
        Ucom = QS.XY2uv(im, mean(pts), doubleCovered, 1., 1.) ;
        comV2D0 = barycentricMap2d(submF, submU, V2Dr, Ucom) ;
        inBox = find(abs(bcV2D(:, 1) - comV2D0(1)) <imW*0.5 &  ...
            abs(bcV2D(:, 2) - comV2D0(2)) <imH*0.5);
        gF = zeros(size(inBox, 1), 2, 2) ;
        for cfid = 1:size(inBox, 1)
            gF(cfid, :, :) = gg{inBox(cfid)} ;
        end
    end
    
    % Check that we are looking in the box
    if preview
        figure(3);
        triplot(triangulation(submF, V2Dr))
        hold on; 
        plot(bcV2D(inBox, 1), bcV2D(inBox, 2), 'o')
        if scaleByMetricComponents
            title('Triangles whose metric we force to be approx isothermal')
        else
            title('Triangles whose metric we scale to be approx det=1')
        end
        waitfor(gcf)
    end


    % Scale by det(g) -- either mean dilation 
    if scaleByMetricComponents
        % Scale each dim separately
        dilX = mean(gF(:, 1, 1)) ;
        dilY = mean(gF(:, 2, 2)) ;
        V2Dr(:, 1) = V2Dr(:, 1) * sqrt(dilX) ;
        V2Dr(:, 2) = V2Dr(:, 2) * sqrt(dilY) ;
        subm.V2Dr_scaled_g11g22 = V2Dr ;

        % plot preview of the scaling and its variation
        if preview
            close all
            subplot(1, 2, 1)
            trisurf(submF(inBox, :), V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                gF(:, 1, 1), 'edgecolor', 'none')
            view(2) ; grid off ;
            axis equal
            cb = colorbar ;
            ylabel(cb, '$g_{11}$', 'interpreter', 'latex')

            subplot(1, 2, 2)
            trisurf(submF(inBox, :), V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                gF(:, 2, 2), 'edgecolor', 'none')
            view(2) ; grid off ;
            axis equal
            cb = colorbar ;
            ylabel(cb, '$g_{22}$', 'interpreter', 'latex')
            sgtitle('scaling in ARAP mapping', 'interpreter', 'latex')
            set(gcf, 'color', 'white')
            pause(5) ;
            clf
        end
    else
        cellg = [mean(gF(:, 1, 1)), mean(gF(:, 1, 2));...
            mean(gF(:, 2, 1)), mean(gF(:, 2, 2))] ;
        celldetg = sqrt(det(cellg)) ;
        V2Dr = V2Dr * sqrt(celldetg) ;
        subm.V2D_scaled_detg = V2Dr ;

        % plot preview of the scaling and its variation
        if preview
            clf
            dets = zeros(size(submF, 1), 1) ;
            for cfid = 1:size(submF, 1)
                dets(cfid) = sqrt(det(gg{cfid})) ;
            end
            trisurf(submF, V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                dets, 'edgecolor', 'none')
            view(2)
            axis equal
            cb = colorbar ;
            ylabel(cb, '$\sqrt(\det g)$', 'interpreter', 'latex')
            title('scaling in ARAP mapping', 'interpreter', 'latex')
            pause(5) ;
            clf
        end
    end
end

% Supply new points (advected onto V2Dr)
if nargout > 1
    newpts = barycentricMap2d(submF, submU, V2Dr, ptsU) ;
end