function [mesh, out_of_frame] = medianFilterMongeMesh(mesh, gridOptions, planeOptions, preview)
% For a mesh that is protruding from a surface parametrizable by z(x,y),
% take median z for each bin in binned (x,y) grid and define new surface.
% Also allows for local max or local min instead of local median (see
% takeMedianMaxMin in gridOptions)
%
% Parameters
% ----------
% gridOptions : struct with fields
%   uminmax, vminmax, nU, nV
%   axisOrder : 'xyz', 'yxz', etc
%   resampleMethod : 'median', 'min', 'max' 'mean' 'interpolate' (default='interpolate')
%       method for finding z values in resampled grid
%   erosionRadius : int (default=0) (used to be=5)
%       after resampling, take local min by eroding with strel of this 
%       radius
%   dilationRadius : int (default=0)
%       after resampling and erosion, take local max by dilating with strel
%       of this radius
%   gaussSigma : float (default=0)
%       after resampling, erosion, and dilation, smooth z(x,y)
%   zOffset : float (default=0)
%       offset in z to artificially subtract
% planeOptions : struct with fields
%   maxIsPlane, 
%   thresPlane, 
%   Z0
%
% Returns
% -------
%
%
% NPMitchell 2021

%% Default Options

% Grid Options
axisOrder = 'xyz' ;
resampleMethod = 'interpolate'; % could be 'min', 'max', 'median' 'mean' 'interpolate'
zOffset = 0 ;
% number of x bins
nU = 150 ;
uminmax = [min(mesh.v(:, 1)), max(mesh.v(:, 1))] ;
% number of y bins
nV = 150 ;
vminmax = [min(mesh.v(:, 2)), max(mesh.v(:, 2))] ;
szX = uminmax(2) - uminmax(1) ;
szY = vminmax(2) - vminmax(1) ;
erosionRadius = 0 ;
dilationRadius = 0 ;
gaussSigma = 0 ;

% Plane options
maxIsPlane = true ;  % Monge form is protruding from plane 
thresPlane = eps ; % used to be 0.1 
% may also supply Z0 (or 'z0')


%% Supply Options
if nargin > 1 && isfield(gridOptions, 'uminmax')
    uminmax = gridOptions.uminmax ; 
end
if nargin > 1 && isfield(gridOptions, 'vminmax')
    vminmax = gridOptions.vminmax ; 
end
if nargin > 1 && isfield(gridOptions, 'nU')
    nU = gridOptions.nU ; 
end
if nargin > 1 && isfield(gridOptions, 'nV')
    nV = gridOptions.nV ; 
end
if nargin > 1 && isfield(gridOptions, 'szX')
    szX = gridOptions.szX ; 
end
if nargin > 1 && isfield(gridOptions, 'szY')
    szY = gridOptions.szY ; 
end
if nargin > 1 && isfield(gridOptions, 'axisOrder')
    axisOrder = gridOptions.axisOrder ; 
end
if nargin > 1 && isfield(gridOptions, 'zOffset')
    zOffset = gridOptions.zOffset ; 
end
if nargin > 1 && isfield(gridOptions, 'erosionRadius')
    erosionRadius = gridOptions.erosionRadius ; 
end
if nargin > 1 && isfield(gridOptions, 'dilationRadius')
    dilationRadius = gridOptions.dilationRadius ; 
end

% Plane options
if nargin > 2 && isfield(planeOptions, 'maxIsPlane')
    maxIsPlane = planeOptions.maxIsPlane ;
end
if nargin > 3 && isfield(planeOptions, 'thresPlane')
    thresPlane = planeOptions.thresPlane ;
end

% Get zdim
axisOrder = lower(axisOrder) ;
[~, zdim] = find(axisOrder == 'z');

% Determine z0, the reference plane's z value (if a reference plane exists) 
if nargin > 3 && isfield(planeOptions, 'Z0')
    Z0 = planeOptions.Z0 ;
elseif nargin > 3 && isfield(planeOptions, 'z0')
    Z0 = planeOptions.z0 ;
else
    if maxIsPlane
        Z0 = max(mesh.v(:, zdim)) ;
    else
        Z0 = min(mesh.v(:, zdim)) ;
    end
end

%% Remove floor of mesh if present (if closed sphere, for ex, remove "back")
% Floor is max(x)
if maxIsPlane
    % used to be max(mesh.v(:, zdim)) instead of z0
    plane = find(mesh.v(:, zdim) > Z0 - thresPlane) ;
else
    plane = find(mesh.v(:, zdim) < Z0 + thresPlane) ;    
end

%% Subsample the surface 
[face, vertex] = remove_vertex_from_mesh(mesh.f, mesh.v, plane) ;
[ surfF, surfV, oldVertexIDx ] = ...
     remove_isolated_mesh_components( face, vertex, 1000) ;

% check it
if preview
    trisurf(triangulation(surfF, surfV), surfV(:, 3), 'edgecolor', 'none')
    axis equal;
    view(2)
end

% surfV = mesh.v(surfId, :) ;
if contains(axisOrder, 'zyx')
    uvz = [surfV(:, 3), surfV(:, 2), surfV(:, 1)] ;
elseif contains(axisOrder, 'yxz')
    uvz = [surfV(:, 2), surfV(:, 1), surfV(:, 3)] ;
elseif contains(axisOrder, 'xyz')
    uvz = [surfV(:, 1), surfV(:, 2), surfV(:, 3)] ;
    
    if preview
        scatter3(uvz(:, 1), uvz(:, 2), uvz(:, 3), 5, uvz(:, 3))
        axis equal
        xlabel('x'); ylabel('y'); zlabel('z') ;
    end
else
    error('handle axisOrder here')
end

% subsample by binning and finding medians or interpolating
uvec = linspace(uminmax(1), uminmax(2), nU) ;
vvec = linspace(vminmax(1), vminmax(2), nV) ;
[xx, yy] = meshgrid(uvec, vvec) ;
xx = xx';
yy = yy';
if strcmpi(resampleMethod, 'interpolate')
    
    % intp = scatteredInterpolant(uvz(:,1), uvz(:, 2), uvz(:, 3), 'natural', 'none') ;
    % zvec = intp(ugrid(:), vgrid(:)) ;
    % med2 = reshape(zvec, [nU, nV]) ;
    
    zvec = interpolate2Dpts_3Dmesh(surfF, surfV(:, [1,2]), surfV, [xx(:), yy(:)]) ;
    medz = zvec(:, 3) ;
    medz = reshape(medz, [nU, nV]) ;
else
    % Take local mean/median/min/max
    [zmeans, counts, zs, xidx, yidx] = ...
        binData2dGrid(uvz, uminmax, vminmax, nU, nV, false) ;
    medz = zeros(nU, nV) ;
    for pp = 1:size(zs, 1)
        for qq = 1:size(zs, 2)
            if isempty(zs{pp, qq})
                medz(pp, qq) = NaN ;
            else
                if strcmpi(resampleMethod, 'median')
                    medz(pp, qq) = median(zs{pp, qq}) - zOffset ;
                elseif strcmpi(resampleMethod, 'max')
                    medz(pp, qq) = max(zs{pp, qq}) - zOffset ;
                elseif strcmpi(resampleMethod, 'min')
                    medz(pp, qq) = min(zs{pp, qq}) - zOffset ;
                end
            end
        end
    end
    
    % Create x and y array
    [yy, xx] = meshgrid(linspace(uminmax(1), uminmax(2), nU), ...
        linspace(vminmax(1), vminmax(2), nV)) ;
    xx = xx';
    yy = yy';
end
% NOTE: This should look correct but have NaNs
if preview
    imagesc( vvec, uvec, medz)
    axis equal ;
    title('resampled surface')
    pause(1)
end    

out_of_frame = find(isnan(medz)) ;
medz(out_of_frame) = Z0 ;

% local minimum filter (lower z is further from midsagittal plane)
if erosionRadius > 0
    medz = imerode(medz, true(erosionRadius)) ;
end

% local minimum filter (lower z is further from midsagittal plane)
if dilationRadius > 0
    medz = imdilate(medz, true(dilationRadius)) ;
end

% local smoothing filter 
if gaussSigma > 0
    medz = imgaussfilt(medz, gaussSigma) ;
end

% NOTE: This should look correct -- now smoothed
if preview
    imagesc( vvec, uvec, medz)
    axis equal ;
    title('resampled surface, smoothed')
    pause(1)
end    

% Create ring of zeros around shape
% se = strel('disk', 2);
% keep = find(imdilate(~isnan(zmeans), se)) ;
% assert(all(size(medz) == size(zmeans)))
 
xr = xx(:) ;
yr = yy(:) ;
zr = medz(:) ;
% scatter3(xr, yr, zr, 4, zr); view(2); axis equal

% xr = xr(keep) ;
% yr = yr(keep) ;
% zr = zr(keep) ;
% xr = [xr; 0; sz1; sz1; 0] ;
% yr = [yr; 0; 0; sz2; sz2] ;
% zr = [zr; Z0; Z0; Z0; Z0] ;

% Inspect cross-section (this is turned off for now)
if preview && false
    % Make this number larger to sample more of the nearby mesh
    width = 4 ;
    leaves = 1:50:szY ;
    % Show the cross-section
    for leaf = leaves
        inds = find(abs(yy - leaf) < width) ;
        clf
        if strcmpi(erase(axisOrder, 'c'), 'xyz')
            im = squeeze(IV{2}(:, leaf, :))' ;
        elseif strcmpi(erase(axisOrder, 'c'), 'yxz')
            im = squeeze(IV{2}(leaf, :, :))' ;
        end
        imshow(mat2gray(im, [0, double(rms1d(im(:)))]))
        hold on; 
        if any(inds)
            hold on;
            plot(xx(inds), medz(inds), 'c.')
        end
        pause(0.1)
    end
end

% Triangulate the result
faces = delaunay(xr, yr) ;
mesh = struct() ;
mesh.f = faces ;
mesh.v = [xr, yr, zr] ;

