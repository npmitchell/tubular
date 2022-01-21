function [xyzlim_raw, xyzlim, xyzlim_um, xyzlim_um_buff] = ...
    measureXYZLims(QS, options)
% Redo calculation of XYZ limits
%
% Parameters
% ----------
% 
%
% NPMitchell 2020

% Unpack options
overwrite = false ;
if nargin > 1
    if isfield(options, 'overwrite')
        overwrite = true ;
    end
end

%% Unpack QS
timePoints = QS.xp.fileMeta.timePoints ;

% Data file names
xyzlimname_raw = QS.fileName.xyzlim_raw ;
xyzlimname_pix = QS.fileName.xyzlim_pix ;
xyzlimname_um = QS.fileName.xyzlim_um ;
xyzlimname_um_buff = QS.fileName.xyzlim_um_buff ;

% on_disk = exist(xyzlimname_raw, 'file') && ...
%      exist(xyzlimname_pix, 'file') &&  exist(xyzlimname_um, 'file') && ...
%       exist(xyzlimname_um_buff, 'file') ; 
   
disp('Extracting xyzlimits for raw meshes...')
disp([' ... since ' xyzlimname_raw ' does not exist'])
for tidx = 1:length(timePoints)
    tt = timePoints(tidx) ;
    if mod(tidx, 10) == 1
        disp(['tt = ', num2str(tt)])
    end
    
    % Get the timestamp string from the name of the mesh
    mesh = read_ply_mod(sprintf(QS.fullFileBase.mesh, tt)) ;
    
    % Could load APDV from disk
    % meshAPDVfn = sprintf(QS.fullFileBase.alignedMesh, tt) ; 
    % if exist(meshAPDVfn, 'file')
    %     meshAPDV = read_ply_mod(sprintf(QS.fullFileBase.mesh, tt)) ;
    % else
    % end
    
    % Option 2: rotate/scale/translate in-place
    meshAPDV = struct() ;
    meshAPDV.f = mesh.f ;
    meshAPDV.v = QS.xyz2APDV(mesh.v) ;
    
    % Check 
    % trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), ...
    %     mesh.v(:, 3), mesh.v(:, 3), 'edgecolor', 'none')
    % trisurf(meshAPDV.f, meshAPDV.v(:, 1), meshAPDV.v(:, 2), ...
    %     meshAPDV.v(:, 3), meshAPDV.v(:, 3), 'edgecolor', 'none')
    
    % Get limits
    minx = min(mesh.v) ;
    maxx = max(mesh.v) ;
    if tt == timePoints(1)
        xmin = minx(1) ;
        ymin = minx(2) ;
        zmin = minx(3) ;
        xmax = maxx(1) ;
        ymax = maxx(2) ;
        zmax = maxx(3) ;
    else
        xmin= min(xmin, minx(1)) ;
        ymin = min(ymin, minx(2)) ;
        zmin = min(zmin, minx(3)) ;
        xmax = max(xmax, maxx(1)) ;
        ymax = max(ymax, maxx(2)) ;
        zmax = max(zmax, maxx(3)) ;
    end 
    
    % Get APDV limits
    minxAPDV = min(meshAPDV.v) ;
    maxxAPDV = max(meshAPDV.v) ;
    if tt == timePoints(1)
        xminrs = minxAPDV(1) ;
        yminrs = minxAPDV(2) ;
        zminrs = minxAPDV(3) ;
        xmaxrs = maxxAPDV(1) ;
        ymaxrs = maxxAPDV(2) ;
        zmaxrs = maxxAPDV(3) ;
    else
        xminrs = min(xminrs, minxAPDV(1)) ;
        yminrs = min(yminrs, minxAPDV(2)) ;
        zminrs = min(zminrs, minxAPDV(3)) ;
        xmaxrs = max(xmaxrs, maxxAPDV(1)) ;
        ymaxrs = max(ymaxrs, maxxAPDV(2)) ;
        zmaxrs = max(zmaxrs, maxxAPDV(3)) ;
    end 
    
end

% Todo: save raw xyzlim in full resolution pixels but not rotated/scaled
% Save xyzlim_raw
if overwrite || ~exist(xyzlimname_raw, 'file')
    disp('Saving rot/trans mesh xyzlimits for plotting')
    header = 'xyzlimits for raw meshes in units of full resolution pixels' ;
    xyzlim_raw = [xmin, xmax; ymin, ymax; zmin, zmax] ;
    write_txt_with_header(xyzlimname_raw, xyzlim_raw, header) ;
else
    xyzlim_raw = [xmin, xmax; ymin, ymax; zmin, zmax] ;
end

% Save xyzlimits 
resolution = QS.APDV.resolution ;
if overwrite || ~exist(xyzlimname_pix, 'file')
    disp('Saving rot/trans mesh xyzlimits for plotting')
    header = 'xyzlimits for rotated translated meshes in units of full resolution pixels' ;
    xyzlim = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] / resolution;
    write_txt_with_header(xyzlimname_pix, xyzlim, header) ;
else
    xyzlim = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] / resolution;
end

% Save xyzlimits in um
if overwrite || ~exist(xyzlimname_um, 'file')
    disp('Saving rot/trans mesh xyzlimits for plotting, in microns')
    header = 'xyzlimits for rotated translated meshes in microns' ;
    xyzlim_um = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] ;
    write_txt_with_header(xyzlimname_um, xyzlim_um, header) ;
else
    xyzlim_um = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] ;
end

% Save buffered xyzlimits in um
if overwrite || ~exist(xyzlimname_um_buff, 'file')
    disp('Saving rot/trans mesh xyzlimits for plotting, in microns')
    header = 'xyzlimits for rotated translated meshes in microns, with padding (buffered)' ;
    xyzlim_um = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] ;
    xyzlim_um_buff = xyzlim_um + 2 * abs(QS.normalShift) * resolution * [-1, 1] ;
    write_txt_with_header(xyzlimname_um_buff, xyzlim_um_buff, header) ;
else
    xyzlim_um_buff = xyzlim_um + 2 * abs(QS.normalShift) * resolution * [-1, 1] ;
end