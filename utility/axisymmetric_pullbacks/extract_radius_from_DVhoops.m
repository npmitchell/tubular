%% Extract Pullback Radius from DV Hoops Pipeline =========================
% Measure radius of the pullback images. Here, use
% 'axisymmetric' pullbacks of the Drosophila midgut
% NOTE: Unlike in Generate_Axisymmetric_Pullbacks_Orbifold, this code uses
% finite vertices rather than generating face-spanning curves in DV.
% NOTE: Use Generate_Axisymmetric_Pullbacks_Orbifold.m instead!
%
% Returns
% -------
% Saves segID for each meshpoint in cutMesh
%
% Execute from the projectDir, where the data is. 
% By NPMitchell 2019, outdated -- now use instead:
% Generate_Axisymmetric_Pullbacks_Orbifold.m 
%
% This was used before it became clear that the original centerline from
% fast marching is useful only as a step towards the avgpts centerline.
% This digitizes the pullbacks of the original mesh (not the resampled DV 
% meshes) and plots the radius of those, using tricks to find the relevant
% portion of the original centerline.
%==========================================================================

clear; close all; clc;
cd /mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/

%% Options

%% Parameters
overwrite = false ;
save_ims = true ;
normal_shift = 10 ;
a_fixed = 2 ;
patch_width = 30 ;
preview = false ;
washout2d = 0.5 ;
washout3d = 0.5 ;
colorwheel_position = [.8 .01 .15 .15] ;
meshorder = 'zyx' ;  % ordering of axes in loaded mesh wrt iLastik output
axorder = [2, 1, 3] ;  % axis order from ilastik

%% Add paths
% Add some necessary code to the path (ImSAnE should also be setup!) ------
addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath(genpath('/mnt/data/code/gut_matlab/geometry'));
addpath(genpath('/mnt/data/code/gut_matlab/curve_functions'));
addpath(genpath('/mnt/data/code/gut_matlab/axisymmetric_pullbacks'));
addpath(genpath('/mnt/data/code/gut_matlab/TexturePatch'));
addpath(genpath('/mnt/data/code/gut_matlab/polarity'));
addpath(genpath('/mnt/data/code/gut_matlab/PeakFinding'));
addpath_recurse('/mnt/data/code/imsaneV1.2.3/external/') ;
addpath_recurse('/mnt/data/code/gut_matlab/plotting/') ;
% addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));

%% Define some colors
blue = [0 0.4470 0.7410] ;
orange = [0.8500 0.3250 0.0980] ;
yellow = [0.9290, 0.6940, 0.1250] ;
purple = [0.4940, 0.1840, 0.5560] ;
green = [0.4660, 0.6740, 0.1880] ;
sky = [0.3010, 0.7450, 0.9330] ;
red = [0.6350, 0.0780, 0.1840] ;
brick = [0.800000 0.250000 0.330000] ;
light_green =[0.560000 0.930000 0.560000] ;
light_gray = [0.830000 0.830000 0.830000] ;
bwr = diverging_cmap([0:0.01:1], 1, 2) ;

%% Initialize ImSAnE Project ==============================================
dataDir = pwd ;
projectDir = dataDir ;
cd( projectDir );
% A filename base template - to be used throughout this script
fileNameBase = 'Time_%06d_c1_stab';


%% Initialize Some Directory Definitions ==================================
% The top level data directory
meshDir = fullfile(dataDir, 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1') ;

% The file name base for the full meshes
fullMeshBase = fullfile( meshDir, 'mesh_apical_stab_%06d.ply' );

% The file name base for the cylinder meshes
cylinderMeshBase = fullfile( meshDir, ...
    'cylindercut/mesh_apical_stab_%06d_cylindercut.ply' );

% The file constaing the AD/PD points
dpFile = fullfile( meshDir, ...
    'cylindercut/ap_boundary_dorsalpts.h5' );

% The dataset name base for the AD points
ADBase = '/mesh_apical_stab_%06d/adorsal';
% The dataset name base for the PD points
PDBase = '/mesh_apical_stab_%06d/pdorsal';

% The folder where the pullback images will be saved
nshift = strrep(sprintf('%03d', normal_shift), '-', 'n') ;
imFolder = fullfile(meshDir, ['PullbackImages_' nshift 'step'] ) ;
imFolder_e = [imFolder '_extended' filesep] ; % debug
imFolder_r = [imFolder '_relaxed' filesep] ; % debug
% imFolder_re = [imFolder '_relaxed_extended' filesep] ;
% imFolder_es = [imFolder '_extended_shifted' filesep] ;
% pivDir = fullfile(meshDir, 'piv') ;
% polDir = fullfile(meshDir, ['polarity' filesep 'radon' ]) ;
cylcutDir = [fullfile(meshDir, 'cylindercut') filesep ];
cylcutMeshFileName = 'mesh_apical_stab_%06d_cylindercut.ply' ;
meshFileName = 'mesh_apical_stab_%06d.ply' ;
% The extensile scale factor in x for relaxing the mesh
arfn = fullfile(imFolder_r, 'ar_scalefactors.h5') ;
centerlineDir = fullfile(meshDir, 'centerline') ;
cntrsFileName = fullfile(centerlineDir, 'mesh_apical_stab_%06d_centerline_scaled_exp1p0_res1p0.txt') ;
cntrFileName = fullfile(centerlineDir, 'mesh_apical_stab_%06d_centerline_exp1p0_res1p0.txt') ;
radiusDir = fullfile(meshDir, 'radiusDVhoops') ;
radiusImDir = fullfile(radiusDir, 'images') ;
 
tomake = {radiusDir, radiusImDir} ;
for ii = 1:length(tomake)
    dir2make = tomake{ii} ;
    if ~exist( dir2make, 'dir' )
        mkdir(dir2make);
    end
end

%% Load Pullback Mesh Stack ===============================================
% Check if cutmeshes already saved
mstckfn = fullfile(meshDir, 'meshStack_orbifold.mat') ; % debug
if exist(mstckfn, 'file') 
    load(mstckfn)
else
    msg = ['Did not find ', mstckfn] ;
    msg = [msg '--> Run Generate_Axisymmetric_Pullbacks_Orbifold.m first'];
    error(msg)
end
disp('done loading meshStack')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get radius of each vertex from nearest centerline point shared by hoop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load rotation, translation, resolution
rot = dlmread(fullfile(meshDir, 'rotation_APDV.txt')) ;
xyzlim = dlmread(fullfile(meshDir, 'xyzlim.txt'), ',', 1, 0) ;
xyzlimAPDV = dlmread(fullfile(meshDir, 'xyzlim_APDV_um.txt'), ',', 1, 0) ;
trans = dlmread(fullfile(meshDir, 'translation_APDV.txt'));
resolution = dlmread(fullfile(meshDir, 'resolution.txt'), ',', 1, 0) ;
if resolution(:) == resolution(1)
    resolution = resolution(1) ;
else
    error('Have not handled case for anisotropic resolution')
end
% load meshidx and time
load(fullfile(meshDir, 'timestamps_orbifold.mat'))

eps = 1e-6 ;
Radii = cell(length(time), 1) ;
RadiiFromMean = cell(length(time), 1) ;
cIDs = cell(length(time), 1) ;
cIDList = cell(length(time), 1) ;

%% Compute the radius of points with approx equal x values from cntrline
cmap = colormap ;
pmap = phasemap ;
for ii=1:length(time)
    % Consider mesh for this timepoint
    tidx = meshidx(ii) ;
    disp(['Considering mesh ' num2str(tidx)])
    cutMesh = meshStack{tidx} ;  

    % Generate Tiled Orbifold Triangulation ------------------------------
    tileCount = [1 1];  % how many above, how many below
    v2d = cutMesh.u ;
    v3d = cutMesh.v ;
    meshf = cutMesh.f ;
    
    % [ TF, TV2D, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount );
    
    % Find the edge indices
    meshTri = triangulation( meshf, v2d );
    % The vertex IDs of vertices on the mesh boundary
    bdyIDx = meshTri.freeBoundary;
    % Consider all points on the left free boundary between y=(0, 1)
    bdLeft = bdyIDx(v2d(bdyIDx(:, 1), 1) < eps, 1) ;
    bdLeft = bdLeft(v2d(bdLeft, 2) < 1+eps & v2d(bdLeft, 2) > -eps) ;
    % Find matching endpoint on the right
    rightmost = max(v2d(:, 1));
    bdRight = bdyIDx(v2d(bdyIDx(:, 1), 1) > rightmost - eps) ;
    
    % Load centerline in raw units
    cntrfn = sprintf(cntrsFileName, time(ii)) ;
    cline = dlmread(cntrfn, ',') ;
    ss = cline(:, 1) ;
    cline = cline(:, 2:end) ;
    
    % Rotate and translate TV3D
    % cline = ((rot * cline')' + trans) * resolution  ;
    v3d = ((rot * v3d')' + trans) * resolution  ;
    % plot3(cline(:, 1), cline(:, 3), cline(:, 2), 'k-')
    % set(gcf, 'visible', 'on')
    % error('break')
    
    % Find segment of centerline to use
    % grab "front"/"start" of centerline nearest to bdLeft
    % distance from each point in bdLeft to this point in cntrline
    Adist = zeros(length(cline), 1) ;
    for kk = 1:length(cline)
        Adist(kk) = mean(vecnorm(v3d(bdLeft, :) - cline(kk, :), 2, 2)) ;
    end
    [~, acID] = min(Adist) ; 

    % grab "back"/"end" of centerline nearest to bdRight
    Pdist = zeros(length(cline), 1) ;
    for kk = 1:length(cline)
        Pdist(kk) = mean(vecnorm(v3d(bdRight, :) - cline(kk, :), 2, 2)) ;
    end
    [~, pcID] = min(Pdist) ;
    cseg = cline(acID:pcID, :) ;
    
    % Now make circumferential stripes of similar x value
    [bID, edges] = discretize(v2d(:, 1), 100) ;
    
    % check it 
    if preview
        fig = figure ;
        scatter(v2d(:, 1), v2d(:, 2), 10, bID)
        xlabel('pullback x')
        ylabel('pullback y')
        title('Hoop definitions')
        saveas(gcf, fullfile(radiusImDir, sprintf('hoop_definition_%06d.png', time(ii))))
        waitfor(fig)
        close all
    end
    
    radii = zeros(size(bID)) ;
    cids = zeros(size(bID)) ;
    avgpts = zeros(size(bID, 1), 3) ;
    radii_from_mean = zeros(size(bID)) ;
    cid_list = zeros(max(bID), 1) ;
    % Measure radius to nearest common centerline point for each bin
    for jj = 1:max(bID)
        % Find all vertices that share this bID
        hoop = find(bID == jj) ;
                
        % Compute radii for this bin
        % First find which centerline segment from which to compute radius
        rdist = zeros(length(cseg), 1) ;
        for kk = 1:length(cseg)
            rdist(kk) = mean(vecnorm(v3d(hoop, :) - cseg(kk, :), 2, 2)) ;
        end
        [~, cid] = min(rdist) ;
        
        radii(hoop) = vecnorm(v3d(hoop, :) - cseg(cid, :), 2, 2) ;
        cids(hoop) = acID + cid ;
        avgpts(jj, :) = mean(v3d(hoop, :)) ; 
        radii_from_mean(hoop) = vecnorm(v3d(hoop, :) - avgpts(jj, :), 2, 2) ;
        cid_list(jj) = acID + cid ;
    end
    
    % Save radius and mean radius
    Radii{ii} = radii ;
    RadiiFromMean{ii} = radii_from_mean;
    cIDs{ii} = cids ;
    avgPts{ii} = avgpts ;
    cIDList{ii} = cid_list ;
    
    if save_ims
        % Plot the mesh colored by radius, first in 2D
        close all
        fig = figure('visible', 'off')
        vphi = v2d(:, 2) ;
        vphi(vphi < 0) = vphi(vphi < 0) + 1 ;
        vphi(vphi > 1) = vphi(vphi > 1) - 1 ;
        scatter(v2d(:, 1), vphi, 10, radii)
        xlabel('pullback x')
        ylabel('pullback y')
        title('Radii')
        ylim([0, 1]) ;
        colorbar()
        saveas(gcf, fullfile(radiusImDir, sprintf('radii2d_%06d.png', time(ii))))
        
        % Now plot it in 3D
        close all
        fig = figure('visible', 'off')
        tmp = trisurf(meshf, v3d(:, 1), v3d(:, 2), v3d(:, 3), ...
            radii, 'edgecolor', 'none', 'FaceAlpha', 1.0) ;
        colorbar()
        hold on
        plot3(cline(:, 1), cline(:, 2), cline(:, 3), 'k-')
        plot3(cseg(:, 1), cseg(:, 2), cseg(:, 3), 'k.')
        axis equal
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        xlim(xyzlimAPDV(1, :))
        ylim(xyzlimAPDV(2, :))
        zlim(xyzlimAPDV(3, :))
        title('Radius via DV curves, Dorsal view')
        view(0, 90)
        fn = sprintf('radius_dvhoop_xyD_%06d.png', time(ii)) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radiusImDir, fn)) 
        title('Radius via DV curves, Ventral view')
        view(0, -90)
        fn = sprintf('radius_dvhoop_xyV_%06d.png', time(ii)) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radiusImDir, fn)) 
        title('Radius via DV curves, Posterior view')
        view(90, 0)
        fn = sprintf('radius_dvhoop_yzP_%06d.png', time(ii)) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radiusImDir, fn)) 
        title('Radius via DV curves, Anterior view')
        view(-90, 0)
        fn = sprintf('radius_dvhoop_yzA_%06d.png', time(ii)) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radiusImDir, fn)) 
        title('Radius via DV curves, Lateral view')
        view(0, 0)
        fn = sprintf('radius_dvhoop_xzL_%06d.png', time(ii)) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radiusImDir, fn)) 
        title('Radius via DV curves, Lateral view')
        view(0, 180)
        fn = sprintf('radius_dvhoop_xzR_%06d.png', time(ii)) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radiusImDir, fn)) 
        
        % Now plot it in 1D
        close all
        fig = figure('visible', 'off') ;
        colormap(pmap)
        scatter(ss(cids), radii, 10, vphi)
        xlabel('pathlength, $s$', 'Interpreter', 'Latex')
        ylabel('radius, $R_c$', 'Interpreter', 'Latex')
        title(['Radius from centerline t=', num2str(time(ii))])
        % phasebar %('location', 'best')
        fn = sprintf('radiusc_dvhoop_%06d.png', time(ii)) ;
        saveas(gcf, fullfile(radiusImDir, fn)) 
        close all
        
        fig = figure('visible', 'off') ;
        colormap(pmap)
        scatter(ss(cids), radii_from_mean, 10, vphi)
        xlabel('pathlength, $s$', 'Interpreter', 'Latex')
        ylabel('radius, $R_h$', 'Interpreter', 'Latex')
        title(['Radius from hoop mean t=', num2str(time(ii))])
        % phasebar %('location', 'best')
        fn = sprintf('radius_from_mean_dvhoop_%06d.png', time(ii)) ;
        saveas(gcf, fullfile(radiusImDir, fn)) 
        close all
        
        colormap(cmap)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REVIEW: 
% To understand alignment of the cutMesh, review how to transform meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% for ii=100
%     % Consider mesh for this timepoint
%     tidx = meshidx(ii) ;
%     disp(['Considering mesh ' num2str(tidx)])
%     cutMesh = meshStack{tidx} ;  
% 
%     cylfn = sprintf(fullfile(cylcutDir, cylcutMeshFileName), time(ii))  ;
%     cylMesh = read_ply_mod(cylfn) ;
%     cylface = cylMesh.f ;
%     cv = cylMesh.v ;
%     cvrs = ((rot * cv')' + trans) * resolution  ;
%     tmp = trisurf(cylface, cvrs(:, 1), cvrs(:, 2), cvrs(:, 3), ...
%         cvrs(:, 2), 'edgecolor', 'none', 'FaceAlpha', 0.3) ;
%     
%     meshfn = sprintf(fullfile(meshDir, meshFileName), time(ii))  ;
%     mesh = read_ply_mod(meshfn) ;
%     if strcmp(meshorder, 'zyx')
%         xs = mesh.v(:, 3) ;
%         ys = mesh.v(:, 2) ;
%         zs = mesh.v(:, 1) ; 
%         vn = [mesh.vn(:, 3), mesh.vn(:, 2), mesh.vn(:, 1)] ;
%     end
%     ov = [xs, ys, zs] ;
%     ovrs = ((rot * ov')' + trans) * resolution  ;
%     oface = mesh.f ;
%     % Plot the original mesh
%     hold on
%     tmp = trisurf(oface, ovrs(:, 1), ovrs(:, 2), ovrs(:, 3), ...
%         ones(size(ovrs(:, 3))), 'edgecolor', 'none', ...
%             'FaceColor', sky, 'FaceAlpha', 0.1) ;
%     
%     
%     % Generate Tiled Orbifold Triangulation ------------------------------
%     tileCount = [1 1];  % how many above, how many below
%     [ TF, TV2D, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount );
%     v3d = TV3D ;
%     size(v3d)
%     v3d = ((rot * v3d')' + trans) * resolution  ;
%     tmp = trisurf(TF, v3d(:, 1), v3d(:, 2), v3d(:, 3), ...
%         v3d(:, 2), 'edgecolor', 'none', 'FaceAlpha', 0.3) ;
%     
%     % Load centerline in raw units
%     cntrfn = sprintf(cntrFileName, time(ii)) ;
%     cline = dlmread(cntrfn, ',') ;
%     % ss = cline(:, 1) ;
%     % cline = cline(:, 2:end) ;
%     cline = ((rot * cline')' + trans) * resolution  ;
% 
%     % Plot the centerline on top
%     plot3(cline(:, 1), cline(:, 2), cline(:, 3), 'k-');
%     view(2)
%     axis equal
% end

