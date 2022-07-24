function plotPathlineStrain(tubi, options)
% plotPathlineStrain(QS, options)
%   Plot the strain (from piv pathlines) along pathlines as kymographs and 
%   correlation plots. These are computed via finite differencing of 
%   pathline mesh metrics. That is, vertices are advected along piv flow 
%   and projected into 3D from pullback space to embedding, then the
%   metrics g_{ij} of those meshes are measured and differenced. 
% 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs         : bool
%   plot_kymographs_cumsum  : bool
%   plot_correlations       : bool 
%   plot_fold_strainRate    : bool
%   plot_lobe_strainRate    : bool
%   plot_fold_strain        : bool
%   plot_lobe_strain        : bool
% 
% Returns
% -------
% <none>
%
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_strain3d = false ;
maxWFrac = 0.03 ;  % Maximum halfwidth of feature/fold regions as fraction of ap length
t0 = tubi.t0set() ;
t0Pathline = t0 ;
viewAngles = [-20,25] ;
clim_tre = 2 ;
clim_dev = 5 ;

strain_trace_label = '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon]$';
strain_trace_label_norm = '$||\mathrm{Tr} [\varepsilon]||$';
strain_deviator_label = ...
    '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$' ;
strain_deviator_label_short = ...
    '$||$Dev$\left[\varepsilon\right]||$' ;
% strain_theta_label = '$\theta_{\mathrm{Dev}[{\varepsilon}]}$'; 

%% Parameter options
lambda = tubi.smoothing.lambda ;
lambda_mesh = tubi.smoothing.lambda_mesh ;
nmodes = tubi.smoothing.nmodes ;
zwidth = tubi.smoothing.zwidth ;
climit = 0.2 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
%% parameter options
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
end

% smoothing parameter options
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'nmodes')
    nmodes = options.nmodes ;
end
if isfield(options, 'zwidth')
    zwidth = options.zwidth ;
end

% plotting and sampling
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'climitWide')
    climitWide = options.climitWide ;
else
    climitWide = climit * 3; 
end
if isfield(options, 'maxWFrac')
    maxWFrac = options.maxWFrac ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_strain3d')
    plot_strain3d = options.plot_strain3d ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
else
    error("Could not parse samplingResolution: set to '1x' ")
end

%% Unpack QS
tubi.getXYZLims ;
xyzlim = tubi.plotting.xyzlim_um ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
bbr256 = blueblackred(256) ;
% clf
% set(gcf, 'visible', 'off')
% imagesc([-1, 0, 1; -1, 0, 1])
% caxis([0, 1])
% pos256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;
bluecolor = tubi.plotting.colors(1, :) ;
orangecolor = tubi.plotting.colors(2, :) ;
yellowcolor = tubi.plotting.colors(3, :) ;
purplecolor = tubi.plotting.colors(4, :) ;
greencolor = tubi.plotting.colors(5, :) ;
graycolor = tubi.plotting.colors(8, :) ;
browncolor = tubi.plotting.colors(9, :) ;

%% Choose colors
trecolor = yellowcolor ;
Hposcolor = greencolor ;
Hnegcolor = purplecolor ;
Hsz = 3 ;  % size of scatter markers for mean curvature

%% load from QS
if doubleResolution
    nU = tubi.nU * 2 - 1 ;
    nV = tubi.nV * 2 - 1 ;
else
    nU = tubi.nU ;
    nV = tubi.nV ;    
end

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = tubi.xp.fileMeta.timePoints(1:end-1) - t0 ;

% Output directory is inside piv/pathlines/t0_XXXX dir
mKPDir = sprintf(tubi.dir.pathlines.strain, t0Pathline) ;
datdir = mKPDir ;
% Data for kinematics on meshes (defined on vertices) [not needed here]
% mdatdir = fullfile(mKDir, 'measurements') ;

% Unit definitions for axis labels
unitstr = [ '[unitless]' ];
tidx0 = tubi.xp.tIdx(tubi.t0set()) ;

%% Plot strains in 3d
if plot_strain3d
    pb = tubi.getPullbackPathlines([], 'vertexPathlines3d') ;
    refMesh = pb.refMesh ;
    % Set the rotated scaled vertices as field "v" for
    % inducedStrainPeriodicMesh() call
    refMesh.v = refMesh.vrs ;
    fa = doublearea(refMesh.vrs, refMesh.f) * 0.5 ;
    
    for tidx = tidx0 + 90 ; %1 :length(tps)
        tp = tubi.xp.fileMeta.timePoints(tidx) ;
        tubi.setTime(tp) ;
        
        outfn = fullfile(sprintf(tubi.dir.pathlines.strain_images, t0Pathline),...
            sprintf('smoothing_example_strain_%06d.png', tp)) ;
        
        strain = tubi.getCurrentPathlineStrain(tubi.t0, 'strain') ;
        % tre = strain.strain.tre(strain.strain.faceIDs) ;
        % dev = strain.strain.dev(strain.strain.faceIDs) ;
        tre = strain.strain.tre_sm ;
        fra = strain.strain.fractionalAreaChange_sm ;
        dev = strain.strain.dev_sm ;

        % positions of the pathline vertices
        xx = pb.vertices3d.vXrs(tidx, :, :) ;
        yy = pb.vertices3d.vYrs(tidx, :, :) ;
        zz = pb.vertices3d.vZrs(tidx, :, :) ;
        v3d = [ xx(:), yy(:), zz(:) ] ;
        mesh = struct() ;
        mesh.f = refMesh.f ;
        mesh.v = v3d ;
        mesh.u = refMesh.u ;
        mesh.pathPairs = refMesh.pathPairs ;
        mesh.nU = tubi.nU ;
        mesh.nV = tubi.nV ;
        
        [strainTensor, tre2, dev2, theta, outputStruct] = ...
            inducedStrainPeriodicMesh(mesh, refMesh, options) ;
        tre2 = tre2(strain.strain.faceIDs) ;
        dev2 = dev2(strain.strain.faceIDs) ;
        
        
        ax1 = subtightplot(2, 2, 1) ;
        % tre2 = (doublearea(mesh.v, mesh.f) * 0.5 - fa) ./ fa ;
        trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', tre2, 'edgecolor', 'none')
        colorbar()
        axis equal
        view(viewAngles)
        caxis(clim_tre*[-1,1])
        colormap(ax1, brewermap(256, '*RdBu'))
        title('tre')
        grid off
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        
        ax2 = subtightplot(2, 2, 2) ;
        trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', fra, 'edgecolor', 'none')
        colorbar()
        axis equal
        view(viewAngles)
        caxis(clim_tre*[-1,1])
        colormap(ax2, brewermap(256, '*RdBu'))
        title('\deltaA/A_0')
        grid off
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        
        ax3 = subtightplot(2, 2, 3) ;
        trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', dev2, 'edgecolor', 'none')
        colorbar()
        axis equal
        view(viewAngles)
        caxis([0,clim_dev])
        colormap(ax3, batlowk)
        grid off
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        
        ax4 = subtightplot(2, 2, 4) ;
        trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', dev, 'edgecolor', 'none')
        colorbar()
        axis equal
        view(viewAngles)
        caxis([0,clim_dev])
        colormap(ax4, batlowk)
        grid off
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        
        set(gcf, 'color', 'w')
        export_fig(outfn, '-nocrop', '-r600')
        
    end
end


%%
if plot_kymographs

    %% Compute or load all timepoints
    apSKymoFn = fullfile(datdir, 'apKymographLagrangianMetric.mat') ;
    lSKymoFn = fullfile(datdir, 'leftKymographLagrangianMetric.mat') ;
    rSKymoFn = fullfile(datdir, 'rightKymographLagrangianMetric.mat') ;
    dSKymoFn = fullfile(datdir, 'dorsalKymographLagrangianMetric.mat') ;
    vSKymoFn = fullfile(datdir, 'ventralKymographLagrangianMetric.mat') ;
    files_exist = exist(apSKymoFn, 'file') && ...
        exist(lSKymoFn, 'file') && exist(rSKymoFn, 'file') && ...
        exist(dSKymoFn, 'file') && exist(vSKymoFn, 'file') ;
    if files_exist
        load(apSKymoFn, 'str_apM', 'sdv_apM', 'sth_apM')
        load(lSKymoFn, 'str_lM', 'sdv_lM', 'sth_lM')
        load(rSKymoFn, 'str_rM', 'sdv_rM', 'sth_rM')
        load(dSKymoFn, 'str_dM', 'sdv_dM', 'sth_dM')
        load(vSKymoFn, 'str_vM', 'sdv_vM', 'sth_vM')
    else
        error('First run tubi.measurePathlineStrain()')
    end


    %% Store kymograph data in cell arrays
    trsK = {0.5*str_apM, 0.5*str_lM, 0.5*str_rM, 0.5*str_dM, 0.5*str_vM} ;
    dvsK = {sdv_apM, sdv_lM, sdv_rM, sdv_dM, sdv_vM} ;
    thsK = {sth_apM, sth_lM, sth_rM, sth_dM, sth_vM} ;

    %% Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKPDir, 'images', 'avgDV') ;
    lDir = fullfile(mKPDir, 'images', 'avgLeft') ;
    rDir = fullfile(mKPDir, 'images', 'avgRight') ;
    dDir = fullfile(mKPDir, 'images', 'avgDorsal') ;
    vDir = fullfile(mKPDir, 'images', 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;

    %% Now plot different measured quantities as kymographs
    if plot_kymographs
        titleadd = {': circumferentially averaged', ...
            ': left side', ': right side', ': dorsal side', ': ventral side'} ;

        for qq = 1:length(outdirs)
            % Prep the output directory for this averaging
            odir = outdirs{qq} ;
            if ~exist(odir, 'dir')
                mkdir(odir)
            end

            % Unpack what to plot (averaged kymographs, vary averaging region)
            trK = trsK{qq} ;
            dvK = dvsK{qq} ;
            thK = thsK{qq} ;

            titles = {['dilation, ' strain_trace_label],...
                ['shear, ', strain_deviator_label]} ;
            labels = {[strain_trace_label ' ' unitstr], ...
                [strain_deviator_label ' ' unitstr]} ;
            names = {'strain_dilation', 'strain_deviator'} ;

            %% Plot strainRate DV-averaged Lagrangian pathline kymographs 
            % Check if images already exist on disk

            % Consider both wide color limit and narrow
            zoomstr = {'_wide', '', '_zoom'} ;
            climits = {climit, 0.5*climit, 0.25*climit} ;

            for pp = 1:length(climits)
                %% Plot STRAIN traceful DV-averaged pathline kymograph
                % Check if images already exist on disk
                fn = fullfile(odir, [ names{1} zoomstr{pp} '.png']) ;
                if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite
                    close all
                    set(gcf, 'visible', 'off')
                    imagesc((1:nU)/nU, tps, trK)
                    caxis([-climits{pp}, climits{pp}])
                    colormap(bbr256)


                    % Titles 
                    title([titles{1}, titleadd{qq}], 'Interpreter', 'Latex')
                    ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex')
                    xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                    cb = colorbar() ;

                    % title and save
                    ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
                    disp(['saving ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')   

                    clf
                end

                %% DEVIATOR -- strain
                fn = fullfile(odir, [ names{2} zoomstr{pp} '.png']) ;
                if ~exist(fn, 'file') || overwrite
                    close all
                    set(gcf, 'visible', 'off')
                    % Map intensity from dev and color from the theta
                    indx = max(1, round(mod(2*thK(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
                    colors = pm256(indx, :) ;
                    devKclipped = min(dvK / climits{pp}, 1) ;
                    colorsM = devKclipped(:) .* colors ;
                    colorsM = reshape(colorsM, [size(dvK, 1), size(dvK, 2), 3]) ;
                    imagesc((1:nU)/nU, tps, colorsM)
                    caxis([0, climits{pp}])

                    % Add folds to plot
                    hold on;

                    % Colorbar and phasewheel
                    colormap(gca, phasemap)
                    phasebar('colormap', phasemap, ...
                        'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic') ;
                    ax = gca ;
                    get(gca, 'position')
                    cb = colorbar('location', 'eastOutside') ;
                    drawnow
                    axpos = get(ax, 'position') ;
                    cbpos = get(cb, 'position') ;
                    set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
                    set(ax, 'position', axpos) 
                    hold on;
                    caxis([0, climits{pp}])
                    colormap(gca, gray)

                    % title and save
                    title([titles{2}, titleadd{qq}], 'Interpreter', 'Latex')
                    ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex')
                    xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                    ylabel(cb, labels{2}, 'Interpreter', 'Latex')  

                    % title and save
                    ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                    disp(['saving ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')   
                    clf
                end
            end
        end
    end
end

disp('done')



