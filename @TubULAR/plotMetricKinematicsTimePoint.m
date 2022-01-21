function plotMetricKinematicsTimePoint(QS, tp, options)
% 
% Plot the metric Kinematics, either for all timepoints or for a single 
% timepoint. This function is called by measureMetricKinematics()
% TODO: also allow the function to be a standalone method to plot all
% timepoints if tp is not passed to options.
% 
% NPMitchell 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global operations (timepoint non-specific)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unpack options
% Operational options
overwrite = options.overwrite ;
plot_flows = options.plot_flows ;
plot_Hgdot = options.plot_Hgdot ;
plot_factors = options.plot_factors ;

% parameter options
doubleResolution = options.doubleResolution ;
lambda = options.lambda ;
lambda_err = options.lambda_err ;
lambda_mesh = options.lambda_mesh ;
nmodes = options.nmodes ;
zwidth = options.zwidth ;
H2vn2d = options.H2vn2d ;
divv2d = options.divv2d ;
gdot2d = options.gdot2d ;
veln2d = options.veln2d ;
radi2d = options.radi2d ;
H2d = options.H2d ;
H2vn3d = options.H2vn3d ;
divv3d = options.divv3d ; 
gdot3d = options.gdot3d ;
veln3d = options.veln3d ;
radi3d = options.radi3d ;
H3d = options.H3d ;
cutMesh = options.cutMesh ;  % can be empty to save computation time if not overwrite
mesh = options.mesh ;        % can be empty to save computation time if not overwrite
climit = options.climit ;
climit_err = options.climit_err ;
climit_H = options.climit_H ;
climit_veln = options.climit_veln ;

if doubleResolution
    sresStr = 'doubleRes_' ;
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    sresStr = '' ;
    nU = QS.nU ;
    nV = QS.nV ;
end

%% Unpack QS
% Load time offset for first fold, t0
QS.t0set() ;
tfold = QS.t0 ;
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
mKDir = fullfile(QS.dir.metricKinematics.root, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
    lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
dimDirs = {fullfile(mKDir, 'images_2d'), ...
           fullfile(mKDir, 'images_3d')} ;
% Make sure the directories exist
for qq = 1:length(dimDirs)
    if ~exist(dimDirs{qq}, 'dir')
        mkdir(dimDirs{qq})
    end
end
       
% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];
Hunitstr = [ '[1/' QS.spaceUnits ']' ];
vunitstr = [ '[' QS.spaceUnits '/' QS.timeUnits ']' ];
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Timepoint specific operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr = ['$t=$' sprintf('%03d', tp - tfold) ' ' QS.timeUnits ] ;

% Load meshes if not supplied    
if isempty(mesh) 
    % Check if ALL plot files already exist
    % FLOWS
    fn2 = fullfile(dimDirs{1}, sprintf('incompr_%1dd_%06d.png', 2, tp)) ;
    fn3 = fullfile(dimDirs{2}, sprintf('incompr_%1dd_%06d.png', 3, tp)) ;
    redo_prediction = ~exist(fn2, 'file') || ...
                      ~exist(fn3, 'file') || overwrite ;
    % GDOT_H
    gdot_H_dir = fullfile(mKDir, 'gdot_vs_H') ;
    fn_gdot = fullfile(gdot_H_dir, sprintf('gdot_vn_H_2d_%06d.png', tp)) ;
    % FACTORS
    factorsDir = fullfile(mKDir, 'factors') ;
    fn_factors = fullfile(factorsDir, sprintf('factors_2d_%06d.png', tp)) ;

    % If they don't all already exist, load the meshes
    if (plot_flows && redo_prediction) || ...
            (plot_Hgdot && (~exist(fn_gdot, 'file') || overwrite)) || ...
            (plot_factors && (~exist(fn_factors, 'file') || overwrite))

        % Load current mesh
        if doubleResolution
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC2x, tp)) ;
            mesh = tmp.spcutMeshSmRSC2x ;
        else
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
            mesh = tmp.spcutMeshSmRSC ;
        end
        
        % Smooth the mesh with lambda_mesh
        if lambda_mesh > 0 
            disp('smoothing mesh vertices before computations')
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
        end

        % Load cutMesh also
        if isempty(cutMesh)
            % Load cutMesh
            if doubleResolution
                tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS2x, tp)) ;
                cutMesh = tmp.spcutMeshSmRS2x ;
                clearvars tmp
            else
                tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
                cutMesh = tmp.spcutMeshSmRS ;
                clearvars tmp
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the prediction, the measurement, and the difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors2d = {H2vn2d, divv2d, gdot2d} ; 
colors3d = {H2vn3d, divv3d, gdot3d} ;
set(gcf, 'visible', 'off') ;
% colormap bwr
colormap(twilight_shifted_mod(256))
% Check if files already exist
fn2 = fullfile(dimDirs{1}, sprintf('incompr_%1dd_%06d.png', 2, tp)) ;
fn3 = fullfile(dimDirs{2}, sprintf('incompr_%1dd_%06d.png', 3, tp)) ;
fns = {fn2, fn3} ;
redo_prediction = ~exist(fn2, 'file') || ...
                  ~exist(fn3, 'file') || overwrite ;
if plot_flows && redo_prediction
    for dim = 2:3    
        % Create all panels in 2d or 3d
        for row = 1:2  % row
            for col = 1:3  % column
                % create panel
                subtightplot(3, length(colors2d), col + (row - 1) * 3)

                % If 2d, plot in pullback space
                % if 3d, plot in embedding space
                if dim == 2 && row == 2
                    ss = cutMesh.u(:, 1) ;
                    ssvals = QS.a_fixed * ss / max(ss) ;
                    trisurf(cutMesh.f, ssvals, ...
                        cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
                        colors2d{col}, 'edgecolor', 'none')
                    xlim([0, QS.a_fixed])
                    ylim([0, 1]) 
                    caxis([-climit, climit]) 
                    axis off
                else
                    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                        colors3d{col}, 'edgecolor', 'none')
                    axis equal
                    axis off
                    caxis([-climit, climit]) 

                    % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
                    % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
                    % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
                    xlim(xyzlim(1, :))
                    ylim(xyzlim(2, :))
                    zlim(xyzlim(3, :))
                end

                % Set title and colorbars
                if col == 1 && row == 1
                    % title(['$v_n 2 H$', ...
                    %     ', $t=$' sprintf('%03d', tp - tfold)], ...
                    %     'Interpreter', 'Latex')
                    view(0, 0)
                    caxis([-climit, climit]) 
                elseif col == 2 && row == 1
                    % title(['$\nabla \cdot \bf{v}_\parallel$', ...
                    %     ', $t=$' sprintf('%03d', tp - tfold)], ...
                    %     'Interpreter', 'Latex')
                    title(titlestr, 'Interpreter', 'Latex')
                    view(0, 0)
                    caxis([-climit, climit]) 
                elseif col == 3 && row == 1
                    % title(['$\textrm{Tr}[g^{-1} \dot{g}]$', ...
                    %     ', $t=$' sprintf('%03d', tp - tfold)], ...
                    %     'Interpreter', 'Latex')
                    view(0, 0)
                    caxis([-climit_err, climit_err]) 
                elseif col == 1 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    % set(cb, 'position',[.17 .1 .2 .03])
                    set(cb, 'position',[.165 .2 .15 .03], ...
                        'XTick', [-climit, 0, climit])     
                    ylabel(cb, ['$v_n 2H$ ' unitstr], ...
                        'Interpreter', 'Latex')                    
                    caxis([-climit, climit]) 
                elseif col == 2 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    % set(cb, 'position',[.63 .1 .2 .03])
                    set(cb, 'position',[.445 .2 .15 .03], ...
                        'XTick', [-climit, 0, climit])
                    ylabel(cb, ...
                        ['$\nabla \cdot \bf{v}_\parallel$ ' unitstr ], ...
                        'Interpreter', 'Latex')
                    caxis([-climit, climit]) 
                elseif col == 3 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.725 .2 .15 .03], ...
                        'XTick', [-climit_err, 0, climit_err])
                    ylabel(cb, ...
                        ['local area change ' unitstr], ...
                        'Interpreter', 'Latex')
                    % xticks(cb, [-climit_err, climit_err]) 
                    caxis([-climit_err, climit_err]) 
               end
            end
        end

        % Save the plot
        fn = fns{dim-1} ;
        disp(['saving ', fn])
        set(gcf, 'Color', 'white')
        export_fig(fn, '-png', '-nocrop', '-r200') 
    end
end

%% Plot the residual with separate factors of vn, H, div, and gdot
gdot_H_dir = fullfile(mKDir, 'gdot_vs_H') ;
fn = fullfile(gdot_H_dir, sprintf('gdot_vn_H_2d_%06d.png', tp)) ;
if plot_Hgdot && (~exist(fn, 'file') || overwrite)
    % Plot mean curvature, normal velocity, divv and error
    colors2d = {veln2d, H2d, divv2d, gdot2d} ; 
    colors3d = {veln3d, H3d, divv3d, gdot3d} ; 
    climits = [climit_veln, climit_H, climit, climit_err] ;
    if ~exist(gdot_H_dir, 'dir')
        mkdir(gdot_H_dir) ;
    end
    for row = 1:2  % row
        for col = 1:4  % column
            % create panel
            subplot(3, length(colors2d), col + (row - 1) * 4)

            % plot in pullback space
            if row == 2
                ss = cutMesh.u(:, 1) ;
                ssvals = QS.a_fixed * ss / max(ss) ;
                trisurf(cutMesh.f, ssvals, ...
                    cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
                    colors2d{col}, 'edgecolor', 'none')
                xlim([0, QS.a_fixed])
                ylim([0, 1]) 
                caxis([-climits(col), climits(col) ]) 
                axis off
            else
                trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                    colors3d{col}, 'edgecolor', 'none')
                axis equal
                axis off
                caxis([-climits(col), climits(col)])

                % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
                % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
                % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
                xlim(xyzlim(1, :))
                ylim(xyzlim(2, :))
                zlim(xyzlim(3, :))
            end

            % Set title and colorbars
            if col == 2 && row == 1
                % 3d view with title
                title(titlestr, 'Interpreter', 'Latex')
                view(0, 0)
            elseif row == 1
                % 3d view
                view(0, 0)
            elseif col == 1 && row == 2
                view(0, 270)
                cb = colorbar('south') ;
                set(cb, 'position',[.15 .2 .12 .03], 'Xtick', ...
                    [-climits(col), 0, climits(col)]) 
                ylabel(cb, ['$v_n$ ' vunitstr], 'Interpreter', 'Latex')
            elseif col == 2 && row == 2
                view(0, 270)
                cb = colorbar('south') ;
                set(cb, 'position',[.355 .2 .12 .03], 'Xtick', ...
                    [-climits(col), 0, climits(col)])
                ylabel(cb, ['$H$ ' Hunitstr], ...
                    'Interpreter', 'Latex')
            elseif col == 3 && row == 2
                view(0, 270)
                cb = colorbar('south') ;
                set(cb, 'position',[.565 .2 .12 .03], 'Xtick', ...
                    [-climits(col), 0, climits(col)])
                ylabel(cb, ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
                    'Interpreter', 'Latex')
            elseif col == 4 && row == 2
                view(0, 270)
                cb = colorbar('south') ;
                set(cb, 'position',[.765 .2 .12 .03], 'Xtick', ...
                    [-climits(col), 0, climits(col)])
                ylabel(cb, ['$\textrm{Tr}[g^{-1} \dot{g}]$ ' unitstr], ...
                    'Interpreter', 'Latex')
           end
           caxis([-climits(col), climits(col)])
        end
    end
    disp(['saving ', fn])
    set(gcf, 'Color', 'white')
    export_fig(fn, '-png', '-nocrop', '-r200')    

end

%% Plot the factors separately of vn, 2H, and vn*2H
factorsDir = fullfile(mKDir, 'factors') ;
fn = fullfile(factorsDir, sprintf('factors_2d_%06d.png', tp)) ;
if plot_factors && (~exist(fn, 'file') || overwrite)
    % ensure output directory
    if ~exist(factorsDir, 'dir')
        mkdir(factorsDir)
    end

    % Plot mean curvature, normal velocity, and product
    colors2d = {veln2d, H2d, H2vn2d} ; 
    colors3d = {veln3d, H3d, H2vn3d} ; 
    climits = [climit_veln, climit_H, climit] ;
    for row = 1:2  % row
        for col = 1:3  % column
            % create panel
            subplot(3, length(colors2d), col + (row - 1) * 3)

            % plot in pullback space
            if row == 2
                ss = cutMesh.u(:, 1) ;
                ssvals = QS.a_fixed * ss / max(ss) ;
                trisurf(cutMesh.f, ssvals, ...
                    cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
                    colors2d{col}, 'edgecolor', 'none')
                xlim([0, QS.a_fixed])
                ylim([0, 1]) 
                caxis([-climits(col), climits(col) ]) 
                axis off
            else
                trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                    colors3d{col}, 'edgecolor', 'none')
                axis equal
                axis off
                caxis([-climits(col), climits(col)])

                % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
                % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
                % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
                xlim(xyzlim(1, :))
                ylim(xyzlim(2, :))
                zlim(xyzlim(3, :))
            end

            % Set title and colorbars
            if col == 2 && row == 1
                % 3d view with title
                title(titlestr, 'Interpreter', 'Latex')
                view(0, 0)
            elseif row == 1
                % 3d view
                view(0, 0)
            elseif col == 1 && row == 2
                view(0, 270)
                cb = colorbar('south') ;
                set(cb, 'position',[.165 .2 .15 .03])     
                ylabel(cb, ['$v_n$ ' vunitstr], 'Interpreter', 'Latex')
            elseif col == 2 && row == 2
                view(0, 270)
                cb = colorbar('south') ;
                set(cb, 'position',[.445 .2 .15 .03])
                ylabel(cb, ['$H$ ' Hunitstr ], ...
                    'Interpreter', 'Latex')
            elseif col == 3 && row == 2
                view(0, 270)
                cb = colorbar('south') ;
                set(cb, 'position',[.725 .2 .15 .03])
                ylabel(cb, ['$v_n 2H$ ' unitstr], ...
                    'Interpreter', 'Latex') 
           end
           caxis([-climits(col), climits(col)])
        end
    end

    % Save figure
    disp(['saving ', fn])
    set(gcf, 'Color', 'white')
    export_fig(fn, '-png', '-nocrop', '-r200')    
end
