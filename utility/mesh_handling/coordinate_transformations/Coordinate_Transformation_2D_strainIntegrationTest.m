%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Expanding mesh strain rate integration test -- 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 'rotate' ;  % stretchX, rotate
pausetime = 0.1 ;    % amount of time to pause to inspect each timepoint

% Magnitude of strain rate to accumulate
txy = 0.01 ;  % shear strain
sxx = 0.02 ;      % xx strain
syy = 0 ;      % yy strain
r = 0.02 ; % rate of deformation of the mesh per timepoint 

% Size of mesh
nU = 10 ; 
nV = 5 ;

% Define mesh vertices
[YY, XX] = meshgrid(1:nV, 1:nU) ;
XY = [XX(:), YY(:)] ;
ff = defineFacesRectilinearGrid(XY, nU, nV) ;
mm.f = ff ; 
mm.u = XY ;

% Check the mesh
trisurf(mm.f, mm.u(:, 1), mm.u(:, 2), 0*mm.u(:, 2), 'faceColor', 'none')
view(2)
pause(2)

% Colors
close all
pm256 = phasemap(256) ;
ntps = 15 ;
clim = ntps * r ;

% Axes
ax1 = subplot(2, 1, 1) ;
ax2 = subplot(2, 1, 2) ;

for tp = 1:ntps
    %% Accumulate strain rate into STRAIN            
    strainrate = [sxx, txy; txy, syy]  ;
    
    % Warp the previous and next mesh
    mesh0 = mm ;
    mesh1 = mm ;        
    switch lower(method)
        case 'stretchx'
            mesh0.u(:, 1) = mm.u(:, 1) * ((1+r)^tp) ;
            mesh1.u(:, 1) = mm.u(:, 1) * ((1+r)^(tp+1)) ;
        case 'rotate'
            mesh0.u(:, 1) = mm.u(:, 1) * cos(r*tp) - ...
                            mm.u(:, 2) * sin(r*tp) ;
            mesh0.u(:, 2) = mm.u(:, 2) * cos(r*(tp+1)) + ...
                            mm.u(:, 1) * sin(r*(tp+1)) ;    
            mesh1.u(:, 1) = mm.u(:, 1) * cos(r*tp) - ...
                            mm.u(:, 2) * sin(r*tp) ;
            mesh1.u(:, 2) = mm.u(:, 2) * cos(r*(tp+1)) + ...
                            mm.u(:, 1) * sin(r*(tp+1)) ;      
    end
    
    % Check deformation
    % hold on;
    % trisurf(mesh1.f, mesh1.u(:, 1), mesh1.u(:, 2), 5+0*mesh1.u(:, 2), 'faceColor', 'none')
    % drawnow
    
    % Define positions of Lagrangian pathlines in domain of
    % parameterization at this time
    uv = barycenter(mesh1.u, mesh1.f) ;

    % Transform the previous strain rate into current basis
    if tp > 1
        % Time increment
        dt = 1;

        % Construct the Jacobian matrix on each mesh face
        J01 = jacobian2Dto2DMesh(mesh1.u, mesh0.u, mesh1.f) ;
        % J10 = jacobian2Dto2DMesh(mesh0.u, mesh1.u, F);

        % Find which jacobians to use for each pathline point
        tri = triangulation(mesh1.f, mesh1.u) ;
        uv(:, 1) = max(min(mesh1.u(:, 1)) + eps, uv(:, 1)) ;
        uv(:, 2) = max(min(mesh1.u(:, 2)) + eps, uv(:, 2)) ;
        uv(:, 1) = min(max(mesh1.u(:, 1)) - eps, uv(:, 1)) ;
        uv(:, 2) = min(max(mesh1.u(:, 2)) - eps, uv(:, 2)) ;
        fieldfaces = pointLocation(tri, uv) ;

        % Accumulate the strain via implicit Euler (backward Euler
        % scheme)
        % Transform as a (0,2)-tensor (NOTICE THE MATRIX INVERSES)
        for qq = 1:size(uv,1)
            strain0 = [strain(qq, 1), strain(qq, 2); ...
                       strain(qq, 3), strain(qq, 4)] ;
            try
                J01f = J01{fieldfaces(qq)} ;
                strainM{qq} = inv(J01f) * strain0 * (inv(J01f).') +...
                    dt * strainrate ;
            catch
                nanid = find(isnan(fieldfaces)) ;
                hold on; plot(uv(nanid, 1), uv(nanid, 2), '.')
                axis equal
                assert(all(min(fieldfaces) > 0))
                error('Ensure that all uv lie in the mesh.u')
            end
        end

        % Convert from cell (needed for face Jacobians) to array
        strain = zeros(size(strainrate)) ;
        for qq = 1:size(uv, 1)
            strain(qq, 1) = strainM{qq}(1, 1) ;
            strain(qq, 2) = strainM{qq}(1, 2) ;
            strain(qq, 3) = strainM{qq}(2, 1) ;
            strain(qq, 4) = strainM{qq}(2, 2) ; 
        end
    else
        strain = zeros(size(uv, 1), 4) ;
        for qq=1:size(uv, 1)
            strain(qq, :) = strainrate(:) ;
        end
    end
    
    %% Decide on metric for this timepoint
    gxx = 1 * ones(size(uv(:), 1), 1) ;
    gxy = 0 * ones(size(uv(:), 1), 1) ;
    gyx = 0 * ones(size(uv(:), 1), 1) ;
    gyy = 1 * ones(size(uv(:), 1), 1) ;

    %% Trace/Determinant of strain 
    strain_tr = zeros(size(strain, 1), 1) ;
    strain_dv = zeros(size(strain, 1), 1) ;
    strain_theta = zeros(size(strain, 1), 1) ;
    for qq = 1:size(strain, 1)
        eq = [strain(qq, 1), strain(qq, 2); ...
              strain(qq, 3), strain(qq, 4)] ;
        gq = [gxx(qq), gxy(qq); ...
              gyx(qq), gyy(qq)] ;
        [strain_tr(qq), strain_dv(qq), ...
            strain_theta(qq)] = trace_deviator(eq, gq);
    end
    strain_theta= mod(strain_theta, pi) ;
    
    %% Arrange into grid
    % There is one element per face, so we can arrange in (2*nU-2)*(nV-1)
    strain_tr = reshape(strain_tr, [2*nU-2, nV-1]) ;
    strain_dv = reshape(strain_dv, [2*nU-2, nV-1]) ;
    strain_theta = reshape(strain_theta, [2*nU-2, nV-1]) ;
    
    %% PLOT result
    % Map intensity from dvtre and color from the theta
    indx = max(1, round(mod(2*strain_theta(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
    colors = pm256(indx, :) ;
    dvtreKclipped = min(strain_dv / clim, 1) ;
    colorsM = dvtreKclipped(:) .* colors ;
    
    % Show colors as image
    % colorsM = reshape(colorsM, [size(strain_dv, 1), size(strain_dv, 2), 3]) ;
    % imagesc(linspace(0, max(mesh1.u(:, 1)), size(strain_dv, 1)), ...
    %     linspace(0, max(mesh1.u(:, 2)), size(strain_dv, 2)), colorsM)
    
    axes(ax1)
    hh = trisurf(mesh1.f, mesh1.u(:, 1), mesh1.u(:, 2), 0*mesh1.u(:, 2)) ;
    set(hh,'FaceColor','flat',...
       'FaceVertexCData',colorsM,...
       'CDataMapping','scaled')
    view(2) ;
   
    axis equal
    xlim([0, max(mesh1.u(:, 1))])
    caxis([0, clim])
    if tp == 1
        % Colorbar and phasewheel
        colormap(gca, phasemap)
        phasebar('colormap', phasemap, ...
            'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic')
        ax = gca ;
        get(gca, 'position')
        cb = colorbar('location', 'eastOutside') ;
        drawnow
        axpos = get(ax, 'position') ;
        cbpos = get(cb, 'position') ;
        set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
        set(ax, 'position', axpos) 
        hold on;
        caxis([0, clim])
        colormap(gca, gray)
    end
    drawnow 
    
    axes(ax2) 
    h1 = plot(tp, strain(1, 1), 'bo') ;
    h2 = plot(tp, strain(1, 2), 'gs') ;
    h3 = plot(tp, strain(1, 4), 'r^') ;
    hold on;
    
    pause(pausetime)
end

axes(ax1) 
xlim([-1, nU * 1.2])
ylim([-1, nV * 1.5])
title("Strain integration")
axes(ax2) 
legend([h1,h2,h3], {'$\varepsilon_{xx}$', '$\varepsilon_{xy}$', ...
     '$\varepsilon_{yy}$'}, 'Interpreter', 'Latex')
 
saveas(gcf, './strain_integration_test.png')
