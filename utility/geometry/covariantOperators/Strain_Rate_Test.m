%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strain (rate) calculation Test
%
%   Build a tube with a velocity field. Compute the strain(rate) on the
%   mesh. Compare with analytic results.
%   Everything makes sense here except the introduction of factors of 
%   1/2pi in the off-diagonal components of the strain rate during twisting
%   deformations. I suppose the idea is that for unit strain, we must twist
%   by 2*pi*R, so the velocity would have to be 2*pi*r*hat(theta).
%
% Results: 
%   - twisted tube: theta = atan2(1,2*pi*R) = 0.1578 if R=1. Good!
%   - twist + elongation:
%       If twist velocity has magnitude delta so epsilon_rtheta is 
%       delta/2R, then I predict:
%       theta=atan2(delta, 4pi^2R^3+2pi[sqrt(R^2(4pi^2R^4+delta^2)])
%       theta = 0.1258 if epsilon = (1, 1/2R; 1/2R, 0) ;
%       twist and elongation are equal magnitude (1). 
%       Otherwise, if twist is delta/elongation in magnitude, we get a
%       increasing function of delta.
%     --> Numerically, we get theta=0.0657 unless we divide the
%     circumferential flow magnitude by 2*pi.
%   - small twist (delta=0.01) + unit elongation: 
%       theta is increasing function of delta 
%     --> Numerically, we get theta=0.18437
%   - Canted tube: theta = atan2(1, 2*pi*R) or atan2(-1, 2*pi*R),
%   depending on angle around the tube. For <pi/2, first result applies.
%     --> Numerically we get perfect agreement
%   - Contraction returns pi/2 exactly.
%   - Extension returns as predicted
%
% NPMitchell 2020, with Dillon Cislo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
radius = 1 ;     % radius of the tube
pausetime = 1 ;  % seconds to wait while displaying result
nU = 10 ;  % how many faces along the tube's length
nV = 10 ;  % how many faces along the tube's circumference
climit = 1 ;  % color limit for plotting strain rate

%% Colormaps
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;
    
% Consider both a simple tube and a curvy tube
for ii = 1:2
    %% Make a tube 
    [cutMesh, mesh] = buildCylinderMesh(nU, nV) ;
    vtx = mesh.v ;
    phi = atan2(vtx(:, 3), vtx(:, 2)) ;
    
    % On second pass, curve the tube
    if ii == 2
        cutMesh.v(:, 1) = 2 * cutMesh.v(:, 1) - 1 ;
        cutMesh.v(:, 2) = cutMesh.v(:, 2) .* (1 + cutMesh.v(:, 1).^2) ;
        cutMesh.v(:, 3) = cutMesh.v(:, 3) .* (1 + cutMesh.v(:, 1).^2) ;
        mesh.v(:, 1) = 2 * mesh.v(:, 1) - 1 ;
        mesh.v(:, 2) = mesh.v(:, 2) .* (1 + mesh.v(:, 1).^2) ;
        mesh.v(:, 3) = mesh.v(:, 3) .* (1 + mesh.v(:, 1).^2) ;
    end
    
    %% Define velocities on the vertices -- model flows for comparison
    % twisting shear
    vscale = 1 ;  % 2*pi*radius ;
    vs{1} = [0*vtx(:, 1), - vscale * vtx(:, 1) .* sin(phi(:)), ...
        vscale * vtx(:, 1) .* cos(phi(:))] ;
    % twisting shear + elongation
    vs{2} = [vtx(:, 1), -vtx(:, 1) .* sin(phi(:)), vtx(:, 1) .* cos(phi(:))] ;
    % twisting shear + elongation
    vs{2} = [vtx(:, 1), -vtx(:, 1) .* sin(phi(:)), vtx(:, 1) .* cos(phi(:))] ;
    % delta twisting shear + unit elongation
    delta = 1/(2*pi) ; alpha=1;
    vs{3} = [alpha * vtx(:, 1), -delta * vtx(:, 1) .* sin(phi(:)), ...
        delta * vtx(:, 1) .* cos(phi(:))] ;
    % tilting shear
    vs{4} = [vtx(:, 3), 0*vtx(:, 2), 0*vtx(:, 3)] ;
    % contraction -- quadratic
    vx = sign(mean(vtx(:, 1)) - vtx(:, 1)) .* (mean(vtx(:, 1)) - vtx(:, 1)).^2 ;
    vs{5} = [vx, 0*vx, 0*vx] ;
    % elongation -- uniform
    vs{6} = [vtx(:, 1), 0*vx, 0*vx] ;
    % elongation and contraction -- quadratic
    vx = (mean(vtx(:, 1)) - vtx(:, 1)).^2 ;
    vs{7} = [vx, 0 * ones(size(vx)), 0*vx] ;
    % circumferential dilation
    vs{8} = mesh.vn ;
    
    labels = {'twisting shear', ...
        'twisting shear + axial elongation', ...
        'delta twisting shear + unit axial elongation', ...
        'tilting shear', ...
        'quadratic contraction', ...
        'uniform axial elongation', ...
        'quadratic elongation and contraction', ...
        'circumferential dilation'} ;
    
    %% Consider each kind of deformation in turn
    for jj = 1:length(vs)
        
        % Visualize the velocity
        close all; set(gcf, 'visible', 'on')
        subplot(2, 2, 1)
        trisurf(cutMesh.f, cutMesh.u(:, 1), cutMesh.u(:, 2), ...
            0*cutMesh.u(:, 2), 'edgecolor', 'k', 'facecolor', 'none')
        view(2)
        axis equal; xlabel('x'); ylabel('y'); zlabel('z') ;
        title("Mesh vertices in pullback")
        subplot(2, 2, 2)
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            vs{jj}(:, 1), 'edgecolor', 'none')
        title('velocity, $v_x$', 'interpreter', 'latex')
        caxis([-1, 1])
        colormap(bwr256)
        axis equal; xlabel('x'); ylabel('y'); zlabel('z') ;
        subplot(2, 2, 3)
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            vs{jj}(:, 2), 'edgecolor', 'none')
        title('velocity, $v_y$', 'interpreter', 'latex')
        caxis([-1, 1])
        colormap(bwr256)
        axis equal; xlabel('x'); ylabel('y'); zlabel('z') ;
        subplot(2, 2, 4)
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            vs{jj}(:, 3), 'edgecolor', 'none')
        title('velocity, $v_z$', 'interpreter', 'latex')
        caxis([-1, 1])
        colormap(bwr256)
        axis equal; xlabel('x'); ylabel('y'); zlabel('z') ;
        sgtitle(labels{jj})
        pause(pausetime)    
        clf

        %% Obtain mean curvature H for checking against trace(b_ij)
        DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
        H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
        H2d = H3d ;
        H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;

        % Push vectors onto faces
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        vf = V2F * vs{jj} ;

        % Compute divergence and face-based mean curvature Hf
        Hf = V2F * H3d ;
        divv = DEC.divergence(vf) ;
        divf = V2F * divv ;

        % Compute the strain rate tensor
        disp('Computing covariantDerivative')
        [v0n, vf] = resolveTangentNormalVector(cutMesh.f, cutMesh.v, vf) ;
        if suppress_normal{jj}
            v0n = 0*v0n ;
        end
        [~, dvi] = vectorCovariantDerivative(vf, cutMesh.f, cutMesh.v, cutMesh.u) ;
        dvij = cell(size(dvi)) ;
        for qq = 1:length(dvi)
            dvij{qq} = 0.5 * ( dvi{qq} + dvi{qq}' ) ;
        end

        % Compute the second fundamental form
        [gg, ~] = constructFundamentalForms(cutMesh.f, cutMesh.v, cutMesh.u) ;
        [~, bb] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;

        % Strain rate tensor
        strainrate = cell(size(dvi)) ;
        tre = zeros(size(dvi)) ;
        dev = zeros(size(dvi)) ;
        theta = zeros(size(dvi)) ;
        checkTr = zeros(size(dvi)) ;
        checkH = zeros(size(dvi)) ;
        checkdiv = zeros(size(dvi)) ;
        for qq = 1:length(dvi)
            strainrate{qq} = dvij{qq} - v0n(qq) .* bb{qq} ;
            gq = gg{qq} ;
            eq = strainrate{qq} ;
            
            %% traceful component -- 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
            tre(qq) = trace(inv(gq) * eq) ;
            
            %% deviatoric component -- 
            % || epsilon - 1/2 Tr[g^{-1} epsilon] g|| = sqrt(Tr[A A^T]),
            % where A = epsilon - 1/2 Tr[g^{-1} epsilon] g.
            AA = eq - 0.5 * tre(qq) * gq ;
            dev(qq) = sqrt(trace(inv(gq) * (AA * (inv(gq) * AA)))) ;

            %% angle of elongation -- first take eigvectors
            [evec_dev, evals_dev] = eig( inv(gq) * AA ) ;
            % [evec_dev, evals_dev] = eig( AA * inv(gq) ) ;
            [~, idx] = sort(diag(evals_dev)) ;
            evec_dev = evec_dev(:, idx) ;
            pevec = evec_dev(:, end) ;
            theta(qq) = atan2(real(pevec(2)), real(pevec(1))) ;
            
            %% Other checks -- trace, mean curvature, divergence
            checkTr(qq) = trace(inv(gg{qq}) * dvij{qq}) - v0n(qq) .* Hf(qq) ;
            checkH(qq) = trace(inv(gg{qq}) * bb{qq}) ;        % == 2 * Hf(qq) ; 
            checkdiv(qq) = trace(inv(gg{qq}) * dvij{qq}) ;    % == divf(qq) ; 
        end
        % Theta is a nematic -- modulo by pi
        theta = mod(theta, pi) ;
        % 
        % % Check that trace of bij is 2H and that Trace of strain is divergence
        % clf; 
        % set(gcf, 'visible', 'on')
        % subplot(2, 1, 1)
        % plot(checkdiv, divf, '.')
        % hold on;
        % plot(checkdiv, checkdiv, 'k--')
        % axis equal
        % xlabel('Tr$\nabla_i v_j$', 'Interpreter', 'Latex')
        % ylabel('$\nabla \cdot \mathbf{v}$', 'Interpreter', 'Latex')
        % subplot(2, 1, 2)
        % plot(checkH, 2* Hf, 'o') 
        % hold on;
        % plot(checkH, checkH, 'k--') 
        % axis equal
        % xlabel('Tr$b_{ij}$', 'interpreter', 'latex')
        % ylabel('$2H$', 'interpreter', 'latex')
        % sgtitle(labels{jj})
        % pause(pausetime)
        % clf

        %% Check Mean curvature
        clim = max(max(abs(checkH)), 2 * max(abs(Hf))) ;
        subplot(2, 2, 1)
        trisurf(triangulation(mesh.f, mesh.v), checkH, 'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title('Tr$\left[g^{-1} b\right]$', 'interpreter', 'latex')
        subplot(2, 2, 2)
        trisurf(triangulation(mesh.f, mesh.v), 2 * Hf, 'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title('$2H$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        trisurf(triangulation(mesh.f, mesh.v), checkH - 2 * Hf, 'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title('Tr$[g^{-1}b] - 2H$', 'interpreter', 'latex')
        colormap(bwr)
        sgtitle('Mean curvature calculation')
        pause(pausetime)
        clf

        %% Check div(v)
        clf
        clim = max(mean(abs(checkdiv)), mean(abs(divf))) ;
        subplot(2, 2, 1)
        trisurf(triangulation(mesh.f, mesh.v), checkdiv, ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        checkDivStr = '$\frac{1}{2}$Tr$\left[g^{-1} \left(\nabla_i v_j + \nabla_j v_i \right)\right]$' ;
        title(checkDivStr, 'interpreter', 'latex')
        subplot(2, 2, 2)
        trisurf(triangulation(mesh.f, mesh.v), divf, 'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title('$\nabla \cdot \mathbf{v}_{\parallel}$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        trisurf(triangulation(mesh.f, mesh.v), checkdiv - divf, ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title([checkDivStr ' $-\nabla \cdot \mathbf{v}_{\parallel}$'], ...
            'interpreter', 'latex')
        colormap(bwr)
        sgtitle(labels{jj})
        pause(pausetime)
        close all
        
        %% Visualize strain
        clim_trace = min(1, max(abs(tre(:)))) ;
        clim_dev = min(1, max(abs(dev(:)))) ;
        % Panel 1
        subplot(1, 2, 1)
        trisurf(cutMesh.f, ...
            cutMesh.u(:, 1), ...
            cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
            'FaceVertexCData', tre, 'edgecolor', 'none')
        daspect([1,1,1])
        cb = colorbar('location', 'southOutside') ;
        caxis([-clim_trace, clim_trace])
        colormap(bwr256)
        axis off ; view(2)
        title('$\frac{1}{2}$Tr$\left[g^{-1} \varepsilon\right]$', ...
            'interpreter', 'latex')
        % Panel 2 
        subplot(1, 2, 2)
        % Intensity from dvtre and color from the theta
        indx = max(1, round(mod(2*theta, 2*pi)*size(pm256, 1)/(2 * pi))) ;
        colors = pm256(indx, :) ;
        colors = min(dev / clim_dev, 1) .* colors ;
        trisurf(cutMesh.f, cutMesh.u(:, 1), ...
            cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
            'FaceVertexCData', colors, 'edgecolor', 'none')
        daspect([1,1,1])  
        title('$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$', ...
            'interpreter', 'latex') ;
        % Colorbar and phasewheel
        colormap(gca, phasemap)
        phasebar('colormap', phasemap, ...
            'location', [0.82, 0.12, 0.1, 0.135], 'style', 'nematic')
        axis off
        view(2)
        ax = gca ;
        get(gca, 'position')
        cb = colorbar('location', 'southOutside') ;
        drawnow
        axpos = get(ax, 'position') ;
        cbpos = get(cb, 'position') ;
        set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
        set(ax, 'position', axpos) 
        hold on;
        caxis([0, clim_dev])
        colormap(gca, gray)
        sgtitle(labels{jj})
        pause(pausetime)
        clf
        
        %% Plot theta as heatmap
        clf
        trisurf(cutMesh.f, cutMesh.u(:, 1), ...
            cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
            'FaceVertexCData', theta, 'edgecolor', 'none')
        colorbar()
        title('deviatoric component angle, $\theta$', 'interpreter', 'latex')
        view(2)
        pause(pausetime)
        
    end
end
