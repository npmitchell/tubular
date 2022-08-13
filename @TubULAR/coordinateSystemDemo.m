function fig = coordinateSystemDemo(tubi, options)
%COORDINATESYSTEMDEMO(tubi)
%   Draw coordinate system for currentTime stamp, presentation/publication
%   style.
%   Image for publication/presentation on method & coordinate system
%   Create coordinate system charts visualization using smoothed meshes
%
% Parameters
% ----------
% tubi : TubULAR class instance
% options : struct with fields
%   coordSys : char (default='sphism', 'sphi', 'uv')
%   style : 'curves' or 'surface'
%   exten : '.png' or '.pdf' or '.jpg'
%   interpreter : 'latex' or 'default'
%       whether to use the latex interpreter for axes labels
%   normLongitudinal : bool
%   includeCenterline : bool 
%   fillHoops : bool
%
% Returns
% -------
% <none>
%
% Saves to disk
% -------------
% figure saved to fullfile(tubi.dir.spcutMesh, ...
%            'coordinate_system_demo', 'coordinate_system_demo_zeta') ;
% figure saved to fullfile(tubi.dir.spcutMesh, ...
%            'coordinate_system_demo', 'coordinate_system_demo_phi') ;
%
% NPMitchell 2020-2022

    % Default options
    close all
    exten = '.png' ;
    style = 'mesh' ;
    coordSys = 'sphi' ;
    interpreter = 'tex' ;
    normLongitudinal = true ;
    includeCenterline = false ;
    browncolor = [138, 64, 23] ./ 255 ;
    subU = 1 ; % subsampling factor along U
    subV = 1 ; % subsampling factor along V
    fillHoops = false ;
    axisOn = true ;
    
    % Unpack options
    if nargin < 2
        options = struct() ;
    end
    
    if isfield(options, 'exten')
        exten = options.exten ;
    end
    if isfield(options, 'style')
        style = options.style ;
    end
    if isfield(options, 'coordSys')
        coordSys = options.coordSys ;
    end
    if isfield(options, 'interpreter')
        interpreter = options.interpreter ;
    end
    if isfield(options, 'subU')
        subU = options.subU ;
    end
    if isfield(options, 'subV')
        subV = options.subV ;
    end
    if isfield(options, 'includeCenterline')
        includeCenterline = options.includeCenterline ;
    end
    if isfield(options, 'fillHoops')
        fillHoops = options.fillHoops ;
    end
    if isfield(options, 'axisOn')
        axisOn = options.axisOn ;
    end
    
    % add to exten
    if (subU ~= 1 || subV ~=1) && strcmpi(style, 'curves')
        exten = [sprintf('_%03dsU_%03dsV', subU, subV) exten] ;
    end
    if fillHoops && strcmpi(style, 'curves')
        exten = ['_fillHoops' exten] ;
    end
    if includeCenterline
        exten = ['_wcline' exten] ;
    end
    
    % Get mesh
    if strcmpi(coordSys, 'sphism') || strcmpi(coordSys, 'spsm')
        mesh = tubi.getCurrentSPCutMeshSmRS() ;
        coordSys = 'sphism' ;
    elseif strcmpi(coordSys, 'sphi') || strcmpi(coordSys, 'sp')
        mesh = tubi.getCurrentSPCutMesh() ;
        mesh.u = mesh.sphi ;
        mesh.v = tubi.xyz2APDV(mesh.v) ;
        coordSys = 'sphi' ;
    elseif strcmpi(coordSys, 'uv') 
        mesh = tubi.getCurrentUVCutMesh() ;
        mesh.u = mesh.uv ;
        mesh.v = tubi.xyz2APDV(mesh.v) ;
    end
    
    % normalize longitudinal parameterization coordinate
    if normLongitudinal
        mesh.u(:, 1) = mesh.u(:, 1) ./ max(mesh.u(:, 1))  ;                
    end
    
    fig = figure('units', 'centimeters', 'position', [0, 0, 9, 12]) ;
        
    [~,~,~,xyzlims] = tubi.getXYZLims() ;
    
    % plot the mesh    
    if strcmpi(style, 'curves')
        assert(~isempty(tubi.currentTime))
        uu = mesh.u(:, 1) ;
        vv = mesh.u(:, 2) ;
        xx = (mesh.v(:, 1)) ;
        yy = (mesh.v(:, 2)) ;
        zz = (mesh.v(:, 3)) ;
        colorsV = viridis(mesh.nV) ;
        colorsU = viridis(mesh.nU) ;
        nU = mesh.nU ;
        nV = mesh.nV ;

        % AZIMUTH 3D
        subplot(2, 3, [1,2])
        hold off
        for qq = 1:subU:mesh.nU
            plot3(xx(qq:nU:end), yy(qq:nU:end), zz(qq:nU:end), '-', ...
                'color', colorsV(qq, :));
            hold on ;
            if fillHoops
                fill3(xx(qq:nU:end), yy(qq:nU:end), zz(qq:nU:end), ...
                    colorsV(qq, :), 'facealpha', 0.1, 'edgecolor', 'none')
            end
        end
        if includeCenterline
            cline = mesh.avgpts ;
            plot3(cline(:, 1), cline(:, 2), cline(:, 3), ...
                '-', 'linewidth', 2, 'color', browncolor)
        end
        axis equal
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        xlabel([ 'ap position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        ylabel([ 'lateral position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        zlabel([ 'dv position [' tubi.spaceUnits ']'], 'interpreter', interpreter)

        if ~axisOn
            axis off
        end
    
        % LONGITUDE 3D
        subplot(2, 3, [4,5])
        hold off
        for qq = 1:subV:(mesh.nV-1)
            inds = (qq-1)*nU+1:qq*nU ;
            plot3(xx(inds), yy(inds), zz(inds), '-', ...
                'color', colorsU(qq, :));
            hold on ;
        end
        if includeCenterline
            plot3(cline(:, 1), cline(:, 2), cline(:, 3), ...
                '-', 'linewidth', 2, 'color', browncolor)
        end
        axis equal
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        xlabel(['ap position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        ylabel(['lateral position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        zlabel(['dv position [' tubi.spaceUnits ']'], 'interpreter', interpreter)

        if ~axisOn
            axis off
        end

        % LONGITUDE
        subplot(2, 3, 3)
        hold off
        for qq = 1:subU:mesh.nU  
            plot(uu(qq:nU:end), vv(qq:nU:end), '-', ...
                'color', colorsV(qq, :));
            hold on ;
        end        
        
        formatPBSpaceLabels(normLongitudinal, coordSys)
        
        labelAxes_util(coordSys, interpreter, tubi)
        
        % AZIMUTH
        subplot(2, 3, 6)
        hold off
        for qq = 1:subV:mesh.nV
            inds = (qq-1)*nU+1:qq*nU ;
            plot(uu(inds), vv(inds), '-', ...
                'color', colorsU(qq, :));
            hold on ;
        end
        
        % formatting of pullback space
        formatPBSpaceLabels(normLongitudinal, coordSys)
        
    else
        % LONGITUDE 3D
        subplot(2, 3, [1,2])
        hold off
        h1 = trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', ...
            mesh.u(:, 1), 'edgecolor', 'none') ;
        shading interp
        lightangle(-5,30)
        h1.FaceLighting = 'gouraud';
        h1.AmbientStrength = 0.9;
        h1.DiffuseStrength = 0.9;
        h1.SpecularStrength = 0.5;
        h1.SpecularExponent = 5;
        h1.BackFaceLighting = 'unlit';
        colormap viridis
        axis equal
    xlim(xyzlims(1, :))
    ylim(xyzlims(2, :))
    zlim(xyzlims(3, :))
        xlabel([ 'ap position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        ylabel([ 'lateral position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        zlabel([ 'dv position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        grid off

        if ~axisOn
            axis off
        end
    
        % AZIMUTH 3D
        subplot(2, 3, [4,5])
        h2 = trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', ...
            mesh.u(:, 2), 'edgecolor', 'none') ;
        shading interp
        lightangle(-5,30)
        h2.FaceLighting = 'gouraud';
        h2.AmbientStrength = 0.9;
        h2.DiffuseStrength = 0.9;
        h2.SpecularStrength = 0.5;
        h2.SpecularExponent = 5;
        h2.BackFaceLighting = 'unlit';
        colormap viridis
        axis equal
    xlim(xyzlims(1, :))
    ylim(xyzlims(2, :))
    zlim(xyzlims(3, :))
        xlabel([ 'ap position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        ylabel([ 'lateral position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        zlabel([ 'dv position [' tubi.spaceUnits ']'], 'interpreter', interpreter)
        grid off
        if ~axisOn
            axis off
        end

        % LONGITUDE
        subplot(2, 3, 3)
        hold off
        mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
        u3d = [mesh.u, 0*mesh.u(:, 1)] ;
        h3 = trisurf(triangulation(mesh.f, u3d), 'facevertexcdata', ...
            mesh.u(:, 1), 'edgecolor', 'none') ;
        view(2)
        % formatting of pullback space
        formatPBSpaceLabels(normLongitudinal, coordSys)
        
        labelAxes_util(coordSys, interpreter, tubi)

        % AZIMUTH
        subplot(2, 3, 6)
        hold off
        mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
        u3d = [mesh.u, 0*mesh.u(:, 1)] ;
        h4 = trisurf(triangulation(mesh.f, u3d), 'facevertexcdata', ...
            mesh.u(:, 2), 'edgecolor', 'none') ;
        view(2)
        % formatting of pullback space
        formatPBSpaceLabels(normLongitudinal, coordSys)
    end
    
    labelAxes_util(coordSys, interpreter, tubi)
    
    
    set(gcf, 'color', 'w')
    set(gcf, 'position', [0, 0, 9, 12])
    export_fig(fullfile(tubi.dir.uvCoord, [sprintf(...
        'coordSystemDemo_%06d_', tubi.currentTime) style '_' coordSys exten]), ...
        '-nocrop','-r600')

    if nargout > 0
        fig = gcf ;
    end
    
end
            


function labelAxes_util(coordSys, interpreter, tubi)

    if strcmpi(coordSys, 'uv')
        if strcmpi(interpreter', 'latex')
            xlabel('$u$', 'interpreter', interpreter)
            ylabel('$v$', 'interpreter', interpreter, 'rotation', 0)
        else
            xlabel('u', 'interpreter', interpreter)
            ylabel('v', 'interpreter', interpreter, 'rotation', 0)
        end
    else
        if strcmpi(interpreter', 'latex')
            xlabel('$s$', 'interpreter', interpreter)
            ylabel('$\phi$', 'interpreter', interpreter, 'rotation', 0)
        else
            xlabel('s', 'interpreter', interpreter)
            ylabel('\phi', 'interpreter', interpreter, 'rotation', 0)
        end
    end
end


function formatPBSpaceLabels(normLongitudinal, coordSys)
    % formatting of pullback space
    yticks([0, 1])
    if normLongitudinal && ~strcmpi(coordSys, 'uv')    
        xticks([0, 1])
        xticklabels({0, 'L'})
    elseif strcmpi(coordSys, 'uv')
        xticks([0, 1])
    end
    axis square
end
        
%             
% % First make zeta coordinate system chart, then phi
% for zp = 1:2
%     close all
%     set(gcf, 'visible', 'off')
%     
%     if strcmpi(coordSys, 'sphi')
%         spcm = tubi.currentMesh.spcutMesh ;
%     elseif strcmpi(coordSys, 'sphism') || strcmpi(coordSys, 'sphi_sm')
%         spcm = tubi.currentMesh.spcutMeshSm ;
%     elseif strcmpi(coordSys, 'uv')
%         spcm = tubi.getCurrentUVCutMesh ;
%     end
%     
%     % Rotate vertices and prep figure colors
%     xyzrs = xyz2APDV(tubi, spcm.v) ;
%     npanel = 4 ;
%     colormap default ;
%     zcolors = parula(nU) ;
%     % create handles for all subplots
%     h(1) = subplot(1, npanel, 1:2) ;
%     h(2) = subplot(1, npanel, 3) ;
%     h(3) = subplot(1, npanel, 4) ;
%     % add to subplots
%     for qq = 1:nU
%         if zp == 1
%             inds = qq:nU:nU*nV ;
%         else
%             inds = (1:nU) + (qq-1)*nV ;
%         end
%         color = zcolors(qq, :) ;
% 
%         % Make each panel showing coordinate system
%         set(gcf, 'CurrentAxes', h(1))
%         plot3(xyzrs(inds, 1), xyzrs(inds, 2), xyzrs(inds, 3), '-', 'Color', color)
%         hold on;
%         set(gcf, 'CurrentAxes', h(2))
%         plot(spcm.uv(inds, 1) / tubi.a_fixed, spcm.uv(inds, 2), '-', 'Color', color)
%         hold on;
%         set(gcf, 'CurrentAxes', h(3))
%         plot(spcm.sphi(inds, 1), spcm.sphi(inds, 2)*2*pi, '-', 'Color', color)
%         hold on;
%     end
%     % Format axis 1
%     set(gcf, 'CurrentAxes', h(1))
%     [~, ~, xyzlim_um] = getXYZLims(tubi) ;
%     xlim(xyzlim_um(1, :))
%     ylim(xyzlim_um(2, :))
%     zlim(xyzlim_um(3, :))
%     axis equal
%     xlabel('ap position' tubi.spaceUnits], 'Interpreter', 'Latex')
%     ylabel([ 'lateral position' tubi.spaceUnits], 'Interpreter', 'Latex')
%     zlabel([ 'dv position' tubi.spaceUnits], 'Interpreter', 'Latex')
% 
%     % Format axis 2
%     set(gcf, 'CurrentAxes', h(2))
%     pos = get(h(2),'Position');
%     set(h(2),'Position', pos + [.03 0 0 0]);
%     pbaspect([1, 1, 1])
%     yticks([0, 1])
%     ylim([0, 1])
%     xticks([0, 1])
%     xlim([0, 1])
%     xlabel('$u$', 'Interpreter', 'Latex')
%     ylabel('$v$', 'Interpreter', 'Latex', 'rot', 0)
% 
%     % Format axis 3
%     set(gcf, 'CurrentAxes', h(3))
%     pos = get(h(3),'Position');
%     set(h(3),'Position', pos + [.06 0 0 0]);
%     pbaspect([1, 1, 1])
%     yticks([0, 2*pi])
%     ylim([0, 2*pi])
%     yticklabels({'0','2\pi'})
%     xlabel([ '$\zeta$' tubi.spaceUnits], 'Interpreter', 'Latex')
%     ylabel([ '$\phi$', 'Interpreter', 'Latex', 'rot', 0)
%     if zp == 1
%         outfn = fullfile(tubi.dir.spcutMesh, ...
%             'coordinate_system_demo', 'coordinate_system_demo_zeta') ;
%     else
%         outfn = fullfile(tubi.dir.spcutMesh, ...
%             'coordinate_system_demo', 'coordinate_system_demo_phi') ;
%     end
%     reses = [100, 300] ;
%     for rr = 1:length(reses)
%         res = reses(rr) ;
%         exten = ['_res' num2str(res) '_' num2str(tubi.currentTime) ] ;
%         export_fig([outfn exten '.png'], ...
%             '-nocrop', '-transparent', ['-r' num2str(res)])
%     end
% end
% close all