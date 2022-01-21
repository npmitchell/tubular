function coordinateSystemDemo(QS)
%COORDINATESYSTEMDEMO(QS)
%   Draw coordinate system for currentTime stamp, presentation/publication
%   style
%   currently written for spcutMesh --> modify for spcutMeshSm
%
%
% NPMitchell 2020

% Unpack QS
nU = QS.nU ;
nV = QS.nV ;

% First make zeta coordinate system chart, then phi
for zp = 1:2
    close all
    set(gcf, 'visible', 'off')
    
    % Used to say:
    % spcm = QS.currentMesh.spcutMesh ;
    % Now do:
    spcm = QS.currentMesh.spcutMeshSm ;
    
    % Rotate vertices and prep figure colors
    xyzrs = xyz2APDV(QS, spcm.v) ;
    npanel = 4 ;
    colormap default ;
    zcolors = parula(nU) ;
    % create handles for all subplots
    h(1) = subplot(1, npanel, 1:2) ;
    h(2) = subplot(1, npanel, 3) ;
    h(3) = subplot(1, npanel, 4) ;
    % add to subplots
    for qq = 1:nU
        if zp == 1
            inds = qq:nU:nU*nV ;
        else
            inds = (1:nU) + (qq-1)*nV ;
        end
        color = zcolors(qq, :) ;

        % Make each panel showing coordinate system
        set(gcf, 'CurrentAxes', h(1))
        plot3(xyzrs(inds, 1), xyzrs(inds, 2), xyzrs(inds, 3), '-', 'Color', color)
        hold on;
        set(gcf, 'CurrentAxes', h(2))
        plot(spcm.uv(inds, 1) / QS.a_fixed, spcm.uv(inds, 2), '-', 'Color', color)
        hold on;
        set(gcf, 'CurrentAxes', h(3))
        plot(spcm.sphi(inds, 1), spcm.sphi(inds, 2)*2*pi, '-', 'Color', color)
        hold on;
    end
    % Format axis 1
    set(gcf, 'CurrentAxes', h(1))
    [~, ~, xyzlim_um] = getXYZLims(QS) ;
    xlim(xyzlim_um(1, :))
    ylim(xyzlim_um(2, :))
    zlim(xyzlim_um(3, :))
    axis equal
    xlabel('ap position [$\mu$m]', 'Interpreter', 'Latex')
    ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
    zlabel('dv position [$\mu$m]', 'Interpreter', 'Latex')

    % Format axis 2
    set(gcf, 'CurrentAxes', h(2))
    pos = get(h(2),'Position');
    set(h(2),'Position', pos + [.03 0 0 0]);
    pbaspect([1, 1, 1])
    yticks([0, 1])
    ylim([0, 1])
    xticks([0, 1])
    xlim([0, 1])
    xlabel('$u$', 'Interpreter', 'Latex')
    ylabel('$v$', 'Interpreter', 'Latex', 'rot', 0)

    % Format axis 3
    set(gcf, 'CurrentAxes', h(3))
    pos = get(h(3),'Position');
    set(h(3),'Position', pos + [.06 0 0 0]);
    pbaspect([1, 1, 1])
    yticks([0, 2*pi])
    ylim([0, 2*pi])
    yticklabels({'0','2\pi'})
    xlabel('$\zeta$ [$\mu$m]', 'Interpreter', 'Latex')
    ylabel('$\phi$', 'Interpreter', 'Latex', 'rot', 0)
    if zp == 1
        outfn = fullfile(QS.dir.spcutMesh, ...
            'coordinate_system_demo', 'coordinate_system_demo_zeta') ;
    else
        outfn = fullfile(QS.dir.spcutMesh, ...
            'coordinate_system_demo', 'coordinate_system_demo_phi') ;
    end
    reses = [100, 300] ;
    for rr = 1:length(reses)
        res = reses(rr) ;
        exten = ['_res' num2str(res) '_' num2str(QS.currentTime) ] ;
        export_fig([outfn exten '.png'], ...
            '-nocrop', '-transparent', ['-r' num2str(res)])
    end
end
close all