%% Demonstrate the polar Writhe definition on example curves
%
% NPMitchell 2019
%
outdir = 'polarWrithe_demo' ;
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
addpath(fullfile('..','plotting'))

%% First example: half twist, then twist back
t = 0:0.001:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xx = [xx xx(1:end) ];
yy = [yy fliplr(yy(1:end)) ];
zz = [zz fliplr(zz(1:end)) ];
xx = xx(1:5:end) ;
yy = yy(1:5:end) ;
zz = zz(1:5:end) ;
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;
[WrL, wrL] = writheLevitt(xyz, true) ;
[WrG, wrG] = writheGaussIntegral(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:turns(1) ;
seg2 = turns(1):length(xx) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
hold on
plot3(xx(seg2), yy(seg2), zz(seg2), 'o') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr$_p$ = ' sprintf('%0.3f', wr) ];
titletext = [titletext ', Wr$_G$ = ' sprintf('%0.3f', WrG) ] ;
titletext = [titletext ', Wr$_L$ = ' sprintf('%0.3f', WrL) ] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', ...
      'Interpreter', 'Latex') ;
axis off
saveas(gcf, fullfile(outdir, 'twoseg_helix.png'))
close all;
plot(wrL)
hold on
plot(wrG)
plot(wr_local)
legend({'Levitt', 'Gauss', 'polar'})
xlabel('index')
ylabel('Writhe density')

%% Half twist
t = 0:0.002 :1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:length(xyz) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
axis off
saveas(gcf, fullfile(outdir, 'oneseg_helix.png'))
close all


%% Half twist reverse
t = 0:0.01:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xx = [xx ];
yy = [fliplr(yy) ];
zz = [fliplr(zz) ];
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:length(xyz) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
axis off
saveas(gcf, fullfile(outdir, 'oneseg_helix_reverse.png'))
close all


%% half twist, then straight back
t = 0:0.01:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xback = linspace(xx(end), xx(1), 10) ;
yback = linspace(yy(end), yy(1), 10) ;
zback = linspace(zz(end), zz(1), 10) ;
xx = [xx xback ];
yy = [yy yback ];
zz = [zz zback ];
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:turns(1) ;
seg2 = turns(1):length(xx) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
hold on
plot3(xx(seg2), yy(seg2), zz(seg2), 'o') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
axis off
saveas(gcf, fullfile(outdir, 'twoseg_helix.png'))
close all


%% half twist that turns back
t = 0:0.01:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = 1.25 - (t - 0.5).^2 ;
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:turns(1) ;
seg2 = turns(1):length(xx) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
hold on
plot3(xx(seg2), yy(seg2), zz(seg2), 'o') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
 axes( 'Position', [0, 0, 1, 1] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
saveas(gcf, fullfile(outdir, 'benthelix.png'))
close all

%% Incomplete torus wound curve
resl = [0.05, 0.03, 0.02, 0.01, 0.005, 0.001] ;
Wrs = [] ;
WrLs = [] ;
WrGs = [] ;
resID = 0 ;
for res = resl
    resID = resID + 1;
    t = 0:res:1 ;
    us = {t, t} ;
    vs = {2 * t, 3 * t} ;
    for ii = 1 
        u = us{ii} ;
        v = vs{ii} ;
        yy = cos(2*pi * v) .* (2 + cos(2 * pi * u)) ;
        zz = sin(2*pi * v) .* (2 + cos(2 * pi * u)) ;
        xx = sin(2*pi * u) ;
        xyz = [xx ; yy; zz]';
        [wr, wr_local, wr_nonlocal, turns, segs, segpairs] = polarWrithe(xyz, [], res) ;
        [WrL] = writheLevitt(xyz, false) ;
        [WrG, wrG] = writheGaussIntegral(xyz, []) ;

        Wrs(resID) = wr ;
        WrLs(resID) = WrL ; 
        WrGs(resID) = WrG ;

        colors = define_colors(length(segs)) ;
        markers = define_markers(length(segs)) ;

        subplot(2, 2, 1)
        % scatter3(xx', yy', zz', 10, yy) 
        cvals = 1:length(xx) ;
        scatter3(xx', yy', zz', 10, cvals) 
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis equal
        title('Curve')

        % Plot the segments colored
        subplot(2,2,2)
        for jj = 1:length(segs)
            plot3(xx(segs{jj}), yy(segs{jj}), zz(segs{jj}), '.', 'Color', colors(jj, :)) ;
            hold on
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis equal
        title('Segments')
        % Plot the local writhe
        subplot(2, 2, 3)
        plot(t, wr_local)
        title('local writhe')
        xlabel('position, s')
        ylabel('local writhe')
        % Plot the nonlocal writhe
        subplot(2, 2, 4)
        for jj = 1:length(wr_nonlocal)
            segments = segpairs{jj} ;
            segi = segments(1) ;
            segj = segments(2) ;
            plot(segi, wr_nonlocal(jj), markers{segj}, 'Color', colors(segj, :))
            hold on
        end
        title('nonlocal writhe')
        xlabel('segment pair index')
        ylabel('nonlocal writhe')
        axes( 'Position', [0, 0, 1, 1] ) ;
        set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
        titletext = ['Wr = ' num2str(wr)] ;
        text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
              'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
        axis off
        name = strrep(['torusknot_' num2str(ii) '_' sprintf('%0.3f', res)], '.', 'p') ;
        saveas(gcf, fullfile(outdir, [name '.png']))
        close all
    end
end

close all
h1 = plot(resl, Wrs, [markers{1} '-'], 'Color', colors(1, :)) ; 
hold on;
h2 = plot(resl, WrLs, [markers{2} '-'], 'Color', colors(2, :)) ;
% h3 = plot(resl(WrGs < 10), WrGs(WrGs < 10), [markers{3} '-'], 'Color', colors(3, :)) ; 
legend([h1, h2], {'polar', 'Levitt'}, 'location', 'best')
title('Writhe comparison')
xlabel('resolution')
ylabel('Writhe') 
saveas(gcf, fullfile(outdir, ['writhecomparison_curv' num2str(ii) '.png']))


%% Torus wound open curve 
close all
res = 0.005 ;
resl = [.1:0.05:1.0] ;
Wrs = [] ;
WrLs = [] ;
WrGs = [] ;
resID = 0 ;
for leng = resl
    resID = resID + 1;
    t = 0:res:leng ;
    us = {t} ;
    vs = {2 * t} ;
    for ii = 1:length(us) 
        u = us{ii} ;
        v = vs{ii} ;
        yy = cos(2*pi * v) .* (2 + cos(2 * pi * u)) ;
        zz = sin(2*pi * v) .* (2 + cos(2 * pi * u)) ;
        xx = sin(2*pi * u) ;
        xyz = [xx ; yy; zz]';
        [wr, wr_local, wr_nonlocal, turns, segs, segpairs] = polarWrithe(xyz, [], res) ;
        [WrL] = writheLevitt(xyz, true) ;
        [WrG, wrG] = writheGaussIntegral(xyz, []) ;

        Wrs(resID) = wr ;
        WrLs(resID) = WrL ; 
        WrGs(resID) = WrG ;

        colors = define_colors(length(segs)) ;
        markers = define_markers(length(segs)) ;

        subplot(2, 2, 1)
        % scatter3(xx', yy', zz', 10, yy) 
        cvals = 1:length(xx) ;
        scatter3(xx', yy', zz', 10, cvals) 
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis equal
        title('Curve')

        % Plot the segments colored
        subplot(2,2,2)
        for jj = 1:length(segs)
            plot3(xx(segs{jj}), yy(segs{jj}), zz(segs{jj}), '.', 'Color', colors(jj, :)) ;
            hold on
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis equal
        title('Segments')
        % Plot the local writhe
        subplot(2, 2, 3)
        plot(zz, wr_local)
        title('local writhe')
        xlabel('position, z')
        ylabel('local writhe')
        % Plot the nonlocal writhe
        subplot(2, 2, 4)
        for jj = 1:length(wr_nonlocal)
            segments = segpairs{jj} ;
            segi = segments(1) ;
            segj = segments(2) ;
            plot(segi, wr_nonlocal(jj), markers{segj}, 'Color', colors(segj, :))
            hold on
        end
        title('nonlocal writhe')
        xlabel('segment pair index')
        ylabel('nonlocal writhe')
        axes( 'Position', [0, 0, 1, 1] ) ;
        set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
        titletext = ['Wr = ' num2str(wr)] ;
        text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
              'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
        axis off
        name = strrep(['torusincomplete_' num2str(ii) '_' sprintf('%0.3f', leng)], '.', 'p') ;
        saveas(gcf, fullfile(outdir, [name '.png']))
        close all
    end
end

close all
h1 = plot(resl, Wrs, [markers{1} '-'], 'Color', colors(1, :)) ; 
hold on;
h2 = plot(resl, WrLs, [markers{2} '-'], 'Color', colors(2, :)) ;
% h3 = plot(resl(WrGs < 10), WrGs(WrGs < 10), [markers{3} '-'], 'Color', colors(3, :)) ; 
legend([h1, h2], {'polar', 'Levitt'}, 'location', 'best')
title('Writhe comparison')
xlabel('length of open curve')
ylabel('Writhe') 
saveas(gcf, fullfile(outdir, 'writhecomparison_incompletetorus.png'))

%% Helix segment
close all
Wrs = [] ;
WrLs = [] ;
WrGs = [] ;
resID = 0 ;
resl = [0.05:0.05:1] ;
for leng = resl
    resID = resID + 1;
    t = 0:0.002:leng ;
    xx = sin(2*pi*t) ;
    yy = cos(2*pi*t) ;
    zz = t ;
    xyz = [xx ; yy; zz]';

    [wr, wr_local, wr_nonlocal, turns, segs, segpairs] = polarWrithe(xyz, [], res) ;
    [WrL] = writheLevitt(xyz, true) ;
    [WrG, wrG] = writheGaussIntegral(xyz, []) ;

    Wrs(resID) = wr ;
    WrLs(resID) = WrL ; 
    WrGs(resID) = WrG ;

    colors = define_colors(length(segs)) ;
    markers = define_markers(length(segs)) ;

    subplot(2, 2, 1)
    % scatter3(xx', yy', zz', 10, yy) 
    cvals = 1:length(xx) ;
    scatter3(xx', yy', zz', 10, cvals) 
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    title('Curve')

    % Plot the segments colored
    subplot(2,2,2)
    for jj = 1:length(segs)
        plot3(xx(segs{jj}), yy(segs{jj}), zz(segs{jj}), '.', 'Color', colors(jj, :)) ;
        hold on
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    title('Segments')
    % Plot the local writhe
    subplot(2, 2, 3)
    plot(zz, wr_local)
    title('local writhe')
    xlabel('position, z')
    ylabel('local writhe')
    % Plot the nonlocal writhe
    subplot(2, 2, 4)
    for jj = 1:length(wr_nonlocal)
        segments = segpairs{jj} ;
        segi = segments(1) ;
        segj = segments(2) ;
        plot(segi, wr_nonlocal(jj), markers{segj}, 'Color', colors(segj, :))
        hold on
    end
    title('nonlocal writhe')
    xlabel('segment pair index')
    ylabel('nonlocal writhe')
    axes( 'Position', [0, 0, 1, 1] ) ;
    set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
    titletext = ['Wr = ' num2str(wr)] ;
    text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    axis off
    name = strrep(['torusincomplete_' num2str(ii) '_' sprintf('%0.3f', leng)], '.', 'p') ;
    saveas(gcf, fullfile(outdir, [name '.png']))
    close all
end

close all
h1 = plot(resl, Wrs, [markers{1} '-'], 'Color', colors(1, :)) ; 
hold on;
h2 = plot(resl, WrLs, [markers{2} '-'], 'Color', colors(2, :)) ;
% h3 = plot(resl(WrGs < 10), WrGs(WrGs < 10), [markers{3} '-'], 'Color', colors(3, :)) ; 
legend([h1, h2], {'polar', 'Levitt'}, 'location', 'best')
title('Writhe comparison')
xlabel('length of open curve')
ylabel('Writhe') 
saveas(gcf, fullfile(outdir, 'writhecomparison_incompletetorus.png'))
