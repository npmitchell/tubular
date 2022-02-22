function plotBonds(xy, BL, plotOpts)

Xs = zeros(size(BL, 1), 1) ;
Ys = Xs ;
Us = Xs ;
Vs = Xs ;
for qq = 1:length(BL)            
    Xs(qq) = xy(BL(qq, 1), 1) ;
    Ys(qq) = xy(BL(qq, 1), 2) ; 
    Us(qq) = xy(BL(qq, 2), 1) - xy(BL(qq, 1), 1) ;
    Vs(qq) = xy(BL(qq, 2), 2) - xy(BL(qq, 1), 2) ; 
end
% plot(seg2d.vdat.v(:, 1), seg2d.vdat.v(:, 2), '.')
hold on;
if exist('plotOpts', 'var')
    q = quiver(Xs,Ys, Us, Vs, 0, {plotOpts});
else
    q = quiver(Xs,Ys, Us, Vs, 0);
end

axis equal
q.ShowArrowHead = 'off';