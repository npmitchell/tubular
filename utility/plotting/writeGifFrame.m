function writeGifFrame(fig, fname, isFirst)
% writeGifFrame(fig, fname, isFirst)
% Saving gif file frame by frame
% 
% Parameters
% ----------
% fig : figure handle
% fname : char, filename for the gif
% isFirst : bool, whether this frame is the first one
% 
    [imind, cm] = rgb2ind(frame2im(getframe(fig)), 256);
    if isFirst
        imwrite(imind, cm, fname, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
    else
        imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
    end
end
