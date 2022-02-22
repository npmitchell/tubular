function trisurf2d(faces, xy, varargin)

if nargin < 3
    trisurf(faces, xy(:, 1), xy(:, 2), 0*xy(:, 1))
else
    trisurf(faces, xy(:, 1), xy(:, 2), 0*xy(:, 1), varargin{:})    
end
view(2)