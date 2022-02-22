function [intx, specialCase] = linePlaneIntersection(lineDirec, linePt, normalToPlane, ptInPlane)
% Find intersection between a plane and a line.
%
% This solves Descartes' plane equation by parameterizing a line in tt.
%
% Parameters
% ----------
% lineDirec : 3 x 1 float
%       vector pointing along the line (either direction)
% linePt : 3 x 1 float 
%       A point on the line
% normalToPlane : 3 x 1 float
%       normal vector to the plane
% ptInPlane : 3 x 1 float
%       a point in the plane
%
% Returns
% -------
% intx : 3 x 1 float
%       intersection in 3D between line and plane
% specialCase : bool

lineDirec = lineDirec(:) ;
linePt = linePt(:) ;
normalToPlane = normalToPlane(:) ;
ptInPlane = ptInPlane(:) ;

% Find intersection
dd = -dot(normalToPlane, ptInPlane) ;

% Check if line is in plane or parallel to it
if dot(normalToPlane, lineDirec) == 0
    % line is parallel to plane
    if dot(normalToPlane, linePt) == -dd 
        disp('Line is in the plane, intersection is the entire line')
        intx = linePt ;
    else
        disp('Line is parallel to plane without intersection')
        intx = [] ;
    end
    specialCase = true ;
else
    % parameterization of the line
    tt = - (dd + dot(normalToPlane, linePt)) / dot(normalToPlane, lineDirec) ;
    % Intersection is located then at:
    intx = linePt + lineDirec * tt ;
    specialCase = false ;
end

