function [bdyIDx] = freeBoundary(triangulationObj)
%FREEBOUNDARY Find boundary indices of triangulation object
%   
% Parameters
% ----------
% triangulationObj : MATLAB triangulation object
%   mesh with attributes Points and Connectivity List
% 
% Returns
% -------
% bdyIDx : N x 1 int array
%   the indices of triangulationObj.Points that are boundary vertex points
%
% NPMitchell

bdyIDx = triangulationObj.freeBoundary;

end

