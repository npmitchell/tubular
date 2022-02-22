function [fb_pieces, hairpins] = separateFreeBoundaryCurves(fb)
%[fb_pieces, hairpins] = separateFreeBoundaryCurves(fb)
%
%   NOTE: This function is a wrapper/alias for
%   curve_functions/separateClosedCurveComponents.m
%   
%   Separate each topologically distinct closed curve of a free boundary 
%   or of vertex indices of curves in space into a separate list of free 
%   boundary vertex ids. 
%   If a closed curve contains a closed curve (ie a "hairpin" component)
%   that is a subset of the full path, then we return this information if 
%   desired in the second output.
%
%   For example, in the image below, there are two closed curves, one of
%   which contains two pieces, so the output is:
%      fb_pieces = {[ijkljmi], [npqn]}
%      hairpins = {{[ijmi],[jklj]}, [npqn]}
%
%            k             p
%           o             o 
%   i    j / \           / \
%   o-----o---o l     n o---o q
%    \   /
%     \ /
%      o m
%
% Parameters
% ----------
% lsegs : #path points x 2 int
%   linesegment indices, so for the image above the input could be 
%   [i,j; j,k; k,l; l,j; j,m; n,p; p,q; pn] or equivalent. 
%
% NPMitchell 2020

[fb_pieces, hairpins] = separateClosedCurveComponents(fb) ;


