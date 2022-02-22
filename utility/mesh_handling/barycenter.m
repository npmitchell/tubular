function [ BC ] = barycenter(VV, FF)
% BARYCENTER(V,F) Compute the barycenter of triangular faces in a mesh
%
% BC = barycenter(VV,FF)
%
% Parameters
% ----------
% VV : #Vertices x D numeric
%   vertex coordinates
% FF : #Faces x simplex_size integer matrix
%   face connectivity list (indices of triangle corners if simplex_size=3)
% Returns
% -------
% BC : #Faces x D float 
%   barycenters of all faces as 3d positions
% 
% This function shadows gptoolbox's barycenter function, in an effort to
% remove dependencies. The principle of calculation is exactly the same. 
%

BC = zeros(size(FF,1),size(VV,2));
for ii = 1:size(FF,2)
    BC = BC + 1/size(FF,2) * VV(FF(:,ii),:);
end

end

