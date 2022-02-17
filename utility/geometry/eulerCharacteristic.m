function eulerChar = eulerCharacteristic(mesh)
% eulerChar = eulerCharacteristic(mesh)
% 
% Euler characteristic (V-E+F = 2-2*genus) of a mesh with fields f and v 
% for face connectivity and vertex positions in N dimensions, respectively.
% 
% For sphere, eulerChar == 2
% For annulus, eulerChar == 0
% For disk, eulerChar == 1
% For torus, eulerChar == 0
% 
% Parameters
% ----------
% mesh : struct with fields
%   f : (#faces x 3 int array) face connectivity list, as indices into v
%   v : (#vertices x N numeric array) vertex positions in N dimensions
% 
% Returns
% -------
% eulerChar : int
%   the euler characteristic of the mesh
% 
% NPMitchell 2020


    % The #Ex2 edge connectivity list of the mesh
    edgeTri = edges( triangulation( mesh.f, mesh.v ) );
    % Check that the input mesh is a topological cylinderChar
    eulerChar = ( size(mesh.v,1) - length(edgeTri) + size(mesh.f,1) ) ;
end