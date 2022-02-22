function [ P, tri ] = generateSphereMesh( generation, type )
% generateSphereMesh(generation, type)
% Generate a nearly spherical mesh by bisecting a platonic solid
% generation times. 
% 
% Parameters
% ----------
% generation : int
% type : string specifier ('oct', 'tet', 'ico', default='ico')
%   type of Platonic solid to begin with for bisecting
% 
% Returns
% -------
% P : Nx3 float array
%   vertices of approx spherical mesh
% tri : #faces x 3 int array
%   connectivity list indexing into P for each face
% 
% Commented by NPM 2020

    if nargin < 2
        type = 'ico' ;
    end

	if strcmp( type, 'tet' ) % starts with tetrahedron
		[ P, tri ] = getTetrahedralMesh();
	elseif strcmp( type, 'oct' ) % starts with octahedron
		[ P, tri ] = getOctahedralMesh();
	elseif strcmp( type, 'ico' ) % starts with icosahedron
		[ P, tri ] = getIcosahedralMesh();
	else
		fprintf( 'Error: Unrecognized polyhedron type. \n' );
		P = eye( 3 );
		tri = [ 1 2 3 ];
	end

	for n = 1:generation
		[ P, tri ] = refineMesh( P, tri );
	end
	
end
