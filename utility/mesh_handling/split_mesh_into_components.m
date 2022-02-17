function components = split_mesh_into_components( face, ...
    vertex )
%SPLIT_MESH_INTO_COMPONENTS This function splits a composite mesh into its
% isolated connected components for a multiply connected mesh triangulation
%   INPUT PARAMETERS:
%       - face:         NFx3 mesh connectivity list
%       - vertex:       NVxdim mesh vertex coordinate lis
%
%   OUTPUT PARAMETERS:
%       - components:   #components x 1 cell of structs
%                       components{i} contains faces f and vertices v, as
%                       well as oldVertexIDx and newVertexID, which indexes
%                       into the original (composite) mesh.

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% A MATLAB-style representation of the input mesh
tri = triangulation( face, vertex );

% The NEx3 edge connectivity list
edge = tri.edges;

NF = size( face, 1 );
NV = size( vertex, 1 );
NE = size( edge, 1 );

%--------------------------------------------------------------------------
% DETERMINE THE CONNECTED COMPONENTS OF THE INPUT MESH
%--------------------------------------------------------------------------

% Construct an NV x NV vertex adjacency matrix
A = sparse( [ edge(:,1); edge(:,2) ], [ edge(:,2), edge(:,1) ], ...
    ones( 2 * NE, 1 ), NV, NV );

% Construct the Dulmage-Mendelsohn decomposition
[ p, ~, r ] = dmperm( A + speye(size(A)) );

% The number of connected components
S = numel(r)-1;

% A 1xNV row vector. C(v) is the connected component ID of vertex v
C = cumsum( full( sparse( 1, r(1:end-1), 1, 1, size(A,1) ) ) );
C(p) = C;

% Now C tells us which component each vertex belongs to. Use this to split
components = cell(S, 1) ;
% for each connected component, grab the vertices and faces belonging to it
% Initialize a submesh structure, which we will populate for each component
submesh = struct() ;
for ii = 1:S
    rmVIDx = find(C ~= ii) ; 
    % remove everything not in this component
    [ submesh.f, submesh.v, submesh.oldVertexIDx, submesh.newVertexIDx ] = ...
        remove_vertex_from_mesh( face, vertex, rmVIDx ) ;
    components{ii} = submesh ;
end

