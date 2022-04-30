function [ newFace, newVertex, oldVertexIDx, C ] = remove_isolated_mesh_components( face, ...
    vertex, varargin )
%REMOVE_ISOLATED_MESH_COMPONENTS This function removes isolated connected
%components of a multiply connected mesh triangulation
%   INPUT PARAMETERS:
%       - face:         NFx3 mesh connectivity list
%       - vertex:       NVxdim mesh vertex coordinate list
%       - vThresh:      The minimum number of vertices below which
%                       a connected component is removed
%
%   OUTPUT PARAMETERS:
%       - newFace:      NF'x3 output mesh connectivity list
%       - newVertex:    NV'xdim output mesh vertex coordinates list
%       - oldxVertexIDx: NV'x1 list of the vertex IDs of the remaining
%                       vertices in the old vertex list
%       - C:            The number of connected components in the input
%                       mesh

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Set the vertex count threshold
vThresh = -1;
if nargin > 2
    vThresh = varargin{1};
end

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

%--------------------------------------------------------------------------
% CLEAN UP CONNECTED COMPONENTS
%--------------------------------------------------------------------------

if S == 1
    
    % disp( [ 'Mesh contains a single connected component. ', ...
    %     'Returning input mesh.' ] );
    
    newFace = face;
    newVertex = vertex;
    oldVertexIDx = (1:size(vertex,1))';
    C = 1;
    
    return;
    
end

% Find the vIDs of the vertices to remove ---------------------------------

% A 1xS row vector. vC(i) is the number of vertices in the ith connected
% component
vC = histc( C, 1:S );


% If vThresh < 0, remove all components except the largest
% WARNING: If the two largest mesh components have the same number of
% vertices than this method removes one arbitrarily so that only one
% component remains
if vThresh < 0
    
    % The component ID of the largest component
    [ ~, bigC ] = max( vC );
    
    rmVIDx = find( C ~= bigC );
    
else
    
    % The component IDs of the components whose size lies below the
    % threshold
    rmC = find( vC <= vThresh );
    
    rmVIDx = false( size(C) );
    for i = 1:numel(rmC)
        rmVIDx = rmVIDx | ( C == rmC(i) );
    end
    rmVIDx = find( rmVIDx );
    
end

% Remove indicated components and assemble output mesh structure ----------

% Remove the lone vertices at the specified index values
newVertex = vertex;
newVertex(rmVIDx, :) = [];

% Find the new index for each of the new vertices
[~, newVertexIDx] = ismember( vertex, newVertex, 'rows');

% Find any faces that contained the indicated vertices and remove them
newFace = face;

loc_v = false(size(face));
for i = 1:numel(rmVIDx)
    loc_v = loc_v | ( newFace == rmVIDx(i) );
end
loc_v = find( loc_v );

% [~, loc_v] = ismember(rmVIDx, newFace);
[row_v, ~] = ind2sub(size(newFace), loc_v);
row_v = unique(row_v);
newFace(row_v, :) = [];

% Now update the vertex indices in the connectivity list
newFace = newVertexIDx(newFace);

% Find the vIDs of the new vertices in the old vertex list
[ ~, oldVertexIDx] = ismember( newVertex, vertex, 'rows' );

end


