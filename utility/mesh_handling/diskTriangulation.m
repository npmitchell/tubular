function diskTri_Final = diskTriangulation(N,R)
%DISKTRIANGULATION Generates a triangulation of a disk of specified radius
%   Input Parameters:
%       - N:            - A measure of the triangulation resolution.  See
%                         below for more details (default = 16)
%       - R:            - The radius of the disk (default = 1)
%
%   Output Parameters:
%       - diskTri_Final:      - The final triangulation
%
%   By Dillon Cislo

%==========================================================================
% PROCESS USER INPUT
%==========================================================================
if (nargin < 1), N = 16; end
if (nargin < 2), R = 1; end

%==========================================================================
% CONSTRUCT UNREFINED DISK TRIANGULATION
%==========================================================================
% First we create a hexagon inscribed in the unit disk tiled by equilateral
% triangles of side length 1/N
%
% This will create a triangulation with: 3*N^2 + 3*N + 1 vertices
%                                              and 6*N^2 faces
%
% The vertices of the hexagonal triangulation are then stretched radially
% to fit the unit disk
%==========================================================================
Nv = 3*N^2 + 3*N + 1;
vertex = zeros(Nv,1); % The vertex list

% The seed point is a corner vertex of the inscribed hexagon
seed = exp( 1i * ( 2*pi/3+pi/6 ) );

% A fixed translation vector along the edge of a triangular face
transVec = (1/N) .* exp( 1i * pi / 6 );

%--------------------------------------------------------------------------
% Iterate over vertex columns to construct the triangulation vertex list
%--------------------------------------------------------------------------
count = 0;
for col = 1:(N+1)
    
    vertexCol = ( seed + (col-1) * transVec ) .* ones(N+col,1) + ...
        -1i * (1/N) .* (0:(N+col-1))';
    vertexIDx = count + (1:(N+col))';
    
    % Create vertices on the 'left half' of the hexagon -------------------
    vertex( vertexIDx ) = vertexCol;
    
    % Reflect vertices across the imaginary axis to create vertices on the
    % 'right half' of the hexagon -----------------------------------------
    if col < (N+1)
        vertex( flipud( (Nv+1) - vertexIDx ) ) = -conj( vertexCol );
    end
    
    count = count + N + col;
    
end

% Rotate vertices ---------------------------------------------------------
vertex = exp(-1i*(pi/3+pi/6)) .* vertex;

%--------------------------------------------------------------------------
% Construct counter-clockwise oriented face connectivity list
%--------------------------------------------------------------------------

% Construct right-pointing triangles --------------------------------------
RF = (1:(Nv-N-1))';

tmp = [N+(2:(N+1))'; flipud(N+(2:(N+1))')];

rmIDx = tmp; rmIDx(N) = [];
rmIDx = (N+1) .* ones(2*N,1) + [0; cumsum(rmIDx)];

RF(rmIDx) = [];

tmp2 = [N+(0:N)'; flipud(N+(0:(N-1))')]; tmp2(end) = [];

RF = [ RF, RF + 1, RF + repelem( tmp, tmp2 ) ];

% Construct left-pointing triangles ---------------------------------------
LF = ((N+2):Nv)';
LF( rmIDx - N ) = [];

LF = [ LF, LF - 1, LF - repelem( tmp, flipud(tmp2) ) ];

% Combine results ---------------------------------------------------------
face = [ RF; LF ];

%--------------------------------------------------------------------------
% Radially stretch the hexagon to lie on the unit disk
%--------------------------------------------------------------------------
vertex = 2/sqrt(3) .* abs(vertex) .* ...
    sin( 2 * pi/3 - mod(wrapTo2Pi(angle(vertex)),pi/3) ) .* ...
    exp( 1i * angle(vertex) );

%--------------------------------------------------------------------------
% Construct unrefined disk triangulation
%--------------------------------------------------------------------------
diskTri = triangulation( face, [ real(vertex), imag(vertex) ] );

%==========================================================================
% REFINE TRIANGULATION STEP 1
%==========================================================================
% Remove the corner vertices and combine the two faces attached to the
% corner vertices into a single face comprised of the remaning vertex
% neighbors
%
% The final triangulation will have: 3*N^2 + 3*N -5 vertices
%                                     and 6*(N^2-1) faces
%==========================================================================

%--------------------------------------------------------------------------
% Determine vIDs of the six corner vertices
%--------------------------------------------------------------------------
M = round( N * ( 3 * N + 1) / 2 );
cornerIDx = [ 1; N+1; M+1; M+2*N+1; Nv - N; Nv ];

%--------------------------------------------------------------------------
% Determine the fIDs of the 12 faces attached to the corner vertices
%--------------------------------------------------------------------------
cornerFaceIDx = cell2mat( diskTri.vertexAttachments( cornerIDx ) )';
cornerFaceIDx = cornerFaceIDx(:);

%--------------------------------------------------------------------------
% Construct the six new faces which will replace the old corner faces
%--------------------------------------------------------------------------
newFaces = face( cornerFaceIDx, : )';
newFaces = reshape( newFaces, 6, 6 )';

newFaces = cell2mat( cellfun( @unique,  mat2cell(newFaces, ones(6,1)), ...
    'UniformOutput', false ) )';

newFaces = newFaces(:);
newFaces( ismember( newFaces, cornerIDx ) ) = [];
newFaces = reshape( newFaces, 3, 6 )';

%--------------------------------------------------------------------------
% Sort the new faces to be counter-clockwise oriented
%--------------------------------------------------------------------------
for i = 1:size(newFaces,1)
    
    % Shift the centroid of the new face to the origin --------------------
    curVertex  = vertex( newFaces(i,:) );
    curCent = mean( curVertex );
    curVertex = curVertex - curCent;
    
    % The phase angle of each of those vertices defining a face -----------
    theta = wrapTo2Pi( angle( curVertex ) );
    
    % Sort the indices ----------------------------------------------------
    [~, I] = sort( theta );
    newFaces(i,:) = newFaces(i,I);

end

%--------------------------------------------------------------------------
% Compute updated triangulation
%--------------------------------------------------------------------------

% Create invalid triangulation --------------------------------------------
face = [ face; newFaces ];
diskTri = triangulation( face, [ real(vertex), imag(vertex) ] );

% Remove invalid faces and corner vertices --------------------------------
diskTri = remove_vertex( cornerIDx, diskTri );

%==========================================================================
% REFINE TRIANGULATION STEP 2
%==========================================================================
% Redistribute the remaining boundary vertices uniformly against the
% midpoints of the previous radial level.  This produces a more regular
% triangulation at the edge, which can help with numerical accuracy
%==========================================================================

%--------------------------------------------------------------------------
% Update vertex/connectivity lists
%--------------------------------------------------------------------------
face = diskTri.ConnectivityList;
vertex = diskTri.Points;

%--------------------------------------------------------------------------
% Determine the vIDs of boundary vertices
%--------------------------------------------------------------------------
bdyIDx = unique( diskTri.freeBoundary, 'stable' );

%--------------------------------------------------------------------------
% Determine the fIDs of the three faces attached to each boundary vertex
%--------------------------------------------------------------------------
bdyFaceIDx = cell2mat( diskTri.vertexAttachments( bdyIDx ) )';
bdyFaceIDx = bdyFaceIDx(:);

%--------------------------------------------------------------------------
% Extract the vIDs defining those boundary faces
%--------------------------------------------------------------------------
bdyFaceVIDx = diskTri.ConnectivityList( bdyFaceIDx, : );

%--------------------------------------------------------------------------
% Remove all faces from this list with more than two boundary vertices
%--------------------------------------------------------------------------
rmIDx = sum( ismember(bdyFaceVIDx, bdyIDx), 2 );
rmIDx = rmIDx == 2;

bdyFaceVIDx( rmIDx, : ) = [];

%--------------------------------------------------------------------------
% Remove all boundary vIDs from the remaining vertex list to find the
% opposing edge associated with each boundary vertex
%--------------------------------------------------------------------------
bdyFaceVIDx = bdyFaceVIDx';
bdyFaceVIDx = bdyFaceVIDx(:);

bdyFaceVIDx( ismember(bdyFaceVIDx, bdyIDx) ) = [];
bdyFaceVIDx = reshape( bdyFaceVIDx, 2, size(bdyIDx,1) )';

%--------------------------------------------------------------------------
% Find the midpoints of the associated edges
%--------------------------------------------------------------------------
MP = vertex( bdyFaceVIDx(:,2), : ) - vertex( bdyFaceVIDx(:,1), : );
MP = vertex( bdyFaceVIDx(:,1), : ) + MP ./ 2;

%--------------------------------------------------------------------------
% Reposition the boundary vertices
%--------------------------------------------------------------------------
MP = complex( MP(:,1), MP(:,2) );

vertex = complex( vertex(:,1), vertex(:,2) );
vertex( bdyIDx ) = exp( 1i * angle(MP) );

%==========================================================================
% FINAL OUTPUT PROCESSING
%==========================================================================

%--------------------------------------------------------------------------
% Resize unit disk to specified radius if necessary
%--------------------------------------------------------------------------
if R ~= 1
    vertex = R * abs(vertex) * exp( 1i * angle(vertex) );
end

%--------------------------------------------------------------------------
% Compile output
%--------------------------------------------------------------------------
diskTri_Final = triangulation( face, [ real(vertex), imag(vertex) ] );

end


%**************************************************************************
% REMOVE_VERTEX
%**************************************************************************

function [ trNew ] = remove_vertex( vIDx, tr )
%REMOVE_VERTEX Removes specific vertices and any faces containing those
%vertices from a triangulation.
%   Input Parameters:
%       - vIDx:             The list of vertex IDs to be removed
%       - tr:               A triangulation object
%
%   Output Parameters:
%       -trNew:             The updated triangluation

% Remove the vertices at the specified index values
newVertices = tr.Points;
newVertices(vIDx,:) = [];

% Find the new index for each of the remaining vertices
[~, newVertexIndex] = ismember(tr.Points, newVertices, 'rows');

% Find any faces that used the old vertices and remove them
newFaceList = tr.ConnectivityList;
row_v = false(size(newFaceList,1),1);

for i = 1:length(vIDx)
    row_v = row_v | any( newFaceList == vIDx(i), 2 );
end

newFaceList(row_v, :) = [];

% Now update the vertex indices to the new ones
newFaceList = newVertexIndex(newFaceList);

% Remove any unused vertices in the new triangulation
nullVIDx = find( ~ismember( 1:size(newVertices,1), newFaceList ) );
if ~isempty(nullVIDx)
    
    tempVtx = newVertices;
    newVertices(nullVIDx,:) = [];
    [~, newVertexIndex] = ismember(tempVtx, newVertices, 'rows');
    
    newFaceList = newVertexIndex(newFaceList);
    
end

trNew = triangulation(newFaceList, newVertices);


end

