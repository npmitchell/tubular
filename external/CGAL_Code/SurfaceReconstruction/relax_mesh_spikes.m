function [U, Uall] = relax_mesh_spikes( F, V, a, b, L_method, ...
    fixIDx, lambda, method, max_iter, tol)
% RELAX_MESH_SPIKES Attemps to smooth localized spikes in a mesh
% triangulation utilizing a using implicit/explicit Laplacian smoothing
% method with a geometrically modified Laplacian.
%
%	INPUT PARAMETERS:
%
%       - F:            #Fx3 face connectivity list
%
%       - V:            #Vx3 input vertex coordinate list
%
%       - a:            The minimum dihedral angle that will be fixed in
%                       the smoothing
%
%       - b:            The threshold level at which the vertices will be
%                       uniformly smoothed
%
%       - L_method:     Method for basic Laplacian construction
%                       ('uniform', 'cotan')
%
%       - fixIDx:       #FVx1 list of fixed vertex indices
%
%       - lambda:       diffusion speed parameter
%
%       - method:       Method for smoothing
%                       ('implicitx', 'explicit')
%
%       - max_iter:     The maximum number of smoothing iterations
%
%       - tol:          The threshold for insufficient change used for
%                       smoothing termination
%
%   OUTPUT PARAMETERS:
%
%       - U:            #Vx3 output vertex coordinate list
%
%       - Uall:         #Vx3xiters list of vertex coordinate lists for each
%                       smoothing iteration
%
%   by Dillon Cislo 12/08/2021

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

if (~exist('F', 'var'))
    error('Please supply face connectivity list');
end

if (~exist('V', 'var'))
    error('Please supply vertex coordinate list');
end

validateattributes(V, {'numeric'}, {'2d', 'real', 'finite', 'nonnan'});
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'positive', 'integer', '<=', size(V,1)});

TR = triangulation(F,V);
E = edges(TR);

% #Ex2 array of fIDs of the faces attached to a particular edge.
% If an edge is a border edge (i.e., only attached to a single face), then
% that fID is listed twice for dimensional consistency
resizeCell = @(x) repmat( x, 1, 1+mod(numel(x),2) );
edgeFace = edgeAttachments( TR, E );
edgeFace = cell2mat( cellfun( resizeCell, edgeFace, ...
    'UniformOutput', false ) );

% The average edge length in the mesh
h = avgedge(V,F);

if (~exist('a', 'var')), a = 2.5; end
if (~exist('b', 'var')), b = 3; end
if (~exist('L_method','var')), L_method = 'uniform'; end
if (~exist('fixIDx','var')), fixIDx = []; end
if (~exist('lambda','var')), lambda = 0.1; end
if (~exist('method','var')), method = 'implicit'; end
if (~exist('tol','var')), tol = 0.001; end
if (~exist('max_iter','var')), max_iter = 10; end

%--------------------------------------------------------------------------
% Perform Spike Relaxation
%--------------------------------------------------------------------------

% Only construc the uniform Laplacian once
if strcmpi(L_method, 'uniform')
    A = adjacency_matrix(F);
    L = A - diag(sum(A));
end

% Create a copy of the input vertex coordinate list
S = V;

% Build sparse identity matrix
I = speye(size(V,1), size(V,1));

% Place for factorization and symmtery flag used by min_quad_with_fixed
P = [];
% sym = [];

iter = 0;
U = S;
U_prev = S;
if nargout >= 2
    Uall = [];
end

while( iter < max_iter && (iter == 0 || max(abs(U(:)-U_prev(:)))>tol*h))
    
    U_prev = U; % Store a copy of the previous iteration
    
    % Update geometric Laplacian operator
    if strcmpi(L_method, 'cotan')
        L = cotmatrix_embedded(U,F);
    end
    
    % Update relaxing Laplacian operator
    Rcoeff = relaxing_laplacian_coefficient(U, a, b, TR, E, edgeFace);
    R = repmat(Rcoeff, 1, size(V,1)) .* L;
    
    switch method
        
        case 'implicit'
            % Q = (I-lambda*L);
            Q = (I-lambda*R);
            % could prefactor Q for 'uniform' case
            for d = 1:size(S,2)
                [U(:,d),P] = ...
                    min_quad_with_fixed( Q*0.5, -U(:,d), fixIDx, ...
                    S(fixIDx,d), [], [], P);
            end
            
        case 'explicit'
            % Q = (I+lambda*L);
            Q = (I+lambda*R);
            U = Q * U;
            % enforce boundary
            U(fixIDx,:) = S(fixIDx,:);
            
        otherwise
            error(['' method ' is not a supported smoothing method']);
            
    end
    
    if (nargout >= 2)
        Uall = cat(3, Uall, U);
    end
    
    iter = iter + 1;
    
end

end

function Rcoeff = ...
    relaxing_laplacian_coefficient(V, a, b, TR, E, edgeFace)
% Determine the maximum dihedral angle of the edges attached to each vertex

bdyEdge = ((edgeFace(:,1) - edgeFace(:,2)) == 0);

% Face unit normal vectors
fN = TR.faceNormal;
fN1 = fN(edgeFace(:,1), :);
fN2 = fN(edgeFace(:,2), :);

crossN = cross(fN1, fN2, 2);
crossHat = crossN ./ sqrt(sum(crossN.^2, 2));

crossN(bdyEdge, :) = zeros(sum(bdyEdge), 3);
crossHat(bdyEdge, :) = zeros(sum(bdyEdge), 3);

% The absolute dihedral angles of each edge
EAng = abs(2 .* atan2( dot(crossN, crossHat, 2), 1 + dot(fN1, fN2, 2) ));

% The maximum angle of edges attached to each vertex
Rcoeff = a .* ones(size(V,1),1);
for i = 1:size(E,1)
    if (EAng(i) > Rcoeff(E(i,1))), Rcoeff(E(i,1)) = EAng(i); end
    if (EAng(i) > Rcoeff(E(i,2))), Rcoeff(E(i,2)) = EAng(i); end
end

Rcoeff = (Rcoeff - a) / (b-a);

end
