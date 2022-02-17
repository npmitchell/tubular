function wN = annularPathWindingNumber(F, V, P)
%ANNULARPATHWINDINGNUMBER Calculates the winding number of a path along the
%edges of a mesh triangulation with annular topology.  It is assumed that
%the path begins on the inner boundary and terminates on the outer
%boundary.
%
%   INPUT PARAMETERS:
%       - F:    #Fx3 face connectivity list
%       - V:    #Vx2 vertex coordinate list
%       - P:    Nx1 path list of vertex IDs
%
%   OUTPUT PARAMETERS:
%       -wN:    The winding number of the path
%
% by Dillon Cislo 11/21/2019

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

% Verify input triangulation attributes -----------------------------------
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes(V, {'numeric'}, ...
    {'2d', 'ncols', 2', 'finite', 'nonnan'} );

% Verify input mesh topology ----------------------------------------------

% A MATLAB-style triangulation representation
TR = triangulation(F, V);

% The edge connectivity list
E = TR.edges;

% The boundaries of the input trianglation
bdy = DiscreteRicciFlow.compute_boundaries(F);

if ( ( size(V,1) - size(E,1) + size(F,1) ) ~= 0 ) || (numel(bdy) ~= 2)
    error('Input triangulation is not a topological annulus');
end

% Verify input mesh boundary properties -----------------------------------

% Determine which boundary is the outer boundary
if all(inpolygon(V(bdy{2},1), V(bdy{2},2), V(bdy{1},1), V(bdy{1},2)))
    outerIDx = bdy{1};
    innerIDx = bdy{2};
elseif all(inpolygon(V(bdy{1},1), V(bdy{1},2), V(bdy{2},1), V(bdy{2},2)))
    outerIDx = bdy{2};
    innerIDx = bdy{1};
else
    error('Invalid mesh boundaries supplied');
end

% Check that the inner boundary has (0,0) in its interior
if ~inpolygon(0, 0, V(innerIDx,1), V(innerIDx,2))
    
    % Try scaling the outer boundary to the unit disk
    V = ( V - min(V, [], 1) );
    V = V ./ max(V);
    V = 2 .* ( V - 0.5 );
    
    % Then shift COM of the inner boundary to zero
    COM = mean( V(innerIDx,:), 1 );
    V = V - COM;
    
    % Still possible that the inner boundary does not contain (0,0) if it
    % is a non-convex polygon
    if ~inpolygon(0, 0, V(innerIDx,1), V(innerIDx,2))
        error(['Inner boundary must contain the point ' ...
            '(0,0) within its interior']);
    end
    
end

% Verify path properties --------------------------------------------------
validateattributes(P, {'numeric'}, ...
{'2d', 'vector', 'integer', 'positive'});

% Ensure that P is a column vector
if (size(P,2) ~= 1), P = P.'; end

% Verify the path origin and termination ----------------------------------
if ~ismember(P(1), innerIDx)
    P = flipud(P); % Try reversing the path
    if ~ismember(P(1), innerIDx)
        error('Path must originate on the inner boundary');
    end
end

if ~ismember(P(end), outerIDx)
    error('Path must terminate on the outer boundary');
end

% Verify that the path runs along mesh edges ------------------------------
PE = [P(1:(end-1)), P(2:end)];
if any( ~ismember(sort(PE,2), sort(E,2), 'rows') )
    error('Path must run along mesh edges');
end

% Verify that the path does not contain any boundary edges ----------------
bdyE = [ innerIDx', circshift(innerIDx, -1)'; ...
    outerIDx', circshift(outerIDx, -1)' ];
bdyE = sort(bdyE, 2);
if any( ismember( sort(PE,2), bdyE, 'rows' ) )
    error('Path cannot contain any boundary edges');
end

% Rotate the origin point of the mesh to lie on the real axis -------------
VC = complex( V(:,1), V(:,2) );
VC = exp(-1i * angle(VC(P(1)))) .* VC;
V = [real(VC) imag(VC)];

%--------------------------------------------------------------------------
% Calculate the Winding Number
%--------------------------------------------------------------------------
% The winding number is determined by the number of signed crossings of the
% positive real axis

% Do not consider the first segment of the path
PE = PE(2:end, :);

% Calculate the number of positive crossings
numPos = (V(PE(:,1),2) <= 0) & (V(PE(:,2),2) >= 0) & (V(PE(:,2),1) >= 0);
numPos = sum(numPos);

% Calculate the number of negative crossings
numNeg = (V(PE(:,1),2) >= 0) & (V(PE(:,2),2) <= 0) & (V(PE(:,2),1) >= 0);
numNeg = sum(numNeg);

wN = numPos - numNeg;
wN = wN - sign(wN); % A single crossing corresponds to wN = 0;

end

