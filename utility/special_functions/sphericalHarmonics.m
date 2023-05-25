function YLM = sphericalHarmonics(P, L, M, SHType)
%SPHERICALHARMONICS A vectorized calculation of the spherical harmonics
%evaluated at a set of points on the unit sphere
%
%   INPUT PARAMETERS:
%
%       - P:        #Px3 point coordinate list
%       - L:        #Nx1 vector of spherical harmonic l-indices
%       - M:        #Nx1 vector of spherical harmonic m-indices
%       - SHType:   Real or imaginary convention
%
%   OUTPUT PARAMETERS:
%
%       - YLM:      #PxN array of spherical harmonic valued for each set of
%                   (l,m)-indicies
%
%   by Dillon Cislo 04/21/2023

if (nargin < 4), SHType = 'real'; end

validateattributes(P, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real'});
validateattributes(L, {'numeric'}, ...
    {'vector', 'nonnegative', 'integer', 'finite', 'real'});
validateattributes(M, {'numeric'}, ...
    {'vector', 'integer', 'finite', 'real'});
validateattributes(SHType, {'char'}, {'vector'});

SHType = lower(SHType);
assert(ismember(SHType, {'real', 'imaginary'}), ['Invalid spherical ' ...
    'harmonic convention. Please choose real or imaginary']);

if (size(L,2) ~= 1), L = L.'; end
if (size(M,2) ~= 1), M = M.'; end
assert( all(abs(M) <= L), 'Please enforce l >= |m| for all (l,m)');

% Normalize the input point to lie on the unit sphere for good measure
numPts = size(P,1);
P = P ./ repmat(sqrt(sum(P.^2, 2)), 1, 3);
phi = atan2(P(:,2), P(:,1));

% Evaluate the spherical harmonics
YLM = zeros(numPts, numel(L));

if strcmpi(SHType, 'real')
    
    for i = 1:numel(L)
        
        l = L(i); m = M(i);
        
        % Calculates Plm for abs(m)
        Plm = legendre( l, P(:,3) );
        Plm = reshape( Plm(abs(m)+1,:), [numPts, 1] );
        
        a = (2*l+1) * factorial(l-abs(m));
        b = 4 * pi * factorial(l+abs(m));
        C = (-1)^m * sqrt(a/b);
        
        if m < 0
            YLM(:,i) = sqrt(2) * C .* Plm .* sin( abs(m) .* phi );
        elseif m > 0
            YLM(:,i) = sqrt(2) * C .* Plm .* cos( m .* phi );
        else
            YLM(:,i) = C .* Plm;
        end
        
    end
    
else
    
    for i = 1:numel(L)
        
        l = L(i); m = M(i);
        
        % Calculates Plm for abs(m)
        Plm = legendre( l, P(:,3) );
        Plm = reshape( Plm(abs(m)+1,:), [numPts, 1] );
        
        a = (2*l+1) * factorial(l-abs(m));
        b = 4 * pi * factorial(l+abs(m));
        C = (-1)^m * sqrt(a/b);
        
        YLM(:,i) = C .* Plm .* exp(1i * m * phi);
        
    end
    
end

end

