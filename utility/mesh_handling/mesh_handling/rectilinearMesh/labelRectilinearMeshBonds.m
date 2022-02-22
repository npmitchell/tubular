function [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(mesh, options)
% metricSPhiGridMesh(mesh, varargin)
%   Identify the bond on each face of a rectilinear grid mesh that is along
%   x and y (z and phi). Assumes input mesh is rectilinear with indices
%   increasing first along x, then along y. Mesh may be periodic in either
%   or both directions.
%
% Note: Allows passing topologicalStructureTools instead of recomputing
%
% Parameters
% ----------
% mesh : struct with fields
%   v : N x 3 float
%       mesh vertices
%   f : M x 3 int
%       mesh faces, as defined by indexing into vertices. Assumes
%       consistent face orientation throughout the mesh
%   nU : int
%       number of mesh vertices along each row of the rectilinear mesh
%   nV : int
%       number of mesh vertices along each column of the mesh
% options : struct with fields (optional)
%   eIDx : #bonds x 2 int
%       bond vertex IDs, bond vertices from topologicalStructureTools()
%   feIDx : #faces x 3 int
%       face bond IDs, face bonds from topologicalStructureTools()
%   bulkEdgeIDx : #bonds x 1 int
%       label of whether in bulk (0) or on edge (1)
%   preview : bool
%       show intermediate results
%
% Returns
% -------
% labels : struct with fields
%   fe_is_u : #faces x 3 int array matching feIDx (face-Edge IDs)
%       whether a bond in the triangulation is along u (row)
%   fe_is_v : #faces x 3 int array match feIDx (face-Edge IDs)
%       whether a bond in the triangulation is along v (column)
% dbonds : struct with fields
%   If mesh has no baseSpace (no mesh.u field), then output has fields:
%       u : #faces x dim float array
%           directed bond vector oriented along the 'u' direction of the
%           rectilinear mesh (rows, connecting (i, i+1) vertices)
%       v : #faces x dim float array
%           directed bond vector oriented along the 'v' direction of the
%           rectilinear mesh (columns, connecting (i, i+nU) vertices)
%   Otherwise, if mesh DOES have base space, there is a separate field for
%   each set of directed bond vectors
%       baseSpace : struct with fields u,v as above
%           directed bonds in domain of parameterization
%       realSpace : struct with fields u,v as above
%           directed bonds in coordinate space of mesh vertices
%   If the input mesh is a closed cylinder, cutMesh is an additional field
%       cutMesh : struct with fields
%           f : (nU-1)*2*(nV-1) x 3 int array, face connectivity list
%           u : #vertices x 2 float array, pullback vertex positions
%           v : #vertices x 3 float array, embedding vertex positions
%           vn : #vertices x 3 float array, vertex normals
%           nU : int, number of vertices along pullback U direction
%           nV : int, number of vertices along pullback V direction
%       
% topStructTools : optional cell array output, with elements
%   eIDx : #bonds x 2 int
%       bond vertex IDs, bond vertices from topologicalStructureTools()
%   feIDx : #faces x 3 int
%       face bond IDs, face bonds from topologicalStructureTools()
%   bulkEdgeIDx : #bonds x 1 int
%       label of whether in bulk (0) or on edge (1)
%
% NPMitchell 2020

% Store mesh vertices and faces
VV = mesh.v ;
FF = mesh.f ;
nU = mesh.nU ;

%% Construct Triangulation
tri = triangulation(FF, VV) ;

%% Unpack or compute topological structure tools
load_tstools = true ;
preview = false ;
if nargin > 1
    try
        eIDx = options.eIDx ;
        feIDx = options.feIDx ;
        load_tstools = false ;
    catch
        disp(['labelRectilinearMeshBonds: ', ...
            'topological structure tools could not be parsed'])
    end
    if isfield(options, 'preview')
        preview = options.preview ;
    end
end
if load_tstools
    [eIDx, feIDx, ~] = topologicalStructureTools(tri) ;
    if nargout > 2
        topStructTools = {eIDx, feIDx, } ;
    end
end
% Note: eIDx is the bond list

%% Handle two cases separately: is mesh closed in V or open?
edge = edges( tri );
eulerChi = size(VV, 1) - size(edge,1) + size(FF,1);
if eulerChi == 1 
    mesh_closed = false ;
elseif eulerChi == 0
    mesh_closed = true ;
else
    error('Input rectilinear mesh is neither topological cylinder nor disk. Handle case here')
end

%% Count which bonds connect distant/adjacent/modular indices for 0/1/2
% Initial bonds define (s,phi) coordinate system
% Are they shat (1), phihat (2), or diagonal (0)?
% Build sorp==1 is shat, sorp==2 is phihat, sorp==0 is diagonal
sorp = zeros(size(eIDx(:, 1))) ;
% shat will differ by one
bondDX = abs(eIDx(:, 2) - eIDx(:, 1)) ;
sorp(bondDX == 1) = 1 ;
% phihat will differ by multiple of nU
sorp(mod(bondDX, nU) == 0) = 2 ;

% fe_is_s is boolean, true where feIDx element is along s
fe_is_u = sorp(feIDx) == 1 ;    
% fe_is_phi is boolean, true where feIDx element is along phi
fe_is_v = sorp(feIDx) == 2 ;

% check that there is one uhat vector in each triangle
assert(all(any(fe_is_u, 2)))  
% check that there is one vhat vector in each triangle
assert(all(any(fe_is_v, 2)))  

labels.fe_is_u = fe_is_u ;
labels.fe_is_v = fe_is_v ;

if nargout > 1
    % Directed edge vectors in mesh configuration space
    eij = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);

    % Return directed bond vectors for bonds along u and along v as struct
    % with fields u, v
    dbond_u3d = eij(sum(fe_is_u .* feIDx, 2), :) ;
    dbond_v3d = eij(sum(fe_is_v .* feIDx, 2), :) ;
    
    % Check for a base space of the mesh
    if isfield(mesh, 'u')
        % if mesh is closed, adjust for periodic vertices in base space
        if mesh_closed
            disp('cutting closed input mesh for baseSpace dbonds')
            cutMesh = cutRectilinearCylMesh(mesh) ;
            cutTri = triangulation(cutMesh.f, cutMesh.v) ;
            [eIDx, feIDx, ~] = topologicalStructureTools(cutTri) ;
            UU = cutMesh.u ;
            dbonds.cutMesh = cutMesh ;
        else        
            UU = mesh.u ;
        end
        
        % Directed edge vectors in mesh base space
        eij = UU(eIDx(:,2), :) - UU(eIDx(:,1), :);
        
        % Return directed bond vectors for bonds along u and along v
        dbond_u = eij(sum(fe_is_u .* feIDx, 2), :) ;
        dbond_v = eij(sum(fe_is_v .* feIDx, 2), :) ;
        dbonds_baseSpace.u = dbond_u ;
        dbonds_baseSpace.v = dbond_v ;
        
        % Return dbonds as struct with fields realSpace and baseSpace
        dbonds_realSpace.u = dbond_u3d ;
        dbonds_realSpace.v = dbond_v3d ;
        
        % Store in master struct
        dbonds.realSpace = dbonds_realSpace ;
        dbonds.baseSpace = dbonds_baseSpace ;
    else
        dbonds.u = dbond_u3d ;
        dbonds.v = dbond_v3d ;
    end
    
end

%% preview results if desired
if preview
    clf
    v01 = VV(eIDx(sorp == 0, 1), :) ;
    v02 = VV(eIDx(sorp == 0, 2), :) ;
    plot3([v01(:, 1), v02(:, 1)]', [v01(:, 2), v02(:, 2)]', ...
        [v01(:, 3), v02(:, 3)]', 'k')
    axis equal
    view(2)
    % xlim([0, 100])
    ylim([-10, Inf])
    zlim([0, Inf])
    hold on; 
    pause(1)
    v11 = VV(eIDx(sorp == 1, 1), :) ;
    v12 = VV(eIDx(sorp == 1, 2), :) ;
    plot3([v11(:, 1), v12(:, 1)]', [v11(:, 2), v12(:, 2)]', ...
        [v11(:, 3), v12(:, 3)]', 'r') 
    pause(1)
    v21 = VV(eIDx(sorp == 2, 1), :) ;
    v22 = VV(eIDx(sorp == 2, 2), :) ;
    plot3([v21(:, 1), v22(:, 1)]', [v21(:, 2), v22(:, 2)]', ...
        [v21(:, 3), v22(:, 3)]', 'c')
    
end
