function sf = interpolateOntoPullbackXY(QS, XY, scalar_field, options)
% interpolateOntoPullbackXY(QS, XY, scalar_field, options)    
%   Interpolate scalar field defined on current mesh onto 
%   supplied pullback coordinates. Input scalar_field may be
%   defined on vertices or faces. Automatically detects whether
%   field is at 1x mesh resolution or higher. The pixel locations of the
%   pullback are computed from scratch via options.Lx and options.Ly, not
%   passed directly.
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   
% XY : Nx2 numeric array
%   positions at which to evaluate scalar field
% scalar_field : Nx1 float array
%   the field to interpolate, defined on either vertices or faces of the
%   current mesh 
% options : struct with fields
%   imCoords : str pullback specifier (default=QS.piv.imCoords)
%       coordinate system of the pullback space
%   sfLocation : 'vertices' or 'faces' (default='vertices')
%       evaluation location of the scalar field on currentMesh
%   iMethod : 'linear', 'cubic', or 'nearest' (default='linear')
%       interpolation method for scalar field
%   Lx : int or float (optional, saves time if supplied)
%       extent of pullback image space in horizontal dimension
%   Ly : int or float (optional, saves time if supplied)
%       extent of pullback image space in vertical dimension
%   
% Returns
% -------
% sf : Nx1 float array
%   scalar field evaluated onto pullback coordinates

% default options
coordSys = QS.piv.imCoords;  % coodinate system for pullback, default is (s,phi) smoothed extended
sfLocation = 'vertices' ;    % whether scalar field is on vertices or faces
iMethod = 'linear' ;         % interpolation method for scalar field
% Unpack options
if isfield(options, 'imCoords')
    coordSys = options.coordSys ;
end
% pullback image is a double cover of the mesh (ie 'extended')
if strcmp(coordSys(end), 'e')
    doubleCovered = true ;
end
if isfield(options, 'Lx') && isfield(options, 'Ly')
    Lx = options.Lx ;
    Ly = options.Ly ;
else
    im0 = QS.loadCurrentPullback(coordSys) ;
    [Lx, Ly] = size(im0) ;
end
nU = QS.nU ;
nV = QS.nV ;

% Obtain coordinates from which to interpolate
if strcmp(coordSys, 'sp_sme')
    % mesh = QS.getCurrentSPCutMeshSm() ;                
    % tileCount = [1, 1] ;
    % [TF, TV2D] = tileAnnularCutMesh(mesh, tileCount) ;
    % Mesh uv is in rectilinear grid. Create tiled rectilinear
    % grid with properly increasing Y coordinates. At present
    % 2020-07, the tileAnnulurCutMesh gives the wrong ordering
    % in Y dimension. Assumes that mesh vertices have
    % increasing y values from row i to row i+1.
    uspace = linspace(0, 1, nU) ;
    if doubleCovered
        vspace = linspace(0, 1, 3*nV-2) ;
    else
        vspace = linspace(0, 1, nV) ;
    end
    [TV2Du, TV2Dv] = ndgrid(uspace, vspace) ;
    TV2D = [TV2Du(:), TV2Dv(:)];
    if doubleCovered
        mXY = QS.uv2XY([Lx, 1.5 * Ly], TV2D, false, 1, 1) ;
        mXY(:, 2) = mXY(:, 2) - 0.25 * Ly ;
    else
        error('Make 3 tiles of XY, so should span from approx -Ly to 2Ly.')
        mXY = QS.uv2XY([Lx, 3 * Ly], TV2D, false, 1, 1) ;
    end
else
    error(['handle this coordSys here: ' coordSys])
end

if strcmp(sfLocation, 'vertices')
    % Check that resolution of supplied sf is 1x
    if length(scalar_field(:)) == nU * nV
        sfe = reshape(scalar_field, [nU, nV]) ;
        % Assumes that mesh vertices are indexed as 
        % row1, row2, ... row_nV, with each row containing nU
        % vertices.
        % Extend the scalar field with above assumption in Y
        sfe = cat(2, sfe(:, 1:nV-1), sfe, sfe(:, 2:end)) ;
        % check sizes
        assert(length(sfe(:)) == size(TV2D, 1))
        % reshape into tiled grid [+1,-1]
        xx = reshape(mXY(:, 1), [nU, 3*nV-2]) ;
        yy = reshape(mXY(:, 2), [nU, 3*nV-2]) ;
        Fsf = griddedInterpolant(xx, yy, sfe,...
            iMethod, 'nearest') ;
        sf = Fsf(XY(:, 1), XY(:, 2)) ;
    else
        disp('Size of scalar field = ')
        disp(size(scalar_field))
        disp(['nU*nV = ', num2str(nU*nV)])
        error('Handle higher resolution scalar_field here')
    end
elseif strcmp(sfLocation, 'faces')
    % Scalarfield is on barycentric coordinates of faces
    bcXY = barycenter(mXY, mesh.f) ;
    Fsf = scatteredInterpolant(bcXY(:, 1), bcXY(:, 2), ...
                scalar_field', iMethod, 'nearest') ;
    sf = Fsf(XY(:, 1), XY(:, 2)) ;
    error('check this case')
else
    error(['Did not recognize sfLocation: ' sfLocation])
end
