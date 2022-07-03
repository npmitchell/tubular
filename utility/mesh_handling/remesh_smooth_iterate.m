function [V,F] = remesh_smooth_iterate(V, F, options)
% [V,F] = remesh_smooth_iterate(V,F, lambda, tar_length, num_iter, protect_constraints)
% Isotropically remesh the surface and smooth it.
%
% Parameters
% ----------
% V : N x 3 float
%   vertices, locations of vertices in 3d before smoothing/remeshing
% F : M x 3
%   faces, face connectivity list indexing input vertices V
% lambda : float 
%   smoothing parameter for laplacian
% tar_length : numeric
%   target edge length for remeshing
% num_iter : int
%   number of iterations for isotropic remeshing
% protect_constraints : bool
%   protect constraints during isotropic remeshing
% enforceQuality : bool
%   error out if there are boundaries, mesh intersections, or if not
%   spherical topology
% maxIterRelaxMeshSpikes : int >=0 
%   maximum number of iterations of relaxing mesh spikes 
% enforceTopology : bool
%   enforce the topology of the mesh to match targetEulerCharacteristic
% targetEulerCharacteristic : int
%   required topology of the mesh if enforceTopology == true    
%
% Returns
% -------
% V : P x 3 float
%   output vertices, smoothed and remeshed locations in 3d
% F : Q x 3 float 
%   output faces, face connectivity list indexing vertices V
%
% DJ Cislo & NPMitchell 2022

% Input processing
lambda = 0.025 ;
tar_length = 1;
num_iter = 5;
protect_constraints = false;
enforceQuality = true ;
maxIterRelaxMeshSpikes = 10 ;
enforceTopology = false ;
targetEulerCharacteristic = 2 ;

if nargin < 3
    options = struct() ;
end

if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'tar_length')
    tar_length = options.tar_length ;
end
if isfield(options, 'num_iter')
    num_iter = options.num_iter ;
end
if isfield(options, 'protect_constraints')
    protect_constraints = options.protect_constraints ;
end
if isfield(options, 'enforceQuality')
    enforceQuality = options.enforceQuality ;
end
if isfield(options, 'maxIterRelaxMeshSpikes')
    maxIterRelaxMeshSpikes = options.maxIterRelaxMeshSpikes ;
end
if isfield(options, 'enforceTopology')
    enforceTopology = options.enforceTopology ;
end
if isfield(options, 'targetEulerCharacteristic')
    targetEulerCharacteristic = options.targetEulerCharacteristic ;
end

if ~enforceQuality
    % Simpler case, with smoothing steps
    
    % Attempt to remove localized mesh spikes by Laplacian relaxation
    if maxIterRelaxMeshSpikes > 0 
        try
            V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
                'uniform', [], 2, 'implicit', maxIterRelaxMeshSpikes);
        catch
            disp('WARNING: could not run relax_mesh_spikes')
        end
    end
    
    tr = triangulation(F, V) ;
    fb = tr.freeBoundary ;
    if ~isempty(fb)
        fb = fb(:, 1) ;
    end
    V = laplacian_smooth(V, F, 'cotan', fb, lambda, 'implicit', V, 10);
else
    
    % Boundary indices
    tr = triangulation(F, V) ;
    fb = tr.freeBoundary ;
    fb = fb(:, 1) ;
    
    % Isotropically remesh the surface
    try
        [F, V, ~, ~] = isotropic_remeshing( F, V, ...
            tar_length, num_iter, protect_constraints);
    catch
         [F, V, ~, ~] = isotropic_remeshing( F, V, ...
                tar_length, num_iter);
    end

    % Attempt to remove localized mesh spikes by Laplacian relaxation
    try 
        V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
            'uniform', fb, lambda, 'implicit', maxIterRelaxMeshSpikes);
    catch
        disp('WARNING: could not run relax_mesh_spikes')
    end

    % Try to remove self-intersections
    [intersects, ~] = mesh_self_intersection_3d(F, V);
    intCount = 0;
    skipRemesh = false ;
    while intersects
        try
            [V, F] = clean_mesh(V, F, 'MinDist', 0, 'MinArea', 0, ...
                'MinAngle', 0, 'SelfIntersections', 'remove', ...
                'SmallTriangles', 'remove');

            % Peform a another isotropic remeshing
            try
                [F, V, ~, ~] = isotropic_remeshing( F, V, ...
                    tar_length, num_iter, protect_constraints);
            catch
                 [F, V, ~, ~] = isotropic_remeshing( F, V, ...
                        tar_length, num_iter);
            end

            % Another round of spike relaxation
            try
                V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
                    'uniform', fb, lambda, 'implicit', maxIterRelaxMeshSpikes);
            catch
                disp('WARNING: could not run relax_mesh_spikes')
            end
            
            [intersects, ~] = mesh_self_intersection_3d(F, V);

            intCount = intCount + 1;
            if intCount > 20
                error('Unable to remove self-intersections');
            end
        catch
            if enforceQuality
                error('Unable to remove self-intersections')
            else
                disp('Unable to remove self-intersections')
            end
            intersects = false ;
            skipRemesh = true ;
        end

    end

    % Another round of isotropic remeshing
    if ~skipRemesh
        try
            [F, V, ~, ~] = isotropic_remeshing(F, V, ...
                tar_length, num_iter, protect_constraints);
        catch
             [F, V, ~, ~] = isotropic_remeshing( F, V, ...
                tar_length, num_iter);
        end
    end

    % Smooth the entire mesh fixing the boundary
    try
        V = laplacian_smooth(V, F, 'cotan', fb, lambda, 'implicit', V, 10);
    catch
        disp('WARNING: Could not smooth additional iteration')
    end

    % Peform a final isotropic remeshing
    if ~skipRemesh
        try
            [F, V, ~, ~] = isotropic_remeshing(F, V, ...
                tar_length, num_iter, protect_constraints);
        catch
            [F, V, ~, ~] = isotropic_remeshing(F, V, ...
                tar_length, num_iter) ;
        end
    end

    % Mesh quality checks ---------------------------------------------
    if enforceTopology
        E = edges(triangulation(F, V));

        numBdy = numel(DiscreteRicciFlow.compute_boundaries(F));
        if (numBdy ~= 0)
            error( ['Mesh has %d ' ...
                'boundary components'], numBdy );
        end
        
        eulerChi = size(F,1) + size(V,1) - size(E,1);
        if (eulerChi ~= targetEulerCharacteristic)
            error( ['Mesh is not ' ...
                'a topological sphere'] );
        end
    end

    [intersects, intx] = mesh_self_intersection_3d(F, V);
    if intersects
        msg = sprintf( ['Mesh contains %d ' ...
            'self-intersections. '], size(intx, 1)) ;
        msg = [msg 'Try adjusting parameter ' ...
            'tubi.xp.detectOptions.target_edgelength, ' ...
            'reduce lambda, or checking that surface training is valid'] ;
        error(msg );
    end
end