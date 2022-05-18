function [V,F] = remesh_smooth_iterate(V,F, lambda,... 
    tar_length, num_iter, protect_constraints, enforceQuality, maxIterRelaxMeshSpikes)
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
if nargin < 3
    lambda = 0.025 ;
end
if nargin  < 4
    tar_length = 1;
end
if nargin < 5
    num_iter = 5;
end
if nargin < 6
    protect_constraints = false;
end

if ~enforceQuality
    % Simpler case, with smoothing steps
    
    % Attempt to remove localized mesh spikes by Laplacian relaxation
    if maxIterRelaxMeshSpikes > 0 
        V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
            'uniform', [], 2, 'implicit', maxIterRelaxMeshSpikes);
    end
    
    V = laplacian_smooth(V, F, 'cotan', [], lambda, 'implicit', V, 10);
else
    % Isotropically remesh the surface
    try
        [F, V, ~, ~] = isotropic_remeshing( F, V, ...
            tar_length, num_iter, protect_constraints);
    catch
         [F, V, ~, ~] = isotropic_remeshing( F, V, ...
                tar_length, num_iter);
    end

    % Attempt to remove localized mesh spikes by Laplacian relaxation
    V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
        'uniform', [], 2, 'implicit', maxIterRelaxMeshSpikes);

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
            [F, V, ~, ~] = isotropic_remeshing(F, V, ...
                tar_length, num_iter, protect_constraints);

            % Another round of spike relaxation
            V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
                'uniform', [], 2, 'implicit', maxIterRelaxMeshSpikes);

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
        V = laplacian_smooth(V, F, 'cotan', [], lambda, 'implicit', V, 10);
    catch
        disp('Could not smooth additional iteration')
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
    E = edges(triangulation(F, V));

    numBdy = numel(DiscreteRicciFlow.compute_boundaries(F));
    if (numBdy ~= 0)
        error( ['Mesh has %d ' ...
            'boundary components'], t, numBdy );
    end

    eulerChi = size(F,1) + size(V,1) - size(E,1);
    if (eulerChi ~= 2)
        error( ['Mesh is not ' ...
            'a topological sphere'] );
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