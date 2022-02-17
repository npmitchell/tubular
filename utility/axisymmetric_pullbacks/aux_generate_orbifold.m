function aux_generate_orbifold(cutMesh, a, IV, imfn, Options, axisorder, save_as_stack)
%AUX_GENERATE_ORBIFOLD(cutMesh, a, IV, imfn)
% called by QuapSlap.generateCurrentPullbacks()
%
% Parameters
% ----------
% cutMesh : struct
%   mesh with fields f (faces), u (2d vertices), and v (3d vertices)
% a : float
%   aspect ratio of width/height of image pullback
% IV : 
%   3d intensity data
% imfn : str
%   path to filename for saving pullback
% Options:  struct
%   Structure containing the standard options for a
%   textured surface patch, such as EdgeColor, EdgeAlpha,
%   etc.  See MATLAB documentation for more information.
%   Additional options as fields are
%       - Options.imSize:       The size of the output image
%       - Options.baseSize:     The side length in pixels of the smallest
%                               side of the output image when using the
%                               tight mesh bounding box
%       - Options.xLim:         The x-bounds of the output image
%       - Options.yLim:         The y-bounds of the output image
%       - Options.pixelSearch:  The method of searching for the faces
%                               containing pixel centers
%                                   - 'AABB' (requires GPToolBox)
%                                   - 'Default' (MATLAB built-ins, faster than AABB)
%       - Options.numLayers:    The number of onion layers to create
%                               Format is [ (num +), (num -) ]
%       - Options.layerSpacing: The spacing between adjacent onion layers
%                               in units of pixels
%       - Options.preSmoothIter:   Number of iterations of Laplacian mesh
%                               smoothing to run on the mesh prior to
%                               vertex normal displacement (requires
%                               GPToolBox) (Default is 0)
%       - Options.smoothIter:   Number of iterations of Laplacian mesh
%                               smoothing to run on the mesh prior to
%                               vertex normal displacement (requires
%                               GPToolBox) (Default is 0)
%       - Options.vertexNormal: User supplied vertex unit normals to the
%                               texture triangulation
%       - Options.Interpolant:  A pre-made texture image volume interpolant
%
% NPMitchell 2020, wrapper for texture_patch_to_image() -- credit to Dillon
% Cislo for building the concepts in this code and for edits.

% Generate Tiled Orbifold Triangulation ------------------------------
tileCount = [1 1];  % how many above, how many below
[ TF, TV2D, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount );

if isfield(Options, 'preSmoothIter')
    preSmoothIter = Options.preSmoothIter ;
else
    preSmoothIter = 0 ;
end

% Before making pullback, smooth surface
if preSmoothIter > 0
    disp(['Smoothing mesh before passing to texture_patch_to_image(): iter=' num2str(preSmoothIter)])
    % Hold the positions of the boundary vertices fixed under smoothing
    % to avoid mesh collapse
    realTri = triangulation(TF, TV3D) ;
    bdyIDx = unique(realTri.freeBoundary);
    
    % check it
    % trisurf(TF, TV3D(:, 1), TV3D(:, 2), TV3D(:, 3), TV3D(:, 2), 'edgecolor', 'none')
    % hold on;
    TV3D = laplacian_smooth( TV3D, TF, 'cotan', ...
        bdyIDx, 0.1, 'implicit', TV3D, preSmoothIter );
    
    % Check smoothed mesh on top
    % trisurf(TF, TV3D(:, 1), TV3D(:, 2), TV3D(:, 3), TV3D(:, 2), 'facecolor', 'none')
    % waitfor(gcf)
end

% View Results -------------------------------------------------------
% patch( 'Faces', TF, 'Vertices', TV2D, 'FaceVertexCData', ...
%     TV3D(:,3), 'FaceColor', 'interp', 'EdgeColor', 'k' );
% axis equal

if nargin < 5
    Options = struct ;
end
if nargin < 6
    axisorder = [1 2 3] ;
end
if nargin < 7
    % Default option for whether to save as a stack
    if isfield(Options, 'numLayers')
        if iscell(Options.numLayers)
            save_as_stack = true ;
        elseif ~isempty(Options.numLayers)
            if any(Options.numLayers)
                save_as_stack = true ;
            else
                save_as_stack = false ;
            end
        else
            save_as_stack = false ; 
        end
    else
        save_as_stack = false ;
    end
end

% Texture image options
if ~isfield(Options, 'PSize')
    Options.PSize = 5;
end
if ~isfield(Options, 'EdgeColor')
    Options.EdgeColor = 'none';
end
if isfield(Options, 'imSize')
    if length(Options.imSize) > 1
        disp("WARNING: Options.imSize is overwritten by parameter 'a'")
    else
        imsz = Options.imSize ;
    end
else
    imsz = 1000 ;
end
if ~isfield(Options, 'yLim')
    Options.yLim = [0 1];
end

Options.imSize = ceil( imsz .* [ 1 a ] ) ;

% profile on
% Create texture image
if any(isnan(TV2D))
    error('here -- check for NaNs case')
end

% IMPORTANT: we want the pullback to be an unwrapping of the surface that
%   respects APDV as viewed from the OUTSIDE of the mesh.
% This has implications in the direction that phi0 runs in real space, so
%   if we adjust phi0 as in phiOffsetsFromPrevMesh() or 
%   fitPhiOffsetsFromPrevMesh(), then we need to ensure the sign is
%   correct. To do that, we set the direction of v earlier based on flipy.

% Now use the axis orders specified
if isfield(Options, 'numLayers')
    if iscell(Options.numLayers)
        if length(Options.numLayers) == 1
            % this is a boring case: the cell just contains one pair of values,
            % so do the usual stack with a single channel or uniform numLayers
            % across all channels
            patchIm = texture_patch_to_image( TF, TV2D,...
                TF, TV3D(:, axisorder), IV, Options );
            % profile viewer
        else
            % There is a different Layer setting for each channel

            % Check that the #(layer settings) == #channels
            try
                assert(iscell(IV))
                assert(length(IV) == length(Options.numLayers))
            catch
                error('If numLayers is a cell, there must be one pair of integers for numLayers per channel of IV')
            end

            % If save_as_stack is true, we need the number of image planes to
            % be the same across channels. 
            if save_as_stack
                try
                    nplanes = zeros(length(Options.numLayers), 1) ;
                    for ch = 1:length(IV)
                        nplanes(ch) = 1 + sum(Options.numLayers{ch}) ;
                    end
                    assert(all(nplanes == nplanes(1)))
                catch
                    error('If saving as stack, #image planes in numLayers must be the same across channels')
                end
            end

            % Obtain each channel's stack, with different numLayers param
            for ch = 1:length(IV)
                disp(['selected layer pullback on channel ' num2str(ch)])
                Options2 = Options ;
                Options2.numLayers = Options.numLayers{ch} ;
                patchIms{ch} = texture_patch_to_image( TF, TV2D,...
                    TF, TV3D(:, axisorder), IV{ch}, Options2 );
            end

            % Prepare colors for each channel
            if isfield(Options, 'falseColors')
                falseColors = Options.falseColors ;
            else
                if length(IV) == 2
                    disp('aux_generate_orbifold: falseColors not supplied --> Using default 2-color channels of red and cyan')
                    falseColors{1} = [1, 0, 0] ;
                    falseColors{2} = [0, 1, 1] ;
                elseif length(IV) == 3
                    falseColors{1} = [1, 0, 0] ;
                    falseColors{2} = [0, 1, 0] ;
                    falseColors{3} = [0, 0, 1] ; 
                else
                    error('aux_generate_orbifold: For texture spaces with more than 3 channels, please supply colors for each channel in Options.falseColors')
                end
            end

            % Apply color to each channel
            if save_as_stack
                patchIm = zeros(size(patchIms{1})) ;
                for ch = 1:length(IV)
                    patchIm(:, :, 1, :) = patchIm(:, :, 1, :) + falseColors{ch}(:, 1) * patchIms{ch} ;
                    patchIm(:, :, 2, :) = patchIm(:, :, 2, :) + falseColors{ch}(:, 2) * patchIms{ch} ;
                    patchIm(:, :, 3, :) = patchIm(:, :, 3, :) + falseColors{ch}(:, 3) * patchIms{ch} ;
                end
            else
                tmp = max(patchIms{1}, [], 3) ;
                patchIm = zeros(size(tmp, 1), size(tmp, 2), 3) ;
                for ch = 1:length(IV)
                    dat = max(patchIms{ch}, [], 3) ;
                    patchIm(:, :, 1) = patchIm(:, :, 1) + falseColors{ch}(:, 1) * dat ;
                    patchIm(:, :, 2) = patchIm(:, :, 2) + falseColors{ch}(:, 2) * dat ;
                    patchIm(:, :, 3) = patchIm(:, :, 3) + falseColors{ch}(:, 3) * dat ;
                end
            end
        end

    else
        patchIm = texture_patch_to_image( TF, TV2D,...
            TF, TV3D(:, axisorder), IV, Options );
        % profile viewer
    end
else    
    disp('no numLayers defined, using default...')
    patchIm = texture_patch_to_image( TF, TV2D,...
        TF, TV3D(:, axisorder), IV, Options );
    % profile viewer
end

fprintf('Done\n');

% View results --------------------------------------------
% imshow( patchIm );
% set( gca, 'YDir', 'Normal' );

% Write figure to file
disp(['Writing ' imfn]) 
if length(size(patchIm)) < 3
    % Image is 2d, save using imwrite
    imwrite( patchIm, imfn, 'TIFF' ) ;
elseif length(size(patchIm))==3 && size(patchIm,3)==3
    % Image is RGB/FalseColor, save using imwrite
    imwrite( patchIm, imfn, 'TIFF' ) ;
elseif save_as_stack
    % image is 3d
    disp('Saving using saveastiff()')
    tiffoptions.overwrite = true ;
    dat = uint8(255 * patchIm) ;
    saveastiff( dat, imfn, tiffoptions) ;
elseif length(size(patchIm)) ==4 && size(patchIm, 3) == 3
    % Image is an RGB/FalseColor stack, take mip
    disp('Saving RGB MIP image')
    dat = max(patchIm, [], 4) ;
    dat = uint8(255 * dat) ;
    imwrite( dat, imfn, 'TIFF') ;
else
    disp('Taking MIP and saving as image')
    dat = max(patchIm, [], 3) ;
    dat = uint8(255 * dat) ;
    imwrite( dat, imfn, 'TIFF') ;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Save extended relaxed image
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% disp('Generating relaxed, extended image...')          
% % Format axes
% xlim([0 ar]); ylim([-0.5 1.5]);
% 
% % Extract image from figure axes
% patchIm_e = getframe(gca);
% patchIm_e = rgb2gray(patchIm_e.cdata);
% 
% % Write figure to file
% imwrite( patchIm_e, ...
%     sprintf( fullfile([imFolder_re, '/', fileNameBase, '.tif']), t ), ...
%     'TIFF' );

% Close open figures
close all