function com = com_region(probability_grid, thres, varargin)
    % com = COM_REGION(probability_grid, thres, gridxyz)
    % Given a probability cloud, find the center of mass of the largest 
    % connected region of the probability cloud above
    % some threshold, with "mass" proportional to probability. 
    % The output is given in units where the left corner of the volume is
    % (0.5, 0.5, 0.5), so its the typical 1-indexed MATLAB output. Note
    % though that the axis order is PRESERVED unlike in other MATLAB
    % functions. That is, probability grid is assumed to be meshgrid-like,
    % not ngrid like, so XYZ --> XYZ. 
    %
    % INPUTS
    % ------
    % probability_grid
    % thres : float
    %   Threshold probability used to segment the probability cloud into
    %   connected regions
    % mesh_vertices : N x 3 float array
    %   Positions of the mesh vertices, as N x 3 array
    % varargin: struct
    %   Options struct, with fields xyzgrid and check
    %   options.xyzgrid: 3d float or int array 
    %       xyzgrid positional values matching the probability_grid
    %   options.check: boolean
    %       whether to display
    %   options.color : colorspec
    %       isosurface color if options.check is true
    %
    % OUTPUTS
    % -------
    % com : 3x1 float 
    %   the position of the cernter of mass of the chunk of probability 
    %   cloud above the supplied threshold thres
    %
    % EXAMPLE USAGE
    % -------------
    % probability_grid = zeros(20,20,20) ;
    % probability_grid(1:12,1:3,1:6) = 1 ; 
    % probability_grid(end-1:end,end-1:end,end-1:end) = 1 ; 
    % options.check = true ;
    % options.color = 'blue'; 
    % com = com_region(probability_grid, 0.5, options) 
    % >> com = [6.5, 2.0, 3.5] 
    % 
    % NPMitchell 2019-2020
    
    
    % Parse arguments
    if nargin > 2
        if isfield(varargin{1}, 'check')
            check = varargin{1}.check ;
        else
            check = false ;
        end
        
        if isfield(varargin{1}, 'waitForPress')
            waitForPress = varargin{1}.waitForPress ;
        else
            waitForPress = false ;
        end
        
        if isfield(varargin{1}, 'clearFigures')
            clearFigures = varargin{1}.clearFigures ;
        else
            clearFigures = true ;
        end
        
        if isfield(varargin{1}, 'check_slices')
            check_slices = varargin{1}.check_slices ;
        else
            check_slices = false ;
        end

        if isfield(varargin{1}, 'color')
            color = varargin{1}.color ;
        else
            color = 'red' ;
        end
        
        if isfield(varargin{1}, 'xyzgrid')
            xyzgrid = varargin{1}.xyzgrid;
        else
            xyzgrid = [] ;
        end
    else
        check = false ;
        check_slices = false;
        color = 'red'; 
        xyzgrid = [] ;
    end
    
    % Find big region, decreasing thres if needed
    bwccOk = false ;
    thresIter = thres ;
    while ~bwccOk
        try
            bwcc = bwconncomp(probability_grid > thresIter) ; 
            npix = cellfun(@numel,bwcc.PixelIdxList);
            [~, indexOfMax] = max(npix); 
            % isolate the largest connected component
            biggest = zeros(size(probability_grid));
            biggest(bwcc.PixelIdxList{indexOfMax}) = 1;
            bwccOk = true ;
        catch
            thresIter = thresIter * 0.9 ;
            disp(['Could not segment big region, decreasing threshold to ' num2str(thresIter)])
            if thresIter < 0.01
                error('Threshold of 0.01 still did not yield a connected component')
            end
        end
    end
    % Now multiply with probability
    mass = probability_grid .* biggest ;
    % Get center of mass. There are two ways
    if ~isempty(xyzgrid)
        % Method 1, using mean
        mass = probability_grid ;
        mean_mass = mean(mass(:)) ;
        comX = mean(mass(:) .* xyzgrid(1, :, :, :)) / mean_mass ;
        comY = mean(mass(:) .* xyzgrid(2, :, :, :)) / mean_mass ;
        comZ = mean(mass(:) .* xyzgrid(3, :, :, :)) / mean_mass ;
        com = [ comX comY comZ ] ; 
    else 
        props = regionprops(true(size(mass)), mass, 'WeightedCentroid');
        com_tmp = props.WeightedCentroid ;
        % Note that we flip XY here since WeightedCentroid treats the 1st
        % dimension as the 2nd (and 2nd as 1st)
        com = [ com_tmp(2) com_tmp(1) com_tmp(3) ] ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check the com
    if check
        if clearFigures
            close all
        end
        
        % Show relative to the mesh
        pprob = permute(probability_grid, [2 1 3]) ;
        iso = isosurface(pprob, 0.5) ;
        patch(iso,'facecolor',color,'facealpha',0.5,'edgecolor','none');
        view(3)
        camlight
        hold on
        scatter3(com(1), com(2), com(3))
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        xlim([0, size(probability_grid, 1)])
        ylim([0, size(probability_grid, 2)])
        zlim([0, size(probability_grid, 3)])
        hold off
        if waitForPress
            title('COM and pointcloud. Click any button to continue')

            % Show the figure
            disp('Showing com and pointcloud. Click any button to continue')
            drawnow;
            % disp(clock) ;
            try
                waitforbuttonpress
                % Close figure or leave it open
                disp('mouse or key pressed')
                % waitfor(gcf)
            catch
                disp('figure closed')
            end
        else
            % Show the figure
            title('COM and pointcloud')
            disp('Previewing COM and pointcloud. Continuing...')
            drawnow;
        end
        % disp(clock)        
    end
    
    % Show mass slices
    if check_slices
        % Plot each section of the intensity data 
        for jj=1:size(mass,1)
            if any(any(mass(jj, :, :)))
                imshow(squeeze(mass(jj,:,:))')
                title(['mass slice ' num2str(jj)])
                pause(0.00001)
            end
        end
        for jj=1:size(probability_grid,1)
            imshow(squeeze(probability_grid(jj,:,:)'))
            title(['prediction slice' num2str(jj)])
            pause(0.00001)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
return