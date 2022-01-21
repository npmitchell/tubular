function measureCellDensity(QS, nuclei_or_membrane, options)
% measureCellDensity(QS, nuclei_or_membrane, options)
%   Measure the cell density over time using relaxed smoothed-mesh 
%   pullbacks.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% nuclei_or_membrane : str specifier ('nuclei' or 'membrane')
%   whether to use membrane or nuclear method to measure
%   density
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   ws_thres_height : float
%       threshold height for watershed algorithm to avoid oversampling
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%   minDist : float, 
%       min distance as fraction of vmax between particles to
%       associate them as identical
%   swapAxes : bool
%       swap iLastik output axes (transpose) to match physical
%
%
% Returns
% -------
% <none>
%
% NPMitchell 2020

%% Unpack options
if isfield(options, 'ws_thres_height')
    ws_thres_height = options.ws_thres_height ;
else
    ws_thres_height = 1.1 ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'preview')
    preview = options.preview ;
else
    preview = false ;
end
if isfield(options, 'minDist')
    minDist = options.minDist ;
else
    minDist = 1e-3 ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
else
    timePoints = QS.xp.fileMeta.timePoints ;
end
if isfield(options, 'swapAxes')
    swapAxes = options.swapAxes ;
else
    swapAxes = true ;
end

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
imDir = fullfile(QS.dir.cellID, 'images') ;
if ~exist(imDir, 'dir')
    mkdir(imDir)
end
[~, ~, ~, xyzlim] = QS.getXYZLims() ;
% Add a bit of axis limit buffer for plotting
xyzlim = xyzlim + [-5, 5] ;
t0 = QS.t0set() ;

%% Extract cell density for each timepoint
if contains(nuclei_or_membrane, 'nucle')
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;        
        datfn = sprintf(QS.fullFileBase.cellID, tp) ;
        figfn = fullfile(imDir, [sprintf(QS.fileBase.cellID, tp) '.png']) ;
        
        if ~exist(datfn, 'file') || overwrite
            tic
            % Obtain current (s,phi) cutMesh
            QS.setTime(tp) ;
            QS.loadCurrentSPCutMeshSm() ;
            mesh = QS.currentMesh.spcutMeshSm ;
            mesh.v = QS.xyz2APDV(mesh.v) ;
            [ TF, TV2D, TV3D, ~ ] = tileAnnularCutMesh(mesh, [1,1]) ;

            % Load the relaxed smoothed pullback image
            fn = sprintf(QS.fullFileBase.im_r_sme, tp) ;
            im = imread(fn) ;
            % se = strel('disk', 15);
            % im2 = imopen(im, se)

            % Load probabilities
            dfn = sprintf(QS.fullFileBase.cellProbabilities, tp) ;
            dat = h5read(dfn, '/exported_data') ;

            % Binarize the data
            bw = squeeze(dat(1, :, :)) > 0.5 ;
            DD = bwdist(~bw) ;

            % Convert probabilities into points
            I2 = imcomplement(DD);
            I2(~bw) = Inf;
            %20 is the height threshold for suppressing shallow minima
            I3 = imhmin(I2, ws_thres_height); 
            disp(['performing watershed: t=' num2str(tp)])
            LL = watershed(I3);

            % find the connected components
            % use bwconncomp to extract the region of the segmented cell
            CC = bwconncomp(LL);
            centroids = regionprops(CC,'Centroid');

            % transform struct to array 
            XY = zeros( length(centroids),2);
            for i = 1:length(centroids)
                XY(i,1) = centroids(i).Centroid(1);
                XY(i,2) = centroids(i).Centroid(2);
            end
            
            % Check for axis permutation / column swapping from iLastik
            if swapAxes
                disp('iLastik data is transposed, swapping axes')
                XY = XY(:, [2, 1]) ;
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                Xsz = size(im, 1) ;
                Ysz = size(im, 2) ;
            end
            
            % Check the cellIDs
            if preview
                close all
                imshow(im); hold on;
                plot(XY(:, 1), XY(:, 2), '.'); 
                waitfor(gcf)
            end
            
            % transform into uv coords
            umax = max(mesh.u(:, 1)) ;
            vmax = max(mesh.u(:, 2)) ;
            doubleCovered = true ;
            xy = QS.XY2uv(im, XY, doubleCovered, umax, vmax) ;

            % Check the cellIDs
            if preview
                imshow(im); hold on;
                plot(XY(:, 1), XY(:, 2), 'o')
                plot((xy(:, 1) * (Xsz-1)) / (1*umax) + 1 , ...
                     (xy(:, 2) * (Ysz-1)) / (2*vmax) + 0.75 ...
                     + size(im,1)*0.25, '.'); 
                waitfor(gcf)
            end
            
            % Limit cells to single cover
            singleCover = (xy(:, 2) > 0) & (xy(:, 2) < vmax) ;
            XY = XY(singleCover, :) ;
            xy = xy(singleCover, :) ;
            % tile xy in as-natural-as-possible coordinates for
            % Delaunay triangulation, with U=(0,aspect) V=(0,1)
            % Note that we multiply by two since the image is double cover
            aspect = 2 * Xsz / Ysz ;
            xy4tri = zeros(size(xy)) ;
            xy4tri(:, 1) = aspect * xy(:, 1) / umax ;
            xy4tri(:, 2) = xy(:, 2) / vmax ;
            % tile xy in the natural coordinates U=(0,aspect) V=(0,1)
            txy = [xy4tri; xy4tri - [0,1]; xy4tri + [0,1]] ;
            [mIDx, mD] = knnsearch(txy, txy, 'K', 2) ;
            if any(mD(:, 2) < minDist)
                repeatID = mD(:, 2) < minDist ;
                mIDx(repeatID, :)
                error(['Handle case of repeated particle here by ', ...
                    'removing one of the near particles for each pair'])
                txy = txy(filtered, :) ;
                xy = txy(in_single, :) ;
                % XY = XY(filtered, :)
                XY = XY(filtered, :) ;
            else
                Ncells = size(xy, 1) ;
            end
            
            % Check 
            if preview
                imshow(im); hold on;
                plot(XY(:, 1), XY(:, 2), 'o')
                % map xy back to image space
                XYp = QS.uv2XY(im, xy, doubleCovered, umax, vmax) ;
                plot(XYp(:, 1), XYp(:, 2), '.'); 
                waitfor(gcf)
            end
            
            % Add extra points along boundaries to eliminate cut-spanning
            % faces at endcaps
            extentX = max(txy(:, 1)) - min(txy(:, 1)) ;
            addL = ones(300, 2) ;
            addR = ones(300, 2) ;
            addL(:, 1) = min(txy(:, 1)) - extentX * 0.01 ;
            addL(:, 2) = linspace(min(txy(:, 2)), max(txy(:, 2)), 300) ;
            addR(:, 1) = max(txy(:, 1)) + extentX * 0.01 ;
            addR(:, 2) = linspace(min(txy(:, 2)), max(txy(:, 2)), 300) ;
            n2keep = size(txy, 1) ;
            txy0 = [txy; addL; addR] ;
            
            % Triangulate points in pullback space 
            faces0 = delaunay(txy0) ; 
            
            % Remove padding points we added on Left and Right
            pt2remove = (n2keep+1):size(txy0, 1) ;
            [faces1, ~, ~] = remove_vertex_from_mesh(faces0, txy0, pt2remove) ;
            % assert(all(txy(:) == txy_check(:)))
            
            % Remove faces that are purely outside of the singleCover
            finc = ismember(faces1, 1:Ncells) ;
            faces1 = faces1(any(finc, 2), :) ;
            % Remove faces with duplicate points
            if any(faces1(:, 2) == faces1(:, 1)) || ...
                     any(faces1(:, 3) == faces1(:, 1)) || ...
                      any(faces1(:, 3) == faces1(:, 2))
                resp = input('faces have duplicate vertices. Proceed?', 'y') ;
                if contains(resp, 'y') || contains(resp, 'Y')
                    disp('Cleaning out repeats')
                    faces1 = faces1(faces1(:, 2) ~= faces1(:, 1), :) ;
                    faces1 = faces1(faces1(:, 3) ~= faces1(:, 1), :) ;
                    faces1 = faces1(faces1(:, 3) ~= faces1(:, 2), :) ;
                else
                    error('Exiting as demanded')
                end
            end
                        
            % Mod out the double cover from faces
            faces = faces1 ;
            faces(faces > Ncells) = faces(faces > Ncells) - Ncells ;
            faces(faces > Ncells) = faces(faces > Ncells) - Ncells ;
            assert(max(faces(:)) == Ncells)
            assert(all(faces(:) > 0))
            
            
            % Check
            if preview
                trimesh(faces0, txy0(:, 1), txy0(:, 2), 0*txy0(:, 2), ...
                    'edgecolor', 'b')
                hold on;
                trimesh(faces1, txy(:, 1), txy(:, 2), 0*txy(:, 2), ...
                    'edgecolor', 'k')
                trimesh(faces, txy(:, 1), txy(:, 2), 0*txy(:, 2), ...
                    'edgecolor', 'g')
                view(2);
                axis equal ;
                waitfor(gcf)
            end
            
            
            % Find positions in 3D of given points
            xyz = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, xy) ;
            
            % Remove faces with duplicate points
            if any(faces(:, 2) == faces(:, 1)) || ...
                     any(faces(:, 3) == faces(:, 1)) || ...
                      any(faces(:, 3) == faces(:, 2))
                
                dupl = faces(:, 2) == faces(:, 1) | ...
                       faces(:, 3) == faces(:, 1) | ...
                       faces(:, 3) == faces(:, 2) ;
                  
                % Show the issue of duplicated vertices  
                subplot(1, 2, 1)
                trimesh(faces(dupl, :), ...
                    xyz(:,1), xyz(:,2), xyz(:,3), xyz(:,2)) ;
                subplot(1, 2, 2)
                trimesh(faces, ...
                    xyz(:,1), xyz(:,2), xyz(:,3), xyz(:,2)) ;
                xlabel('x'); ylabel('y'); zlabel('z')
                set(gcf, 'visible', 'on')
                
                response = input('faces have duplicate vertices. Proceed?', 's') ;
                if ~contains(response, 'y') && ~contains(response, 'Y')
                    error('Exiting as demanded')
                end
            end
                        
            % Check
            if preview                
                trimesh(faces, xyz(:,1), xyz(:,2), xyz(:,3), xyz(:,2))
                xlabel('x'); ylabel('y'); zlabel('z')
                waitfor(gcf)
                trimesh(faces, xyz(:,1), xyz(:,2), xyz(:,3), ...
                    xyz(:,2), 'edgecolor', 'none', 'facecolor', 'interp')
                xlabel('x'); ylabel('y'); zlabel('z')
                waitfor(gcf)
            end
            
            % Compute voronoi areas on curves 3D mesh
            areas = 0.5 * doublearea(xyz, faces) ;
            % Find which areas (ie faces/triangles) are shared by each vertex
            % in xyz
            try
                tri = triangulation(faces, xyz) ;
            catch
                disp('Triangulation not possible: clipping to avoid NaNs')
                eps = 1e-12 ;
                maxvalX = umax - eps ;
                if max(xy(:, 1)) > maxvalX
                    xy(xy(:, 1) > maxvalX, 1) = umax - eps ;
                end
                maxvalY = vmax - eps ;
                if max(xy(:, 2)) > maxvalY
                    xy(xy(:, 2) > maxvalY, 2) = vmax - eps ;
                end
                if min(xy(:, 1)) < eps
                    xy(xy(:, 1) < eps, 1) = eps ;
                end
                if max(xy(:, 2)) < eps
                    xy(xy(:, 2) < eps, 2) = eps ;
                end
                try
                    xyz = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, xy) ;
                    tri = triangulation(faces, xyz) ;
                catch
                    close all
                    set(gcf, 'visible', 'on')
                    problem = find(any(isnan(xyz), 2)) ;
                    trimesh(TF, TV2D(:, 1), TV2D(:, 2), 0*TV2D(:, 1)) ; hold on;
                    plot3(xy(:, 1), xy(:, 2), 0*xy(:, 1), '.')
                    plot3(xy(problem, 1), xy(problem, 2), 0*xy(problem, 1), 'o')
                    title('Problem interpolating into 3D: problem pt is circled')
                    waitfor(gcf)
                    error('not sure what is failing. Debug here')
                end
            end
            
            attachments = vertexAttachments(tri) ;
            vorarea = zeros(size(xyz, 1), 1) ;
            for qq = 1:length(vorarea)
                faceneighbors = attachments{qq} ;
                vorarea(qq) = 1/3 * sum(areas(faceneighbors)) ;
            end
            
            % Smooth the areas just to check that this works
            % lambda = 0.1 ;
            % vorarea_sm = laplacian_smooth(xyz, faces, 'cotan', [],...
            %                         lambda, 'explicit', vorarea) ;
            % if any(isnan(vorarea_sm))
            %     error('some nans in vorarea smoothed')
            % end
                                
            % Compute binned data
            % binData2dGrid(uvz, uminmax, vminmax, nU, nV)
            % bin in u
            du = diff(mesh.u(1:nU, 1)) ;
            % check that all du increments are equal
            assert(all(abs(du - du(1)) / abs(du(1)) < 1e-10)) ;
            du = du(1) ;
            uedges = [mesh.u(1,1) - du*0.5; mesh.u(1:nU, 1) + du*0.5] ;
            uBinID = discretize(xy(:, 1), uedges) ;
            % bin in v
            dv = diff(mesh.u(1:nU:nV*nU, 2)) ;
            % check that all dv increments are equal
            assert(all(abs(dv - dv(1)) / abs(dv(1)) < 1e-10)) ;
            dv = dv(1) ;
            vedges = [mesh.u(1,2) - dv*0.5; mesh.u(1:nU:nU*nV, 2) + dv*0.5] ;
            vBinID = discretize(xy(:, 2), vedges) ;
            
            % Save triple cover
            xytile = [xy; xy-[0,1]; xy+[0,1]] ;
            cells = struct() ;
            cells.tripleCover = struct() ;
            cells.tripleCover.xy = xytile ;
            cells.tripleCover.XY = QS.uv2XY(im, xytile, ...
                                    doubleCovered, umax, vmax) ;
            cells.tripleCover.faces = faces0 ;
            cells.tripleCover.vorarea = [vorarea; vorarea; vorarea] ;
            % Single cover
            cells.XY = XY ;
            cells.vorarea = vorarea ;
            cells.faces = faces ;
            cells.uv = xy ;
            cells.xyz = xyz ;
            cells.umax = umax ;
            cells.vmax = vmax ;
            % Binned data
            cells.uBinID = uBinID ;
            cells.vBinID = vBinID ;

            % Check
            if preview                
                trisurf(faces, xyz(:,1), xyz(:,2), xyz(:,3), ...
                    vorarea, 'edgecolor', 'none')
                xlabel('x'); ylabel('y'); zlabel('z')
                waitfor(gcf)
                
                % Check ubinID
                trisurf(faces, xyz(:,1), xyz(:,2), xyz(:,3), ...
                    uBinID, 'edgecolor', 'none')
                title('uBinID')
                colorbar()
                axis equal 
                xlabel('x'); ylabel('y'); zlabel('z')
                waitfor(gcf)
                
                % Check vbinID
                trisurf(faces, xyz(:,1), xyz(:,2), xyz(:,3), ...
                    vBinID, 'edgecolor', 'none')
                title('vBinID')
                colorbar()
                axis equal 
                xlabel('x'); ylabel('y'); zlabel('z')
                waitfor(gcf)
            end
            
            % Save it
            error('here')
            save(datfn, 'cells')
            cells_in_RAM = true ;
            toc
        else
            cells_in_RAM = false ;
        end
        
        % Check that figures have been saved of the cell density
        if ~exist(figfn, 'file') || overwrite
            if ~cells_in_RAM
                disp(['Loading cells from disk:' datfn])
                load(datfn, 'cells')
                
                % load image too
                fn = sprintf(QS.fullFileBase.im_r_sme, tp) ;
                im = imread(fn) ;
            end
            
            % Save image of density
            xx = cells.xyz ;
            vor = cells.vorarea ;
            ff = cells.faces ;
            close all ;
            set(gcf, 'visible', 'off')
            npanel = 6 ;
            subplot(1, npanel, 4:npanel) ;
            trisurf(ff, xx(:, 1), xx(:, 2), xx(:, 3), vor, ...
                'edgecolor', 'none')
            axis equal
            xlim(xyzlim(1, :)) ;
            ylim(xyzlim(2, :)) ;
            zlim(xyzlim(3, :)) ;
            xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
            ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
            zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
            title(['cell area, $t=$' sprintf('%03d', tp - t0) ' ' QS.timeunits], ...
                 'Interpreter', 'Latex')
            % Colorbar
            cb = colorbar('location', 'southoutside') ;
            ylabel(cb, 'cell area [$\mu$m$^2$]', 'Interpreter', 'Latex') ;
            cpos = get(cb, 'position') ;
            cwidth = 0.25 ;
            cxpos = cpos(1) + 0.5*(cpos(3) - cwidth) - 0.015 ;
            set(cb,'position', [cxpos cpos(2) cwidth .02])
            caxis([0, 150])
            view(0, 0)
            
            % 2D cell density heatmap
            subplot(1, npanel, 1:2)
            imc = imcomplement(im) ;
            imshow(cat(3, imc, imc, imc)); hold on;
            XYtc = cells.tripleCover.XY ;
            disp('Plotting voronoi cells')
            [v, c] = voronoin(XYtc) ;
            for i=1:length(c)
                fill(v(c{i},1), v(c{i},2), cells.tripleCover.vorarea(i),...
                    'FaceAlpha', 0.35)
            end
            caxis([0, 150])
            % Limit plot to the single cover
            xlim([0, size(im, 2)])
            if doubleCovered
                ylim([0.25, 0.75] * size(im, 1))
            else
                ylim([0, size(im, 1)])
            end
            
            % Save the figure
            saveas(gcf, figfn)
            
        end
    end
elseif contains(nuclei_or_membrane, 'mem')
    error('have not considered this case yet')
else
    error(["Parameter nuclei_or_membrane must be string ", ...
        "like 'nuclei' or 'membrane' specifying which method to use"])
end



% Retired code
% 
% % Triangulate points in pullback space 
% faces = delaunay(XY) ;        
% 
% % Find positions in 3D of given points
% umax = max(mesh.u(:, 1)) ;
% vmax = max(mesh.u(:, 2)) ;
% doubleCovered = true ;
% xy = QS.XY2uv(im, XY, doubleCovered, umax, vmax) ;
% xyz = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, xy) ;
% 
% % Compute voronoi areas on curves 3D mesh
% areas = 0.5 * doublearea(xyz, faces) ;
% % Find which areas (ie faces/triangles) are shared by each vertex
% % in xyz
% tri = triangulation(faces, xyz) ;
% attachments = vertexAttachments(tri) ;
% vorarea = zeros(size(xyz, 1), 1) ;
% for qq = 1:length(vorarea)
%     faceneighbors = attachments{qq} ;
%     vorarea(qq) = 1/3 * sum(areas(faceneighbors)) ;
% end
% 
% % Store the details of cell identification results in a struct
% cells = struct() ;
% cells.doubleCover = struct() ;
% cells.doubleCover.vorarea = vorarea ;
% cells.doubleCover.faces = faces ;
% cells.doubleCover.xyz = xyz ;
% cells.doubleCover.XYpix = XY ;
% cells.doubleCover.uv = xy ;
% cells.umax = umax ;
% cells.vmax = vmax ;
% cells.singleCoverId = singleCover ; % boolean to limit cells to single cover
            

