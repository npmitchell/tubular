function aux_plot_constriction_DVhoops_transverse(folds, fold_onset_Idx, outDir, dvexten, ...
    save_ims, overwrite_lobeims, timePoints, spcutMeshBase, ...
    alignedMeshBase, normal_shift, rot, trans, resolution, colors, ...
    xyzlim, flipy, t0, timeInterval, timeUnits, spaceUnits)
%AUX_PLOT_AVGPTCLINE_LOBES auxiliary function for plotting the motion of
%the constrictions between lobes and the centerlines over time
% 
% Parameters
% ----------
% tp : N x 1 int array
%   xp.fileMeta.timePoints - min(fold_onset)
% timePoints : N x 1 int array
%   the timepoints in the experiment (xp.fileMeta.timePoints)
% 
% Returns
% -------
%
% NPMitchell 2020 

preview = false ;
tp = (timePoints - t0) * timeInterval ;
c1 = colors(1, :) ;
c2 = colors(2, :) ;
c3 = colors(3, :) ;
shift_mag = (normal_shift * resolution) + 1 ; 
xlims = xyzlim(1, :) + [-shift_mag, shift_mag] ;
ylims = xyzlim(2, :) + [-shift_mag, shift_mag] ;
zlims = xyzlim(3, :) + [-shift_mag, shift_mag] ;

fold_figfn = fullfile(outDir, 'constriction_transverse_%06d.pdf') ;
fns = dir(strrep(fold_figfn, '%06d', '*')) ;
if save_ims && (~(length(fns)==length(timePoints)) || overwrite_lobeims) || true
    disp('Creating constriction dynamics plots for locally transverse plane...')
    
    disp('Loading hoops for each fold...')
    tidx2do = fold_onset_Idx:10:length(timePoints) ;
    tidx2do = [ tidx2do, setdiff(1:length(timePoints), tidx2do) ];
    for kk = tidx2do
        % Translate to which timestamp
        t = timePoints(kk) ;
        disp(['t = ' num2str(t) '; tp = ' num2str(tp(kk))])
        tp4title = tp(kk) ;
        load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
        mesh = read_ply_mod(sprintf(alignedMeshBase, t)) ;
        % shift to match spcutMesh > note assume inward normal, but shift 
        % is outward, so subtract.
        mesh.v = mesh.v - normal_shift * resolution * mesh.vn ;
        nU = spcutMesh.nU ;
        nV = spcutMesh.nV ;

        % Load the centerline too
        avgpts = spcutMesh.avgpts ;

        % rename the fold indices (in U)
        f1 = folds(kk, 1) ;
        f2 = folds(kk, 2) ;
        f3 = folds(kk, 3) ;
        nfolds = size(folds, 2) ;
        if nfolds ~= 3
            error('Generalize to code to include > or < 3 folds/features')
        end

        % store distance from x axis of folds
        f1pt = avgpts(f1, :) ;
        f2pt = avgpts(f2, :) ;
        f3pt = avgpts(f3, :) ;

        % rotate vertices 
        vrs = ((rot * spcutMesh.v')' + trans) * resolution ;
        if flipy 
            vrs(:, 2) = - vrs(:, 2) ; 
        end
        
        % plot DVhoop
        f1 = double(f1) ;
        f2 = double(f2) ;
        f3 = double(f3) ;
        hoop1 = vrs(f1+nU*(0:(nV-1)), :) - avgpts(f1, :) ;
        hoop2 = vrs(f2+nU*(0:(nV-1)), :) - avgpts(f2, :) ;
        hoop3 = vrs(f3+nU*(0:(nV-1)), :) - avgpts(f3, :) ;
        
        % Flatten each hoop by measuring tangent of centerline at this
        % location
        f1segs = avgpts(f1-1:f1+1, :) ;
        tang1 = diff(f1segs) ;
        % average the tangent vecs on centerline in front and behind fold
        tang1 = mean(tang1 ./ vecnorm(tang1, 2, 2)) ;
        f2segs = avgpts(f2-1:f2+1, :) ;
        tang2 = diff(f2segs) ;
        % average the tangent vecs on centerline in front and behind fold
        tang2 = mean(tang2 ./ vecnorm(tang2, 2, 2)) ;
        f3segs = avgpts(f3-1:f3+1, :) ;
        tang3 = diff(f3segs) ;
        % average the tangent vecs on centerline in front and behind fold
        tang3 = mean(tang3 ./ vecnorm(tang3, 2, 2)) ;
        
        % The "dorsal pointing" vector to take to z axis
        dz1 = hoop1(1, :) ;
        dz2 = hoop2(1, :) ;
        dz3 = hoop3(1, :) ;
        
        
        % THis isn't always great fit, so we can fit to a plane explicitly
        ptCloud1 = pointCloud(hoop1) ;
        ptCloud2 = pointCloud(hoop2) ;
        ptCloud3 = pointCloud(hoop3) ;
        [model1,inlierIndices,outlierIndices] = pcfitplane(ptCloud1,...
            max(abs(hoop1(:))),tang1,45) ;
        [model2,inlierIndices,outlierIndices] = pcfitplane(ptCloud2,...
            max(abs(hoop1(:))),tang2,45) ;
        [model3,inlierIndices,outlierIndices] = pcfitplane(ptCloud3,...
            max(abs(hoop1(:))),tang3,45) ;
        assert(isempty(outlierIndices))
        
        rot1 = rotate3dToAlignAxis(model1.Normal, dz1) ;
        rot2 = rotate3dToAlignAxis(model2.Normal, dz2) ;
        rot3 = rotate3dToAlignAxis(model3.Normal, dz3) ;
        
        % Now actively rotate
        hoop1r = (rot1 * hoop1')' ;
        hoop2r = (rot2 * hoop2')' ;
        hoop3r = (rot3 * hoop3')' ;
        
        
        % preview each fit to plane
        if preview
            clf; 
            subplot(1, 3, 1)
            plot3(hoop1(:, 1), hoop1(:, 2), hoop1(:, 3), '.')
            hold on
            plot3(hoop1r(:, 1), hoop1r(:, 2), hoop1r(:, 3), '.')
            tgv1 = [0,0,0; 50*tang1] ;
            tgv1r = (rot1*tgv1')' ;
            plot3(tgv1(:, 1), tgv1(:, 2), tgv1(:, 3), '-');
            hold on
            plot3(tgv1r(:, 1), tgv1r(:, 2), tgv1r(:, 3), '-')
            axis equal
            
            subplot(1, 3, 2)
            plot3(hoop2(:, 1), hoop2(:, 2), hoop2(:, 3), '.')
            hold on
            plot3(hoop2r(:, 1), hoop2r(:, 2), hoop2r(:, 3), '.')
            tgv2 = [0,0,0; 50*tang2] ;
            tgv2r = (rot2*tgv2')' ;
            plot3(tgv2(:, 1), tgv2(:, 2), tgv2(:, 3), '-');
            hold on
            plot3(tgv2r(:, 1), tgv2r(:, 2), tgv2r(:, 3), '-')
            axis equal
            
            subplot(1, 3, 3)
            plot3(hoop3(:, 1), hoop3(:, 2), hoop3(:, 3), '.')
            hold on
            plot3(hoop3r(:, 1), hoop3r(:, 2), hoop3r(:, 3), '.')
            tgv3 = [0,0,0; 50*tang3] ;
            tgv3r = (rot3*tgv3')' ;
            plot3(tgv3(:, 1), tgv3(:, 2), tgv3(:, 3), '-');
            hold on
            plot3(tgv3r(:, 1), tgv3r(:, 2), tgv3r(:, 3), '-')
            axis equal
            
        end
        
        % Plot each hoop
        close all
        fig = figure('visible', 'off'); 

        sz = 20 ;
        msz = 50 ;
        for qq = 1:nfolds
            subplot(1, 3, qq)
            hold on;
            if qq == 1
                plot(- hoop1r(:, 2), hoop1r(:, 3), '-', 'color', c1);
            elseif qq == 2
                plot(-hoop2r(:, 2), hoop2r(:, 3), '-', 'color', c2);
            else
                plot(-hoop3r(:, 2), hoop3r(:, 3), '-', 'color', c3);
            end
        axis equal
        xlim(ylims)
        ylim(zlims)
        xlabel(['AP position [' spaceUnits ']'], 'interpreter', 'latex')
        ylabel(['lateral position [' spaceUnits ']'], 'interpreter', 'latex')
        zlabel(['DV position [' spaceUnits ']'], 'interpreter', 'latex')
        end
        sgtitle(sprintf(['Constriction dynamics, t=%03d ' timeUnits], tp4title), ...
            'interpreter', 'latex')
        ofn = sprintf(fold_figfn, t) ;
        saveas(fig, ofn) ;
        close all
    end
else
    disp(['Skipping constriction hoop dynamics plots since they exist (' fold_figfn ')...'])
end


