function measureSurfaceAreaVolume(QS, options)
%MEASURESURFACEAREAVOLUME(QS, OPTIONS)
%   Compute the mesh surface area and volume over time
%
% Prerequisites
% -------------
% QS.alignMeshesAPDV()
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   Current QuapSlap object
% options : struct with fields
%   overwrite : bool, default=false
%       overwrite results on disk
%   preivew : bool, default=QS.plotting.preview
%       display intermediate results
%   thres_fraca : float, default=1.0
%       threshold fractional change in area above which we say the area and 
%       volume are NaN 
%   
% Returns
% -------
% <none>
%
% NPMitchell 2020

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
timePoints = QS.xp.fileMeta.timePoints ;
preview = QS.plotting.preview ;
% [rot, trans] = QS.getRotTrans() ;
% resolution = QS.APDV.resolution ;
clineDVhoopBase = QS.fullFileBase.clineDVhoop ;
writheDir = QS.dir.writhe ;
meshDir = QS.dir.mesh ;
% get timestamps to indicate on writhe plot
% [~, ~, xyzlim_um] = getXYZLims(QS) ;
% xlims = xyzlim_um(1, :) ;
% ylims = xyzlim_um(2, :) ;
% zlims = xyzlim_um(3, :) ;
outdir = QS.dir.mesh ;
% Set t0 if possible, otherwise set to zero
t0 = QS.t0set() ;
if isempty(t0)
    t0 = 0 ;
end
try
    load(QS.fileName.fold, 'fold_onset') ;
catch
    fold_onset = [] ;
end
    
         

%% Unpack options
if nargin < 2
    options = struct() ;
end
% overwrite QS preview option if present in local options
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'thres_fraca')
    thres_fraca = options.thres_fraca ;
else
    thres_fracda = 1.0 ;
end

%%
datfn = fullfile(outdir, 'surfacearea_volume_stab.mat') ;
if ~exist(datfn, 'file') || overwrite
    % Prepare for iteration by preallocating memory for volume and area
    % lists
    vvs = zeros(length(timePoints), 1) ;
    aas = zeros(length(timePoints), 1) ;
     %% Compute surface area and volume for each mesh
    for ii = 1:length(timePoints)
        tp = timePoints(ii) ;

        % Load the mesh -- note we load the aligned mesh, before cutting
        % into cylinder
        meshfn = sprintf(QS.fullFileBase.alignedMesh, tp) ;
        disp(['analyzing ', meshfn])
        [tri, pts] = ply_read(meshfn, 'tri') ;

        % View Result --------------------------------------------------------
        if preview
            trisurf(tri, pts(:,1), pts(:,2), pts(:,3));
            axis equal
        end

        % Compute surface area and volume
        [vv, aa] = meshVolumeArea(pts, tri) ;
        if ii > 1
            fracda = abs(abs(aa) - abs(prev_a)) / abs(prev_a) ;
            disp(['fracda = ', num2str(fracda)])
            if fracda > thres_fracda
                vvs(ii) = NaN ;
                aas(ii) = NaN ;
                disp('Fractional change in aa is too great, skipping')
            else
                vvs(ii) = abs(vv) ;
                aas(ii) = aa ;
                prev_a = aa ;
            end
        else
            vvs(ii) = abs(vv) ;
            aas(ii) = aa ;
            prev_a = aa ;
            disp(['First step: aa = ', num2str(aa)])
        end
    end

    %% Take derivatives
    % load(fullfile(outdir, 'surfacearea_volume_stab.mat'))
    % Filter the data
    windowSize = 7; 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    vsm = smoothdata(vvs, 'rlowess', 5) ;
    vsmooth = filter(b, a, vsm) ;
    asm = smoothdata(aas, 'rlowess', 5) ;
    asmooth = filter(b, a, aas) ;
    % optional other smoothing
    % asmooth2 = smoothdata(ass, 'rlowess', 5) ;
    % vsmooth2 = smoothdata(vss, 'rlowess', 11);

    da = gradient(asmooth) ;
    dv = gradient(vsmooth) ;

    %% Save the data
    save(datfn, 'aas', 'vvs', 'da', 'dv')
else
    load(datfn, 'aas', 'vvs')
end

%% Save the image
outfn_pdf = fullfile(outdir, 'area_volume_over_time_stab.pdf') ;
outfn_png = fullfile(outdir, 'area_volume_over_time_stab.png') ;
if ~exist(outfn_png, 'file') || overwrite
    figh = figure();
    hold on;
    t0idx = QS.xp.tIdx(t0) ;
    ah = plot((timePoints - t0) * QS.timeInterval, aas / aas(t0idx)) ;
    vh = plot((timePoints - t0) * QS.timeInterval, vvs / vvs(t0idx)) ; 
    if ~isempty(fold_onset)
        if isnumeric(fold_onset(1)) 
            tidx = QS.xp.tIdx(fold_onset(1)) ;
            plot((fold_onset(1) - t0)* QS.timeInterval, aas(tidx) / aas(t0idx), 's') ;
        end
        if isnumeric(fold_onset(2)) 
            tidx = QS.xp.tIdx(fold_onset(2)) ;
            plot((fold_onset(2) - t0)* QS.timeInterval, aas(tidx) / aas(t0idx), 'o') ;
        end
        if isnumeric(fold_onset(3)) 
            tidx = QS.xp.tIdx(fold_onset(3)) ;
            plot((fold_onset(3) - t0) * QS.timeInterval, aas(tidx) / aas(t0idx), '^') ;
        end
    end
    legend({'area', 'volume'}, 'location', 'northwest')
    title('Surface area and volume')
    xlabel('Time [min]')
    ylabel('Normalized area or volume')
    saveas(figh, outfn_pdf)
    saveas(figh, outfn_png)
end
