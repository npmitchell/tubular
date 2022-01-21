classdef QuapSlap < handle
    % Quasi-Axisymmetric Pipeline for Surface Lagrangian Pullbacks class
    %
    % Coordinate Systems
    % ------------------
    % sphi : (proper length x rectified azimuthal coordinate)
    %       quasi-axisymmetric system in which first coordinate is 
    %       proper length along surface and second is a rectified azimuthal
    %       coordinate. Rectification means that the surface is rotated as
    %       v->phi(s), where each s coordinate is offset by a "rotation"
    %       about the surface along a direction that is perpendicular to s
    %       in pullback space. This rectification may be based on surface
    %       positions in R^3 (geometric) or based on intensity motion in
    %       pullback space (material/Lagrangian) inferred through 
    %       phasecorrelation of tissue strips around discretized s values.
    % uv :  (conformal map onto unit square)
    %       Conformally mapping the cylinderCutMesh onto the unit square 
    %       in the plane results in the instantaneous uv coordinate system. 
    %       Corners of the unit square are taken directly from the cutMesh,
    %       so are liable to include some overall twist.
    % uvprime : (conformal map
    %       [same as uvprime_sm, since uvprime is currently computed via
    %       sphi_sm coordinates]
    % r/spr/sphir :  (proper length x rectified azimuthal coordinate)
    %       same as sphi but with aspect ratio relaxed to minimize isoareal
    %       energy cost --> makes triangles more similar in area by scaling
    %       the s axis (longitudinal axis) by a scalar factor.
    %
    % iLastik's coordinate systems : I recommend using cxyz as the axis 
    %       when outputting from iLastik. For our setup, this will mean
    %       that the meshes are mirrored with respect to the lab frame, so
    %       QS.flipy = true.
    % 
    % PIV measurements fall into two classes: 
    %   - 'piv': principal surface-Lagrangian-frame PIV (sp_sme or up_sme)
    %       --> note that the designation of coordinate system is not
    %       explicitly specified in the filenames:
    %       QS.dir.mesh/gridCoords_nU0100_nV0100/piv/piv3d, etc
    %   - 'piv_uvp_sme': PIV in coordSys for quasiconformal measurements 
    %       --> note that these are less Lagrangian than
    %       sp_sme or up_sme, so they are treated as independent from the 
    %       principal pipeline in which we measure velocities in a
    %       surface-Lagrangian frame (ie sp_sme or up_sme).
    %
    % Properties
    % ----------
    % xyzlim        : 3x2 float, mesh limits in full resolution pixels, in data space
	% xyzlim_um     : 3x2 float, mesh limits in lab APDV frame in microns
    % resolution    : float, resolution of pixels in um
    % rot           : 3x3 float, APDV rotation matrix
    % trans         : 3x1 float, APDV translation 
    % a_fixed       : float, aspect ratio for fixed-width pullbacks
    % phiMethod     : str, '3dcurves' or 'texture'
    % flipy         : bool, APDV coord system is mirrored XZ wrt raw data
    % nV            : int, sampling number along circumferential axis
    % nU            : int, sampling number along longitudinal axis
    % uvexten       : str, naming extension with nU and nV like '_nU0100_nV0100'
    % t0            : int, reference timePoint in the experiment
    % features : struct with fields
    %   folds : #timepoints x #folds int, 
    %       indices of nU sampling of folds
    %   fold_onset : #folds x 1 float
    %       timestamps (not indices) of fold onset
    %   ssmax : #timepoints x 1 float
    %       maximum length of the centerline at each timepoint
    %   ssfold : #timepoints x #folds float
    %       positional pathlength along centerline of folds
    %   rssmax : #timepoints x 1 float
    %       maximum proper length of the surface over time
    %   rssfold : #timepoints x #folds float
    %       positional proper length along surface of folds
    % velocityAverage : struct with fields
    %   vsmM : (#timePoints-1) x (nX*nY) x 3 float array
    %       3d velocities at PIV evaluation coordinates in um/dt rs
    %   vfsmM : (#timePoints-1) x (2*nU*(nV-1)) x 3 float array
    %       3d velocities at face barycenters in um/dt rs
    %   vnsmM : (#timePoints-1) x (nX*nY) float array
    %       normal velocity at PIV evaluation coordinates in um/dt rs
    %   vvsmM : (#timePoints-1) x (nU*nV) x 3 float array
    %       3d velocities at (1x resolution) mesh vertices in um/min rs
    %   v2dsmM : (#timePoints-1) x (nX*nY) x 2 float array
    %       2d velocities at PIV evaluation coordinates in pixels/ min
    %   v2dsmMum : (#timePoints-1) x (nX*nY) x 2 float array
    %       2d velocities at PIV evaluation coordinates in scaled pix/min, but 
    %       proportional to um/min (scaled by dilation of map)
    %
    % APDV coordinate system and Centerline specification
    % ---------------------------------------------------
    % QuapSlap allows the user to designate the AP axis, a DV axis, and a
    % centerline for the segmented object. 
    % To designate APDV coordinate system, use iLastik training on the
    % timepoint t0 (which is an attribute of QS). In particular, train on
    % anterior (A), posterior (P), background (B), and dorsal (D) location 
    % in different iLastik channels. Small blobs of high probability near
    % the anterior end for A, somewhere along the posterior end for P, and 
    % anywhere along dorsal for D are best. Then the centers of mass of
    % thresholded contiguous high probability are computed for each to
    % define the coordinate system. By default, the ilastik results are
    % read as ch1=A, ch2=P, ch3=B, ch4=D, but anteriorChannel, 
    % posteriorChannel, and dorsalChannel specify the iLastik
    % training channel that is used for each specification.
    % Name the h5 file output from iLastik as 
    % ..._Probabilities_APDVcoords.h5
    % For example, dorsal for the gut was chosen at the fused site where 
    % additional 48YGAL4-expressing muscle-like cells form a seam.
    % Posterior is at the rear of the yolk, where the endoderm closes, for 
    % apical surface training. Anterior is at the junction of the midgut 
    % with the foregut.
    % Separately, define the AP points for centerline extraction. For most gut
    % data, the posterior point is different in centerline than it is for AP
    % centerline specification, so we use a different ILP to train for A and P
    % for all timepoints of dynamic data. 
    %
    % 
    properties
        xp                      % ImSAnE experiment class instance
        dynamic                 % true if multiple timepoints, false if fixed
        timeInterval = 1        % increment in time between timepoints with 
                                % indices differing by 1. For example, if
                                % timePoints are [0,1,2,4] and these are 
                                % [1,1,2] minutes apart, then timeInterval 
                                % is 1. 
        timeUnits = 'min'       % units of the timeInterval (ex 'min')
        spaceUnits = '$\mu$m'   % units of the embedding space (ex '$\mu$m')
        dir                     % str, directory where QuapSlap data lives
        dirBase                 % 
        fileName                % fileName
        fileBase
        fullFileBase
        ssfactor                % subsampling factor for probabilities 
        APDV = struct('resolution', [], ...
            'rot', [], ...
            'trans', [])
        flipy                   % whether data is mirror image of lab frame coordinates
        nV                      % sampling number along circumferential axis
        nU                      % sampling number along longitudinal axis
        uvexten                 % naming extension with nU and nV like '_nU0100_nV0100'
        t0                      % reference time in the experiment
        normalShift
        features = struct('folds', [], ...  % #timepoints x #folds int, indices of nU sampling of folds
            'fold_onset', [], ...       % #folds x 1 float, timestamps (not indices) of fold onset
            'ssmax', [], ...            % #timepoints x 1 float, maximum length of the centerline at each timepoint
            'ssfold', [], ...           % #timepoints x #folds float, positional pathlength along centerline of folds
            'rssmax', [], ...           % #timepoints x 1 float, maximum proper length of the surface over time
            'rssfold', []) ;            % #timepoints x #folds float, positional proper length along surface of folds
        a_fixed                         % aspect ratio for fixed geometry pullback meshes
        phiMethod = '3dcurves'          % method for determining Phi map in pullback mesh creation, with 
                                        % the full map from embedding to pullback being [M'=(Phi)o()o()]. 
                                        % This string specifier must be '3dcurves' (geometric phi stabilization) 
                                        % or 'texture' (optical flow phi stabilization)
        endcapOptions
        plotting = struct('preview', false, ... % display intermediate results
            'save_ims', true, ...       % save images
            'xyzlim_um_buff', [], ...   % xyzlimits in um in RS coord sys with buffer
            'xyzlim_raw', [], ...       % xyzlimits in pixels
            'xyzlim_pix', [], ...       % xyzlimits in pixels RS
            'xyzlim_um', [], ...        % xyzlimits in um in RS coord sys
            'colors', [])               % color cycle for QS
        apdvCOM = struct('acom', [], ...
            'pcom', [], ... 
            'acom_sm', [], ...
            'pcom_sm', [], ... 
            'dcom', [], ... 
            'acom_rs', [], ... 
            'pcom_rs', [], ... 
            'dcom_rs', [])
        apdvCOMOptions
        currentTime
        currentMesh = struct('rawMesh', [], ...
            'alignedMesh', [], ...     % APDV rotated and scaled mesh         
            'cylinderMesh', [], ...
            'cylinderMeshClean', [], ...
            'cutMesh', [], ...
            'cutPath', [], ...
            'spcutMesh', [], ...
            'spcutMeshSm', [], ...      
            'spcutMeshSmRS', [], ...    % rectilinear cutMesh in (s,phi) with rotated scaled embedding
            'spcutMeshSmRSC', [], ...   % rectilinear cutMesh as closed cylinder (topological annulus), in (s,phi) with rotated scaled embedding
            'ricciMesh', [], ...        % ricci flow result pullback mesh, topological annulus
            'uvpcutMesh', [])           % rectilinear cutMesh in (u,v) from Dirichlet map result to rectangle 
        currentCline = struct('mss', [], ...
            'mcline', [], ...
            'avgpts', []) ;
        data = struct('adjustlow', 0, ...
            'adjusthigh', 0, ...
            'axisOrder', [1 2 3], ...
            'ilastikOutputAxisOrder', 'cxyz') % options for scaling and transposing image intensity data
        currentData = struct('IV', [], ...
            'adjustlow', 0, ...
            'adjusthigh', 0 )           % image intensity data in 3d and scaling
        currentVelocity = struct('piv3d', struct(), ...
            'piv3d2x', struct()) ;     
        currentVelocityMultiChannel = struct('piv3d', cell(1), ...
            'piv3d2x', cell(1)) ;     
        currentSegmentation = struct(...
            'coordSys', 'spsme', ...    % pullback coordinate system in which segmentation is performed
            'seg2d', [], ...            % 2d pullback segmentation
            'seg3d', [], ...            % segmentation pushed forward into 3d
            'seg2dCorrected', [], ...
            'seg3dCorrected', []) ;     
        piv = struct( ...
            'imCoords', 'sp_sme', ...   % image coord system for measuring PIV / optical flow) ;
            'Lx', [], ...               % width of image, in pixels (x coordinate)
            'Ly', [], ...               % height of image, in pixels (y coordinate)
            'raw', struct(), ...        % raw PIV results from disk/PIVLab
            'smoothed', struct(), ...   % smoothed PIV results after gaussian blur
            'smoothing_sigma', 1 ) ;    % sigma of gaussian smoothing on PIV, in units of PIV sampling grid pixels
        pivMultiChannel = struct( ...   % Multi-channel PIV: one measurement for each image channel (for relative motion)
            'imCoords', 'sp_sme', ...   % image coord system for measuring PIV / optical flow) ;
            'Lx', [], ...               % width of image, in pixels (x coordinate)
            'Ly', [], ...               % height of image, in pixels (y coordinate)
            'raw', cell(1), ...         % raw PIV results from disk/PIVLab
            'smoothed', cell(1), ...    % smoothed PIV results after gaussian blur
            'smoothing_sigma', 1 ) ;    % sigma of gaussian smoothing on PIV, in units of PIV sampling grid pixels
        velocityAverage = struct(...
            'v3d', [], ...              % 3D velocities in embedding space [pix/dt]
            'v2d', [], ...              % 2D tangential velocities in pullback
            'v2dum', [], ...            % 2D tangential velocity scaled by speed in true embedding space
            'vn', [], ...               % normal velocity in spaceUnits per timeInterval timeUnits
            'vf', [], ...               % velocity vielf on face barycenters after Lagrangian avg
            'vv', []) ;                 % velocity field on vertices after Lagrangian avg
        velocityAverage2x = struct(...
            'v3d', [], ...
            'v2d', [], ...
            'v2dum', [], ...
            'vn', [], ...
            'vf', [], ...       
            'vv', []) ;                 % velocity field after Lagrangian avg after triangle subdivision
        velocitySimpleAverage = struct(... % assuming small in-plane motions, averaged veloticites in time at stationary (u,v) coordinates
            'v3d', [], ...              
            'v2d', [], ...
            'v2dum', [], ...
            'vn', [], ...
            'vf', [], ...       
            'vv', []) ;                 % velocity field after in-place (uv) avg assuming
        velocitySimpleAverage2x = struct(... % assuming small in-plane motions, averaged veloticites in time at stationary (u,v) coordinates
            'v3d', [], ...
            'v2d', [], ...
            'v2dum', [], ...
            'vn', [], ...
            'vf', [], ...       
            'vv', []) ;                 % velocity field after in-place (uv) avg
        cleanCntrlines          % centerlines in embedding space after temporal averaging
        % pivPullback = 'sp_sme'; % coordinate system used for velocimetry
        %  ---> instead use QS.piv.imCoords
        smoothing = struct(...
            'lambda', 0.002, ...            % diffusion const for field smoothing on mesh
            'lambda_mesh', 0.001, ...       % diffusion const for vertex smoothing of mesh itself
            'lambda_err', 0.005, ...        % diffusion const for fields inferred from already-smoothed fields on mesh
            'nmodes', 7, ...                % number of low freq modes to keep per DV hoop
            'zwidth', 1) ;                  % half-width of tripulse filter applied along zeta/z/s/u direction in pullback space, in units of du/dz/ds/dzeta
        pathlines = struct('t0', [], ...    % timestamp (not an index) at which pathlines form regular grid in space
            'refMesh', [], ...              % reference mesh for pathline advection
            'piv', [], ...                  % Lagrangian pathlines from piv coords
            'vertices', [], ...             % Lagrangian pathlines from mesh vertices
            'vertices3d', [], ...           % Lagrangian pathlines from mesh vertices
            'faces', [], ...                % Lagrangian pathlines from mesh face barycenters
            'featureIDs', struct(...        % struct with features in pathline coords
                'vertices', [], ...         % longitudinal position of features from pathlines threaded through pullback mesh vertices at t=t0Pathline
                'piv', [], ...              % longitudinal position of features from pathlines threaded through PIV evaluation coordinates at t=t0Pathline
                'faces', []), ...           % longitudinal position of features from pathlines threaded through pullback mesh face barycenters at t=t0Pathline
            'beltrami', ...                 % beltrami coefficient evaluated along pathlines
                struct('mu_material', [], ...
                'mu_material_filtered', [], ...
                'mu_material_vertices', [], ...
                'fitlerOptions', []));     
        pathlines_uvprime = struct('t0', [], ...    % timestamp (not an index) at which pathlines form regular grid in space
            'piv', [], ...                  % Lagrangian pathlines from piv coords
            'vertices', [], ...             % Lagrangian pathlines from mesh vertices
            'faces', [], ...                % Lagrangian pathlines from mesh face barycenters
            'featureIDs', struct(...        % struct with features in pathline coords
                'vertices', [], ...         % longitudinal position of features from pathlines threaded through pullback mesh vertices at t=t0Pathline
                'piv', [], ...              % longitudinal position of features from pathlines threaded through PIV evaluation coordinates at t=t0Pathline
                'faces', []));              % longitudinal position of features from pathlines threaded through pullback mesh face barycenters at t=t0Pathline
        currentStrain = struct(...
            'pathline', ...                 % strain from pathlines
                struct('t0Pathlines', [], ...   % t=0 timepoint for pathlines in question
                'strain', [], ...               % strain evaluated along pathlines
                'beltrami', [])) ;              % beltrami coefficient for pathlines
    end
    
    % Some methods are hidden from public view. These are used internally
    % to the class.
    methods (Hidden)
        function QS = QuapSlap(xp, opts)
            QS.initializeQuapSlap(xp, opts)
        end
        initializeQuapSlap(QS, xp, opts)
        plotSPCutMeshSmSeriesUtility(QS, coordsys, options)
        plotMetricKinematicsTimePoint(QS, tp, options)
        [XX, YY] = pullbackPathlines(QS, x0, y0, t0, options) 
        plotAverageVelocitiesTimePoint(QS, tp, options)
        plotPathlineVelocitiesTimePoint(QS, tp, options)
        plotStrainRateTimePoint(QS, tp, options) 
        plotPathlineStrainRateTimePoint(QS, tp, options)
        plotPathlineStrainTimePoint(QS, tp, options)
    end
    
    % Public methods, accessible from outside the class and reliant on 
    % properties of the class instance
    methods
        function setTime(QS, tt)
            % Set the current time of the dataset and clear current data
            % which was associated with the previously considered time
            %
            % Parameters
            % ----------
            % tt : int or float
            %   timePoint to set to be current, from available times in
            %   QS.xp.fileMeta.timePoints
            %
            if tt ~= QS.currentTime
                QS.clearTime() ;
            end
            QS.currentTime = tt ;
            if ~isa(QS.xp, 'struct')
                % disp('WARNING: QS.xp is interpreted as a struct rather than imsane Experiment() class instance')
                QS.xp.setTime(tt) ;
            end
        end
        
        function [tpTrue, timestr ] = trueTime(QS, tp, div60)
            if nargin < 3
                div60 = false ;
            end
            tpTrue = (tp - QS.t0set()) * QS.timeInterval ;
            if div60
                tpTrue = tpTrue / 60 ;
            end
            
            if nargout > 1
                if div60
                    if contains(QS.timeUnits, 'min')
                        timestr = [num2str(tpTrue) ' hr'] ;
                    elseif contains(QS.timeUnits, 'sec')
                        timestr = [num2str(tpTrue) ' min'] ;
                    end
                else
                    timestr = [sprintf('%03d', tpTrue) ' ' QS.timeUnits] ;
                end
            end
        end
        
        function clearTime(QS)
            % clear current timepoint's data for QS instance
            QS.currentTime = [] ;
            QS.currentMesh.rawMesh = [] ;
            QS.currentMesh.rawMeshRS = [] ;
            QS.currentMesh.cylinderMesh = [] ;
            QS.currentMesh.cylinderMeshClean = [] ;
            QS.currentMesh.cutMesh = [] ;
            QS.currentMesh.cutPath = [] ;
            QS.currentMesh.spcutMesh = [] ;
            QS.currentMesh.cutMesh = [] ;
            QS.currentMesh.spcutMesh = [] ;
            QS.currentMesh.spcutMeshSm = [] ;
            QS.currentMesh.spcutMeshSmRS = [] ;
            QS.currentMesh.spcutMeshSmRSC = [] ;
            QS.currentMesh.uvpcutMesh = [] ;
            QS.currentData.IV = [] ;
            QS.currentData.adjustlow = 0 ;
            QS.currentData.adjusthigh = 0 ;
            QS.currentVelocity.piv3d = struct() ;
            QS.currentVelocity.piv3d2x = struct() ;
            QS.currentVelocityMultiChannel.piv3d = ...
                cell(length(QS.xp.expMeta.channelsUsed), 1) ;
            QS.currentVelocityMultiChannel.piv3d2x = ...
                cell(length(QS.xp.expMeta.channelsUsed), 1);
            QS.currentSegmentation.seg2d = [] ;
            QS.currentSegmentation.seg3d = [] ;
            QS.currentSegmentation.seg2dCorrected = [] ;
            QS.currentSegmentation.seg3dCorrected = [] ;
        end
        
        function t0 = t0set(QS, t0)
            % t0set(QS, t0) Set time offset to 1st fold onset or manually 
            if nargin < 2
                if exist(QS.fileName.fold, 'file')
                    % Note that fold_onset is in units of timepoints, not 
                    % indices into timepoints
                    load(QS.fileName.fold, 'fold_onset') ;
                    QS.t0 = min(fold_onset) ;
                else
                    if numel(QS.xp.fileMeta.timePoints) == 1
                        QS.t0 = QS.xp.fileMeta.timePoints ;
                    else
                        error('No folding times saved to disk')
                    end
                end
            else
                QS.t0 = t0 ;
            end
            t0 = QS.t0 ;
        end
        
        function makeMIPs(QS, dim, pages, timePoints, adjustIV)
            if nargin < 5
                adjustIV = false ;
            end
            if nargin < 4 
                timePoints = QS.xp.fileMeta.timePoints;
            elseif isempty(timePoints)
                timePoints = QS.xp.fileMeta.timePoints;
            end
            if ~iscell(pages)
                pages = {pages} ;
            end
            % create mip directories if needed
            for qq = 1:length(pages)
                outdir = sprintf(QS.dir.mip, dim, ...
                            min(pages{qq}), max(pages{qq})) ;
                if ~exist(outdir, 'dir')
                    mkdir(outdir)
                end
            end
            
            % make the mips
            for tp = timePoints
                for qq = 1:length(pages)
                    im = QS.mip(tp, dim, pages{qq}, adjustIV) ;
                    imfn = sprintf(QS.fullFileBase.mip, dim, ...
                        min(pages{qq}), max(pages{qq}), tp) ;
                    imwrite(im, imfn,'tiff','Compression','none')
                end
            end
            disp('done making mips')
        end
        
        function im = mip(QS, tp, dim, pages, adjustIV)
            if nargin < 5
                adjustIV = false ;
            end
            QS.setTime(tp)
            QS.getCurrentData(adjustIV)
            for qq = 1:length(QS.currentData.IV)
                if dim == 1
                    im = squeeze(max(QS.currentData.IV{qq}(pages, :, :), [], dim)) ;
                elseif dim == 2
                elseif dim == 3
                else
                    error('dim > 3 not understood')
                end
            end
        end
        
        [acom,pcom,dcom] = computeAPDVCoords(QS, opts)
        
        function [acom_sm, pcom_sm] = getAPCOMSm(QS) 
            % Load the anterior and posterior 'centers of mass' ie the
            % endpoints of the object's centerline
            try
                acom_sm = h5read(QS.fileName.apdv, '/acom_sm') ;
                pcom_sm = h5read(QS.fileName.apdv, '/pcom_sm') ;
                assert(size(acom_sm, 1) == length(QS.xp.fileMeta.timePoints))
            catch
                opts = load(QS.fileName.apdv_options) ;
                [acom_sm, pcom_sm] = QS.computeAPDCOMs(opts.apdvOpts) ;
            end
        end
        
        function [rot, trans] = getRotTrans(QS)
            % Load the translation to put anterior to origin and AP axis
            % along x axis 
            if ~isempty(QS.APDV.trans)
                % return from self
                trans = QS.APDV.trans ;
            else
                % load from disk
                trans = importdata(QS.fileName.trans) ;
                QS.APDV.trans = trans ;
            end
            % Load the rotation from XYZ to APDV coordinates
            if ~isempty(QS.APDV.rot)
                rot = QS.APDV.rot ;
            else
                % Load the rotation matrix
                rot = importdata(QS.fileName.rot) ;
                QS.APDV.rot = rot ;
            end
        end
        
        function [xyzlim_raw, xyzlim_pix, xyzlim_um, xyzlim_um_buff] = ...
                getXYZLims(QS)
            %[raw, pix, um, um_buff] = GETXYZLIMS(QS)
            % Grab each xyzlim from self, otherwise load from disk
            % full resolution pix
            if ~isempty(QS.plotting.xyzlim_raw)
                xyzlim_raw = QS.plotting.xyzlim_raw ;
            else
                try
                    xyzlim_raw = dlmread(QS.fileName.xyzlim_raw, ',', 1, 0) ; 
                    QS.plotting.xyzlim_raw = xyzlim_raw ;
                catch
                    [QS.plotting.xyzlim_raw, QS.plotting.xyzlim_pix, ...
                        QS.plotting.xyzlim_um, ...
                        QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_raw = QS.plotting.xyzlim_raw ;
                end
            end
            % rotated scaled in full resolution pix
            if ~isempty(QS.plotting.xyzlim_pix)
                xyzlim_pix = QS.plotting.xyzlim_pix ;
            else
                try
                    xyzlim_pix = dlmread(QS.fileName.xyzlim_pix, ',', 1, 0) ; 
                    QS.plotting.xyzlim_pix = xyzlim_pix ;
                catch
                    [~, QS.plotting.xyzlim_pix, ...
                        QS.plotting.xyzlim_um, ...
                        QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_pix = QS.plotting.xyzlim_pix ;
                end
            end
            % rotated scaled APDV in micron
            if ~isempty(QS.plotting.xyzlim_um)
                xyzlim_um = QS.plotting.xyzlim_um ;
            else
                try
                    xyzlim_um = dlmread(QS.fileName.xyzlim_um, ',', 1, 0) ;
                    QS.plotting.xyzlim_um = xyzlim_um ;
                catch
                    [~, ~, QS.plotting.xyzlim_um, ...
                        QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_um = QS.plotting.xyzlim_um ;
                end
            end
            % rotated scaled APDV in micron, with padding
            if ~isempty(QS.plotting.xyzlim_um_buff)
                xyzlim_um_buff = QS.plotting.xyzlim_um_buff ;
            else
                try
                    xyzlim_um_buff = dlmread(QS.fileName.xyzlim_um_buff, ',', 1, 0) ;
                    QS.plotting.xyzlim_um_buff = xyzlim_um_buff ;
                catch
                    [~, ~, ~, QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_um_buff = QS.plotting.xyzlim_um_buff ;
                end
            end
        end
        
        function features = getFeatures(QS, varargin)
            %GETFEATURES(QS, varargin)
            %   Load features of the QS object (those specied, or all of 
            %   them). Features include {'folds', 'fold_onset', 'ssmax', 
            %   'ssfold', 'rssmax', 'rssfold'}. 
            if nargin > 1
                for qq=1:length(varargin)
                    if isempty(eval(['QS.features.' varargin{qq}]))
                        disp(['Loading feature: ' varargin{qq}])
                        QS.loadFeatures(varargin{qq})
                    end
                end
            else
                QS.loadFeatures() ;
            end
            if nargout > 0
                features = QS.features ; 
            end
        end
        function loadFeatures(QS, varargin)
            % Load all features stored in QS.features
            % 
            % Parameters
            % ----------
            % varargin : optional string list/cell
            %   which specific features to load
            %
            if nargin > 1
                if any(strcmp(varargin, {'folds', 'fold_onset', ...
                    'ssmax', 'ssfold', 'rssmax', 'rssfold'}))
                    
                    % Load all features relating to folds
                    disp('Loading folding features')
                    load(QS.fileName.fold, 'folds', 'fold_onset', ...
                        'ssmax', 'ssfold', 'rssmax', 'rssfold') ;
                    QS.features.folds = folds ;
                    QS.features.fold_onset = fold_onset ; 
                    QS.features.ssmax = ssmax ; 
                    QS.features.ssfold = ssfold ;
                    QS.features.rssmax = rssmax ;
                    QS.features.rssfold = rssfold ;
                else
                    error('Feature not recognized')
                end
            else
                % Load all features
                load(QS.fileName.fold, 'folds', 'fold_onset', ...
                    'ssmax', 'ssfold', 'rssmax', 'rssfold') ;
                QS.features.folds = folds ;
                QS.features.fold_onset = fold_onset ; 
                QS.features.ssmax = ssmax ; 
                QS.features.ssfold = ssfold ;
                QS.features.rssmax = rssmax ;
                QS.features.rssfold = rssfold ;
            end
        end
        
        function data = loadBioFormats(QS, fullFileName)
            r = bfGetReader(fullFileName);
            r.setSeries(QS.xp.fileMeta.series-1);
            nChannelsUsed = numel(QS.xp.expMeta.channelsUsed);
            if QS.xp.fileMeta.swapZT == 0
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeZ(), r.getSizeT()];
            else
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeT(), r.getSizeZ()];
            end
            debugMsg(2, ['stack size (xyzt) ' num2str(stackSize) '\n']);

            xSize = stackSize(1);
            ySize = stackSize(2);
            zSize = stackSize(3);
            
            % number of channels
            nTimePts = stackSize(4);
            
            data = zeros([ySize xSize zSize nChannelsUsed], 'uint16');
            for i = 1:r.getImageCount()

                ZCTidx = r.getZCTCoords(i-1) + 1;
                
                % in the fused embryo data coming out of the python script,
                % Z and T are swaped. In general this isn't the case, thus
                % introduce a file metaField swapZT
                if QS.xp.fileMeta.swapZT == 0
                    zidx = ZCTidx(1);
                    tidx = ZCTidx(3);
                else 
                    zidx = ZCTidx(3);
                    tidx = ZCTidx(1);
                end
                cidx = ZCTidx(2);

                % see above: if there is only one timepoint all the planes
                % should be read, if there are multiple timepoints, only
                % the correct time should be read
                if nTimePts == 1 || (nTimePts > 1 && this.currentTime == tidx-1)
                    
                    debugMsg(1,'.');
                    if rem(i,80) == 0
                        debugMsg(1,'\n');
                    end

                    dataCidx = find(QS.xp.expMeta.channelsUsed == cidx);
                    if ~isempty(dataCidx)
                        data(:,:, zidx, dataCidx) = bfGetPlane(r, i);
                    else
                        disp('skipping channel and z plane')
                    end
                end
            end
        end
        
        function setDataLimits(QS, tp, adjustlow_pctile, adjusthigh_pctile)
            % Use timepoint (tp) to obtain hard values for intensity limits
            % so that data is rescaled to fixed limits instead of
            % percentile. This is useful to avoid flickering of overall
            % intensity in data in which a few voxels vary a lot in
            % intensity.
            QS.xp.loadTime(tp);
            QS.xp.rescaleStackToUnitAspect();
            IV = QS.xp.stack.image.apply() ;
            try
                assert(adjusthigh_pctile > 0 && adjustlow_pctile < 100)
                assert(adjusthigh_pctile > adjustlow_pctile)
            catch
                error('adjustment values must be 0<=val<=100 and increasing')
            end
            adjustlow = prctile(IV{1}(:), adjustlow_pctile) ;
            adjusthigh = prctile(IV{1}(:), adjusthigh_pctile) ;
            QS.data.adjustlow = adjustlow ;
            QS.data.adjusthigh = adjusthigh ;
        end
        
        function IV = loadCurrentData(QS, adjustIV)
            IV = getCurrentData(QS, adjustIV) ;
        end
        
        function IV = getCurrentData(QS, adjustIV)
            % IV = getCurrentData(QS, adjustIV) 
            % Load/return volumetric intensity data for current timepoint
            % Note: axis order permutation is applied here upon loading and
            % assignment to self (current QS instance). 
            % 
            % Parameters
            % ----------
            % QS : QS class instance (self)
            % adjustIV : bool (default=true)
            %   apply the intensity adjustment stored in
            %   QS.data.adjustlow/high
            %
            % Returns
            % -------
            % IV : X*Y*Z intensities 
            %   volumetric intensity data of current timepoint
            %   Note this is also assigned to self as QS.currentData.IV, so
            %   avoiding a duplicate copy in RAM is useful for large data
            %   by not returning an argument and instead calling
            %   QS.currentData.IV when needed.
            if nargin < 2
                adjustIV = true ;
            end
            if isempty(QS.currentTime)
                error('No currentTime set. Use QuapSlap.setTime(timePoint)')
            end
            if isempty(QS.currentData.IV)
                % Load 3D data for coloring mesh pullback
                QS.xp.loadTime(QS.currentTime);
                QS.xp.rescaleStackToUnitAspect();
                IV = QS.xp.stack.image.apply() ;
                if adjustIV
                    adjustlow = QS.data.adjustlow ;
                    adjusthigh = QS.data.adjusthigh ;
                    IVtmp = QS.adjustIV(IV, adjustlow, adjusthigh) ;
                    if ~all(QS.data.axisOrder == [1 2 3])
                        for ch = 1:length(IVtmp)
                            QS.currentData.IV{ch} = permute(IVtmp{ch}, QS.data.axisOrder) ;
                        end
                    else
                        QS.currentData.IV = IVtmp ;
                    end
                else
                    if ~all(QS.data.axisOrder == [1 2 3])
                        for ch = 1:length(IV)
                            QS.currentData.IV{ch} = permute(IV{ch}, QS.data.axisOrder) ;
                        end
                    else                   
                        QS.currentData.IV = IV  ;
                    end
                end
            end
            if nargout > 0
                IV = QS.currentData.IV ;
            end
        end
        
        function IV = adjustIV(QS, IV, adjustlow, adjusthigh, forceValues)
            % Parameters
            % ----------
            % QS : QuapSlap class instance
            % IV : 3D numeric (optional, loads currentData if empty)  
            %   data to adjust
            % adjustlow : numeric or #channels x 1 numeric
            %   minimum intensity or percent intensity to adjust data
            % adjusthigh : numeric or #channels x 1 numeric
            %   maximum intensity or percent intensity to adjust data
            % forceValues : enforce the adjustlow/high to be intensity
            % values, not percentile, even if < 100. 
            if nargin > 2 
                adjustlow = QS.data.adjustlow ;
                adjusthigh = QS.data.adjusthigh ;
            end
            if nargin < 2 
                if ~isempty(QS.currentData.IV) 
                    IV = QS.currentData.IV ;
                end
            end
            if nargin < 5
                forceValues = false ;
            end
            
            % If only one value of intensity limit is supplied, duplicate 
            % for each channel of IV            
            if numel(adjustlow) == 1 && length(IV) > 1
                adjustlow = adjustlow *  ones(size(IV)) ;
            end
            if numel(adjusthigh) == 1 && length(IV) > 1
                adjusthigh = adjusthigh *  ones(size(IV)) ;
            end

            % custom image intensity adjustment
            if all(adjustlow == 0) && all(adjusthigh == 0) && ~forceValues
                disp('Using default limits for imadjustn')
                for ii = 1:length(IV)
                    IV{ii} = imadjustn(IV{ii});
                end
            elseif all(adjustlow < 100) && all(adjusthigh < 100) && ~forceValues
                disp('Taking custom limits for imadjustn as prctile')
                for ii = 1:length(IV)
                    IVii = IV{ii} ;
                    vlo = double(prctile( IVii(:) , adjustlow(ii) )) ;
                    vhi = double(prctile( IVii(:) , adjusthigh(ii))) ;
                    disp(['  raw --> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow(ii)), '/', num2str(adjusthigh(ii))])
                    
                    % Now as fraction for imadjustn
                    vlo = vlo / double(max(IVii(:))) ;
                    vhi = vhi / double(max(IVii(:))) ;
                    disp(['  frac--> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow(ii)), '/', num2str(adjusthigh(ii))])
                    IV{ii} = imadjustn(IVii, [double(vlo); double(vhi)]) ;
                end
            else
                % adjusthigh is > 100 or forceValues, so interpret as an intensity value
                disp('Taking custom limits for imadjustn as direct intensity limit values')

                for ii = 1:length(IV)
                    IVii = IV{ii} ;
                    vlo = double(adjustlow(ii)) ;
                    vhi = double(adjusthigh(ii)) ;
                    disp(['--> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow(ii)), '/', num2str(adjusthigh(ii))])
                    tmp = (double(IVii) - vlo) / (vhi - vlo) ;
                    tmp(tmp > (vhi - vlo)) = 1.0 ;
                    IV{ii} = uint16(2^16 * tmp) ;
                    % cast(tmp, class(IVii)) ;  
                end
            end
            if nargout > 0
                disp('Attributing to self.currentData.IV')
                QS.currentData.IV = IV ;
                QS.currentData.adjustlow = adjustlow ;
                QS.currentData.adjustlow = adjusthigh ;
            else
                disp('WARNING: returning IV instead of attributing to self')
            end
        end
        
        % Get velocity
        function getCurrentVelocity(QS, varargin)
            if isempty(QS.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            if isempty(varargin) 
                do_all = true ;
            else
                do_all = false ;
            end
            
            no_piv3d = isempty(fieldnames(QS.currentVelocity.piv3d)) ;
            if (do_all || contains(varargin, 'piv3d')) && no_piv3d
                % Load 3D data for piv results
                piv3dfn = QS.fullFileBase.piv3d ;
                load(sprintf(piv3dfn, QS.currentTime), 'piv3dstruct') ;
                QS.currentVelocity.piv3d = piv3dstruct ;
            end
            no_piv3d2x = isempty(fieldnames(QS.currentVelocity.piv3d2x)) ;
            if (do_all || contains(varargin, 'piv3d2x')) && no_piv3d2x
                % Load 3D data for piv results
                piv3d2xfn = QS.fullFileBase.piv3d2x ;
                load(sprintf(piv3d2xfn, QS.currentTime), 'piv3dstruct') ;
                QS.currentVelocity.piv3d2x = piv3dstruct ;
            end
        end
        
        function vel = getCurrentVelocityMultiChannel(QS, varargin)
            % By default, load 1x resolution piv3d for each channel
            if isempty(QS.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            if isempty(varargin) 
                varargin = 'piv3d' ;
                do_all = false ;
            else
                do_all = false ;
            end
            
            no_piv3d = isempty(QS.currentVelocityMultiChannel.piv3d) ;
            if ~no_piv3d
                for qq = 1:length(QS.xp.expMeta.channelsUsed)
                    no_piv3d = no_piv3d || ...
                        isempty(QS.currentVelocityMultiChannel.piv3d{qq}) ;
                end
            end
            if (do_all || contains(varargin, 'piv3d')) && no_piv3d
                % Load 3D data for piv results
                piv3dfn = sprintf(QS.fullFileBase.pivMultiChannel.v3d, QS.currentTime) ;
                load(piv3dfn, 'piv3dstruct') ;
                QS.currentVelocityMultiChannel.piv3d = piv3dstruct ;
            end
            
            no_piv3d2x = isempty(QS.currentVelocityMultiChannel.piv3d2x) ;
            if ~no_piv3d2x
                for qq = 1:length(QS.xp.expMeta.channelsUsed)
                    no_piv3d2x = no_piv3d2x || ...
                        isempty(QS.currentVelocityMultiChannel.piv3d2x{qq}) || ...  
                        isempty(fieldnames(QS.currentVelocityMultiChannel.piv3d{qq})) ;
                end
            end
            if (do_all || contains(varargin, 'piv3d2x')) && no_piv3d2x
                % Load 3D data for piv results
                piv3d2xfn = QS.fullFileBase.pivMultiChannel.v3d2x ;
                load(sprintf(piv3d2xfn, QS.currentTime), 'piv3dstruct') ;
                QS.currentVelocityMultiChannel.piv3d2x = piv3dstruct ;
            end
            
            if nargout
                vel = QS.currentVelocityMultiChannel ;
            end
        end
        
        % APDV methods
        [acom_sm, pcom_sm] = computeAPDCOMs(QS, opts)
        function ars = xyz2APDV(QS, a)
            %ars = xyz2APDV(QS, a)
            %   Transform 3d coords from XYZ data space to APDV coord sys
            [rot, trans] = QS.getRotTrans() ;
            ars = ((rot * a')' + trans) * QS.APDV.resolution ;
            if QS.flipy
                ars(:, 2) = - ars(:, 2) ;
            end
        end
        
        function axyz = APDV2xyz(QS, a)
            %ars = xyz2APDV(QS, a)
            %   Transform 3d coords from APDV coord sys to XYZ data space
            [rot, trans] = QS.getRotTrans() ;
            if QS.flipy
                a(:, 2) = - a(:, 2) ;
            end
            invRot = QS.invertRotation(rot) ;
            preRot = a / QS.APDV.resolution - trans ; 
            axyz = (invRot * preRot')' ;
            % Note: ars = ((rot * axyz')' + trans) * QS.APDV.resolution ;
        end
        
        function daxyz = APDV2dxyz(QS, a)
            %ars = xyz2APDV(QS, a)
            %   Transform 3d vectors from APDV coord sys to XYZ data space
            [rot, trans] = QS.getRotTrans() ;
            if QS.flipy
                a(:, 2) = - a(:, 2) ;
            end
            invRot = QS.invertRotation(rot) ;
            preRot = a / QS.APDV.resolution ; 
            daxyz = (invRot * preRot')' ;
            % Note: ars = ((rot * axyz')' + trans) * QS.APDV.resolution ;
        end
        
        function dars = dx2APDV(QS, da)
            %dars = dx2APDV(QS, da)
            %   Transform 3d difference vector from XYZ data space to APDV 
            %   coord sys
            [rot, ~] = QS.getRotTrans() ;
            dars = ((rot * da')') * QS.APDV.resolution ;
            if QS.flipy
                dars(:, 2) = - dars(:, 2) ;
            end
        end
        function setAPDVCOMOptions(QS, apdvCOMOpts)
            QS.apdvCOMOptions = apdvCOMOpts ;
        end
        function apdvCOMOptions = loadAPDVCOMOptions(QS)
            load(QS.fileName.apdvCOMOptions, 'apdvCOMOptions')
            QS.apdvCOMOptions = apdvCOMOptions ;
        end     
        function apdvCOMOptions = saveAPDVCOMOptions(QS)
            apdvCOMOptions = QS.APDVCOMs.apdvCOMOptions ;
            save(QS.fileName.apdvCOMOptions, 'apdvCOMOptions')
        end
        [rot, trans, xyzlim_raw, xyzlim, xyzlim_um, xyzlim_um_buff] = ...
            alignMeshesAPDV(QS, alignAPDVOpts) 
        
        % Plot aligned meshes for pretty presentation 
        % -- for ex, for optogenetic figures
        plotAlignedMeshesPretty(QS, options)

        % Load raw mesh or alignedMesh (rotated & scaled to APDV)
        function rawMesh = loadCurrentRawMesh(QS)
            meshfn = sprintf(QS.fullFileBase.mesh, QS.currentTime) ;
            rawMesh = read_ply_mod(meshfn) ;
            QS.currentMesh.rawMesh = rawMesh ;
        end
        function rawMesh = getCurrentRawMesh(QS)
            if isempty(QS.currentTime)
                error('First set currentTime')
            end
            if isempty(QS.currentMesh.rawMesh)
                QS.loadCurrentRawMesh() ;
            end
            if nargout > 0
                rawMesh = QS.currentMesh.rawMesh ;
            end
        end
        function alignedMesh = loadCurrentAlignedMesh(QS)
            meshfn = sprintf(QS.fullFileBase.alignedMesh, QS.currentTime) ;
            alignedMesh = read_ply_mod(meshfn) ;
            QS.currentMesh.alignedMesh = alignedMesh ;
        end
        function amesh = getCurrentAlignedMesh(QS)
            if isempty(QS.currentTime)
                error('First set currentTime')
            end
            if isempty(QS.currentMesh.spcutMeshSm)
                QS.loadCurrentAlignedMesh() ;
            end
            if nargout > 0
                amesh = QS.currentMesh.alignedMesh ;
            end
        end
        
        % Masked Data
        generateMaskedData(QS, options)
        alignMaskedDataAPDV(QS)
        
        % TexturePatch and demo mesh visualization
        % QS.visualizeMeshEvolution() gives evolution sequence with
        % orthogonalviews of data planes -- aka morphsnakes demo
        visualizeMeshEvolution(QS, options)
        plotSeriesOnSurfaceTexturePatch(QS, overwrite, metadat, ...
                                        TexturePatchOptions)
        
        % Surface Area and Volume over time
        measureSurfaceAreaVolume(QS, options)
        
        % Centerlines & cylinderMesh
        extractCenterlineSeries(QS, cntrlineOpts)
        function setEndcapOptions(QS, endcapOpts)
            QS.endcapOptions = endcapOpts ;
        end        
        function loadEndcapOptions(QS)
            tmp = load(QS.fileName.endcapOptions, 'endcapOptions');
            QS.endcapOptions = tmp.endcapOptions ;
        end        
        function saveEndcapOptions(QS)
            endcapOptions = QS.endcapOptions ;
            save(QS.fileName.endcapOptions, 'endcapOptions')
        end
        sliceMeshEndcaps(QS, endcapOpts, methodOpts)
        generateCleanCntrlines(QS, idOptions)
        function getCleanCntrlines(QS)
            if isempty(QS.cleanCntrlines)
                try
                    tmp = load(QS.fileName.cleanCntrlines, 'cntrlines') ;
                    QS.cleanCntrlines = tmp.cntrlines ;
                    disp('Loaded clean centerlines from disk')
                catch
                    disp('No clean centerlines on disk, generating...')
                    QS.cleanCntrlines = generateCleanCntrlines(QS, idOptions) ;
                end
            end
        end
        function loadCurrentCylinderMeshlean(QS)
            cylmeshfn = ...
                sprintf( QS.fullFileBase.cylinderMesh, QS.currentTime ) ;
            QS.currentMesh.cylinderMesh = read_ply_mod( cylmeshfn );
        end
        function loadCurrentCylinderMeshClean(QS)
            cylmeshfn = ...
                sprintf( QS.fullFileBase.cylinderMeshClean, QS.currentTime ) ;
            disp(['Loading cylinderMeshClean ' cylmeshfn])
            QS.currentMesh.cylinderMeshClean = read_ply_mod( cylmeshfn );
        end
        
        % RawRicci meshes -- can be used to form cutMeshes
        [rawRicciMesh, rawRicciMu] = ...
            generateRawRicciMeshTimePoint(QS, tp, options)
        
        % cutMesh
        generateCurrentCutMesh(QS, options)
        plotCutPath(QS, cutMesh, cutPath)
        function loadCurrentCutMesh(QS)
            if isempty(QS.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            cutMeshfn = sprintf(QS.fullFileBase.cutMesh, QS.currentTime) ;
            cutPfn = sprintf(QS.fullFileBase.cutPath, QS.currentTime) ;
            tmp = load(cutMeshfn, 'cutMesh') ;
            tmp.cutMesh.v = tmp.cutMesh.v + tmp.cutMesh.vn * QS.normalShift ;
            QS.currentMesh.cutMesh = tmp.cutMesh ;
            try
                QS.currentMesh.cutPath = dlmread(cutPfn, ',', 1, 0) ;
            catch
                debugMsg(1, 'Could not load cutPath, cutMesh is limited\n')
                % Wait, isn't cutP a field of cutMesh?
                tmp.cutMesh.cutP
                error('check this here --> is cutP a field?')
            end
        end
        
        % spcutMesh
        generateCurrentSPCutMesh(QS, cutMesh, overwrite)
        function spcutMesh = getCurrentSPCutMesh(QS)
            if isempty(QS.currentTime)
                error('First set currentTime')
            end
            if isempty(QS.currentMesh.spcutMesh)
                QS.loadCurrentSPCutMesh() ;
            end
            if nargout > 0
                spcutMesh = QS.currentMesh.spcutMesh ;
            end
        end
        function spcutMesh = loadCurrentSPCutMesh(QS)
            if isempty(QS.currentTime)
                error('First set currentTime')
            end
            spcutMeshfn = sprintf(QS.fullFileBase.spcutMesh, QS.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMesh') ;
            QS.currentMesh.spcutMesh = tmp.spcutMesh ;
            if nargout > 0
                spcutMesh = QS.currentMesh.spcutMesh ;
            end
        end        
        function spcutMeshSm = getCurrentSPCutMeshSm(QS)
            if isempty(QS.currentTime)
                error('First set currentTime')
            end
            if isempty(QS.currentMesh.spcutMeshSm)
                QS.loadCurrentSPCutMeshSm() ;
            end
            if nargout > 0
                spcutMeshSm = QS.currentMesh.spcutMeshSm ;
            end
        end
        function loadCurrentSPCutMeshSm(QS)
            spcutMeshfn = sprintf(QS.fullFileBase.spcutMeshSm, QS.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMeshSm') ;
            QS.currentMesh.spcutMeshSm = tmp.spcutMeshSm ;
        end
        function spcutMeshSmRS = getCurrentSPCutMeshSmRS(QS)
            if isempty(QS.currentTime)
                error('First set currentTime')
            end
            if isempty(QS.currentMesh.spcutMeshSmRS)
                QS.loadCurrentSPCutMeshSmRS() ;
            end
            if nargout > 0
                spcutMeshSmRS = QS.currentMesh.spcutMeshSmRS ;
            end
        end
        function mesh = loadCurrentSPCutMeshSmRS(QS)
            spcutMeshfn = sprintf(QS.fullFileBase.spcutMeshSmRS, QS.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMeshSmRS') ;
            QS.currentMesh.spcutMeshSmRS = tmp.spcutMeshSmRS ;
            if nargout > 0
                mesh = QS.currentMesh.spcutMeshSmRS ;
            end
        end
        function spcutMeshSmRSC = getCurrentSPCutMeshSmRSC(QS)
            if isempty(QS.currentTime)
                error('First set currentTime')
            end
            if isempty(QS.currentMesh.spcutMeshSmRSC)
                QS.loadCurrentSPCutMeshSmRSC() ;
            end
            if nargout > 0
                spcutMeshSmRSC = QS.currentMesh.spcutMeshSmRSC ;
            end
        end
        function loadCurrentSPCutMeshSmRSC(QS)
            spcutMeshfn = sprintf(QS.fullFileBase.spcutMeshSmRSC, QS.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMeshSmRSC') ;
            QS.currentMesh.spcutMeshSmRSC = tmp.spcutMeshSmRSC ;
        end
        
        
        function [xyzrs, fieldfaces, tri, baryc] = uv2APDV(QS, uv, coordSys, umax, vmax)
            % Convert (u,v) pullback coordinates to xyz coordinates in the
            % APDV frame (rotated and scaled frame aligned with AP axis
            % along x and DV axis along z)
            if nargin < 4
                umax = 1 ;
            end
            if nargin < 5 
                vmax = 1 ;
            end
            if isempty(coordSys)
                coordSys = 'spsmrs' ;
            end
            if strcmpi(coordSys, 'spsmrs') 
                mesh = QS.getCurrentSPCutMeshSmRS() ;
                mesh.u(:, 1) = mesh.u(:, 1) * umax / max(mesh.u(:, 1)) ;
            else
                error('handle coordSys here')
            end
            
            [xyzrs, fieldfaces, tri, baryc] = interpolate2Dpts_3Dmesh(mesh.f, mesh.u,...
                                        mesh.v, uv) ;
            % If the loaded coordinates are not already in APDV, but 
            % instead in the data coordSys, convert to APDV
            if ~contains(coordSys, 'rs')
                xyzrs = xyz2APDV(xyz) ;
            end                       
        end
        
        % Radii in a coordsys
        function [radii, radius_cline] = getRadii(QS, options)
            % Load radii from disk in specified coordinate system
            %
            % Parameters
            % ----------
            % options : struct with fields
            %   coordSys: str specifier for how to return radii
            %       'spcutMeshSmRSC' --> return nU x (nV-1) array for radii
            %       of vertices in spcutMeshSmRSC embedding
            %       <add other options here>
            if nargin < 2
                options = struct() ;
            end
            if isfield(options, 'coordSys')
                coordSys = options.coordSys ;
            else
                if length(QS.xp.fileMeta.timePoints) > 1
                    coordSys = 'spcutMeshSmRSC' ;
                else
                    coordSys = 'spcutMesh' ;
                end
            end
            if strcmpi(coordSys, 'spcutMeshSmRSC')
                QS.loadCurrentSPCutMeshSmRSC() ;
                radii = QS.currentMesh.spcutMeshSmRSC.radius_um ;
                radius_cline = mean(QS.currentMesh.spcutMeshSmRSC.radius_um, 2) ;
            elseif strcmpi(coordSys, 'spcutMesh')
                QS.loadCurrentSPCutMesh() 
                radii = QS.currentMesh.spcutMesh.radii_from_mean_uniform_rs ;
                radius_cline = mean(QS.currentMesh.spcutMesh.radii_from_mean_uniform_rs, 2) ;
            else
                error(['Have not coded for this coord sys yet:' coordSys])
            end
        end
        
        % t0_for_phi0 (uvprime cutMesh)
        function mesh = getCurrentUVPrimeCutMesh(QS)
            if isempty(QS.currentMesh.uvpcutMesh)
                QS.loadCurrentSPCutMeshSm() ;
            end
            mesh = QS.currentMesh.uvpcutMesh ;
        end
        function loadCurrentUVPrimeCutMesh(QS)
            uvpcutMeshfn = sprintf(QS.fullFileBase.uvpcutMesh, QS.currentTime) ;
            tmp = load(uvpcutMeshfn, 'uvpcutMesh') ;
            QS.currentMesh.uvpcutMesh = tmp.uvpcutMesh ;
        end
        measureUVPrimePathlines(QS, options)
        % Note: measureBeltramiCoefficient() allows uvprime
        % coordSys.
        
        % Ricci flow for (r=log(rho), phi) coordinate system
        function mesh = getCurrentRicciMesh(QS)
            if isempty(QS.currentMesh.ricciMesh)
                try
                    QS.loadCurrentRicciMesh() ;
                    mesh = QS.currentMesh.ricciMesh ;
                catch
                    disp('Ricci mesh not on disk! Generating...')
                    mesh = QS.generateRicciMeshTimePoint(QS.currentTime) ;
                end
            else
                mesh = QS.currentMesh.ricciMesh ;
            end
        end
        function ricciMesh = loadCurrentRicciMesh(QS, maxIter)
            if nargin < 2
                maxIter = 200 ;
            end
            ricciMeshfn = sprintf(QS.fullFileBase.ricciMesh, ...
                maxIter, QS.currentTime) ;
            load(ricciMeshfn, 'ricciMesh') ;
            QS.currentMesh.ricciMesh = ricciMesh ;
        end
        measureRPhiPathlines(QS, options)
        % Note: measureBeltramiCoefficient() allows ricci coordSys
        measureBeltramiCoefficient(QS, options)
        function beltrami = getBeltramiCoefficient(QS, options)
            if nargin < 2
                options = struct() ;
            end
            if isempty(QS.pathlines.beltrami)
                QS.pathlines.beltrami = QS.loadBeltramiCoefficient(options) ;
            end
            if nargout > 0 
                beltrami = QS.pathlines.beltrami ;
            end
        end
        function beltrami = loadBeltramiCoefficient(QS, options)    
            if isfield(options, 't0Pathlines')
                t0Pathlines = options.t0Pathlines ;
            else
                t0Pathlines = QS.t0set() ;
            end
            fn = sprintf(QS.fileName.pathlines.quasiconformal, t0Pathlines) ;
            QS.pathlines.beltrami = load(fn)  ;
            if nargout > 0 
                beltrami = QS.pathlines.beltrami ;
            end
        end
        measureBeltramiCoefficientPullbackToPullback(QS, options)
        
        % Radial indentation for pathlines
        function indentation = getPathlineIndentation(QS, options)
            if nargin < 2
                options = struct() ;
            end
            indentation = measurePathlineIndentation(QS, options) ;
        end
        function indentation = measurePathlineIndentation(QS, options)
            overwrite = false ;
            t0p = QS.t0set() ;
            if isfield(options, 'overwrite')
                overwrite = options.overwrite ;
            end
            if isfield(options, 't0Pathline')
                t0p = options.t0Pathline ;
            end
            indentFn = sprintf(QS.fileName.pathlines.indentation, t0p) ;
            if ~exist(indentFn, 'file') || overwrite 
                radFn = sprintf(QS.fileName.pathlines.radius, t0p) ;
                if ~exist(radFn, 'file')
                    disp(['pathline radii not on disk: ' radFn])
                    QS.measurePullbackPathlines(options) ;
                end
                
                load(radFn, 'vRadiusPathlines')
                rad = vRadiusPathlines.radii ;
                nU = size(rad, 2) ;
                indentation = 0 * rad ;
                rad0 = rad(vRadiusPathlines.tIdx0, :, :) ;
                for tidx = 1:length(QS.xp.fileMeta.timePoints)
                    indentation(tidx, :, :) = -(rad(tidx, :, :) - rad0) ./ rad0 ;
                end
                save(indentFn, 'indentation')
                
                % Plot the indentation as a kymograph
                close all
                figfn = fullfile(sprintf(QS.dir.pathlines.data, ...
                    t0p), 'indentation_kymograph.png') ;
                set(gcf, 'visible', 'off')
                indentAP = mean(indentation, 3) ;
                uspace = linspace(0, 1, nU) ;
                imagesc(uspace, QS.xp.fileMeta.timePoints, indentAP)
                xlabel('ap position, $u''/L$', 'interpreter', 'latex')
                ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
                caxis([-max(abs(indentAP(:))), max(abs(indentAP(:)))])
                colormap blueblackred
                cb = colorbar() ;
                ylabel(cb, 'indentation $\delta r/r_0$', 'interpreter', 'latex')
                saveas(gcf, figfn)
                
                % Plot in 3d
                % load reference mesh and pathline vertices in 3d
                load(sprintf(QS.fileName.pathlines.refMesh, t0p), ...
                    'refMesh') ;
                load(sprintf(QS.fileName.pathlines.v3d, t0p), 'v3dPathlines') ;
                indentDir = sprintf(QS.dir.pathlines.indentation, t0p) ;
                if ~exist(indentDir, 'dir')
                    mkdir(indentDir) 
                end
                [~,~,~,xyzlim] = QS.getXYZLims() ;
                for tidx = 1:size(rad, 1)
                    tp = QS.xp.fileMeta.timePoints(tidx) ;
                    fn = fullfile(indentDir, 'indentation_%06d.png') ;
                    if ~exist(fn, 'file') || overwrite
                        close all 
                        fig = figure('visible', 'off') ;
                        opts = struct() ;
                        opts.fig = fig ;
                        opts.ax = gca ;
                        xx = v3dPathlines.vXrs(tidx, :) ;
                        yy = v3dPathlines.vYrs(tidx, :) ;
                        zz = v3dPathlines.vZrs(tidx, :) ;
                        v3d = [ xx(:), yy(:), zz(:) ] ;
                        indent = indentation(tidx,:,:) ;
                        opts.sscale = 0.5 ;
                        opts.axisOff = false ;
                        opts.label = 'constriction, $\delta r/r_0$' ;
                        opts.ax_position = [0.1141, 0.1100, 0.6803, 0.8150] ;
                        scalarFieldOnSurface(refMesh.f, v3d, indent(:), opts) ;
                        view(0, 0)
                        axis equal
                        axis off
                        xlim(xyzlim(1, :))
                        ylim(xyzlim(2, :))
                        zlim(xyzlim(3, :))
                        sgtitle(['constriction, $t=$', ...
                            sprintf('%03d', (tp-t0p)*QS.timeInterval), ...
                            ' ', QS.timeUnits ], 'interpreter', 'latex')
                        saveas(gcf, sprintf(fn, tp)) ;
                        close all
                    end
                end
                
            else
                load(indentFn, 'indentation')
            end
        end
        
        % Radial indentation for UVprime pathlines
        function indentation = measureUVPrimePathlineIndentation(QS, options)
            overwrite = false ;
            t0p = QS.t0set() ;
            if isfield(options, 'overwrite')
                overwrite = options.overwrite ;
            end
            if isfield(options, 't0Pathline')
                t0p = options.t0Pathline ;
            end
            indentFn = sprintf(QS.fileName.pathlines_uvprime.indentation, t0p) ;
            if ~exist(indentFn, 'file') || overwrite
                radFn = sprintf(QS.fileName.pathlines_uvprime.radius, t0p) ;
                if ~exist(radFn, 'file')
                    disp(['pathline radii not on disk: ' radFn])
                    QS.measureUVPrimePathlines(options) ;
                end
                
                load(radFn, 'vRadiusPathlines')
                rad = vRadiusPathlines.radii ;
                nU = size(rad, 2) ;
                indentation = 0 * rad ;
                rad0 = rad(vRadiusPathlines.tIdx0, :, :) ;
                for tidx = 1:length(QS.xp.fileMeta.timePoints)
                    indentation(tidx, :, :) = -(rad(tidx, :, :) - rad0) ./ rad0 ;
                end
                save(indentFn, 'indentation')
                
                % Plot the indentation as a kymograph
                close all
                figfn = fullfile(sprintf(QS.dir.pathlines_uvprime.data, ...
                    t0p), 'indentation_kymograph.png') ;
                set(gcf, 'visible', 'off')
                indentAP = mean(indentation, 3) ;
                uspace = linspace(0, 1, nU) ;
                imagesc(uspace, QS.xp.fileMeta.timePoints, indentAP)
                xlabel('ap position, $u''/L$', 'interpreter', 'latex')
                ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
                caxis([-max(abs(indentAP(:))), max(abs(indentAP(:)))])
                colormap blueblackred
                cb = colorbar() ;
                ylabel(cb, 'indentation $\delta r/r_0$', 'interpreter', 'latex')
                saveas(gcf, figfn)
                
                % Plot in 3d
                % load reference mesh and pathline vertices in 3d
                load(sprintf(QS.fileName.pathlines_uvprime.refMesh, t0p), ...
                    'refMesh') ;
                load(sprintf(QS.fileName.pathlines_uvprime.v3d, t0p), 'v3dPathlines') ;
                indentDir = sprintf(QS.dir.pathlines_uvprime.indentation, t0p) ;
                if ~exist(indentDir, 'dir')
                    mkdir(indentDir) 
                end
                [~,~,~,xyzlim] = QS.getXYZLims() ;
                for tidx = 1:size(rad, 1)
                    tp = QS.xp.fileMeta.timePoints(tidx) ;
                    fn = fullfile(indentDir, 'indentation_%06d.png') ;
                    if ~exist(fn, 'file') || overwrite
                        close all 
                        fig = figure('visible', 'off') ;
                        opts = struct() ;
                        opts.fig = fig ;
                        opts.ax = gca ;
                        xx = v3dPathlines.vXrs(tidx, :) ;
                        yy = v3dPathlines.vYrs(tidx, :) ;
                        zz = v3dPathlines.vZrs(tidx, :) ;
                        v3d = [ xx(:), yy(:), zz(:) ] ;
                        indent = indentation(tidx,:,:) ;
                        opts.sscale = 0.5 ;
                        opts.axisOff = false ;
                        opts.label = 'constriction, $\delta r/r_0$' ;
                        opts.ax_position = [0.1141, 0.1100, 0.6803, 0.8150] ;
                        scalarFieldOnSurface(refMesh.f, v3d, indent(:), opts) ;
                        view(0, 0)
                        axis equal
                        axis off
                        xlim(xyzlim(1, :))
                        ylim(xyzlim(2, :))
                        zlim(xyzlim(3, :))
                        sgtitle(['constriction, $t=$', sprintf('%03d', tp), ...
                            ' ', QS.timeUnits ], 'interpreter', 'latex')
                        saveas(gcf, sprintf(fn, tp)) ;
                        close all
                    end
                end
                
            else
                load(indentFn, 'indentation')
            end
        end
        
        
        % Pullbacks
        generateCurrentPullbacks(QS, cutMesh, spcutMesh, spcutMeshSm, pbOptions)
        function doubleCoverPullbackImages(QS, options)
            % options : struct with fields
            %   coordsys : ('sp', 'uv', 'up')
            %       coordinate system to make double cover 
            %   overwrite : bool, default=false
            %       whether to overwrite current images on disk
            %   histeq : bool, default=true
            %       perform histogram equilization during pullback
            %       extension
            %   ntiles : int, default=50 
            %       The number of bins in each dimension for histogram equilization for a
            %       square original image. That is, the extended image will have (a_fixed *
            %       ntiles, 2 * ntiles) bins in (x,y).
            %   a_fixed : float, default=QS.a_fixed
            %       The aspect ratio of the pullback image: Lx / Ly       
            
            if nargin > 1
                % unpack options
                if isfield(options, 'coordsys')
                    coordsys = options.coordsys ;
                    options = rmfield(options, 'coordsys') ;
                    if strcmp(coordsys, 'sp')
                        imDir = QS.dir.im_sp ;
                        imDir_e = QS.dir.im_spe ;
                        fn0 = QS.fileBase.im_sp ;
                        ofn = QS.fileBase.im_spe ;
                    elseif strcmp(coordsys, 'spsm') || strcmp(coordsys, 'sp_sm')
                        imDir = QS.dir.im_sp_sm ;
                        imDir_e = QS.dir.im_sp_sme ;
                        fn0 = QS.fileBase.im_sp_sm ;
                        ofn = QS.fileBase.im_sp_sme ;
                    elseif strcmp(coordsys, 'spsm2') || ...
                            strcmp(coordsys, 'sp_sm2') || ...
                            strcmp(coordsys, 'spsmLUT') || ...
                            strcmp(coordsys, 'sp_smLUT')
                        % equalize the histogram in patches of the image
                        options.histeq = true ;
                        imDir = QS.dir.im_sp_sm ;
                        imDir_e = QS.dir.im_sp_smeLUT ;
                        fn0 = QS.fileBase.im_sp_sm ;
                        ofn = QS.fileBase.im_sp_smeLUT ;
                    elseif strcmp(coordsys, 'rsm') || strcmp(coordsys, 'r_sm')
                        imDir = QS.dir.im_r_sm ;
                        imDir_e = QS.dir.im_r_sme ;
                        fn0 = QS.fileBase.im_r_sm ;
                        ofn = QS.fileBase.im_r_sme ;
                    elseif strcmp(coordsys, 'uv')
                        imDir = QS.dir.im_uv ;
                        imDir_e = QS.dir.im_uve ;
                        fn0 = QS.fileBase.im_uv ;
                        ofn = QS.fileBase.im_uv_sme ;
                    elseif strcmp(coordsys, 'up')
                        imDir = QS.dir.im_up ;
                        imDir_e = QS.dir.im_upe ;
                        fn0 = QS.fileBase.im_up ;
                        ofn = QS.fileBase.im_up_e ;
                    elseif strcmp(coordsys, 'uvprime')
                        imDir = QS.dir.im_uvprime ;
                        imDir_e = QS.dir.im_uvprime_e ;
                        fn0 = QS.fileBase.im_uvprime ;
                        ofn = QS.fileBase.im_uvprime_e ;
                    elseif strcmp(coordsys, 'r') || strcmp(coordsys, 'spr') || strcmp(coordsys, 'sphir') || strcmp(coordsys, 'rsp')
                        imDir = QS.dir.im_r ;
                        imDir_e = QS.dir.im_re ;
                        fn0 = QS.fileBase.im_r ;
                        ofn = QS.fileBase.im_re ;
                    else
                        error(['did not recognize coordsys: ' coordsys])
                    end
                else
                    % Default value of coordsys = 'sp' ;
                    imDir = QS.dir.im_sp ;
                    imDir_e = QS.dir.im_spe ;
                    fn0 = QS.fileBase.im_sp ;
                    ofn = QS.fileBase.im_spe ;
                end
                
                % Customization of the paths
                if isfield(options, 'subdir')
                    imDir = fullfile(imDir, options.subdir) ;
                end
                if isfield(options, 'fn0')
                    fn0 = options.fn0 ;
                end
                
                % pack options if missing fields
                if ~isfield(options, 'histeq')
                    % equalize the histogram in patches of the image
                    options.histeq = false ;
                end
                if ~isfield(options, 'a_fixed')
                    % Assign the aspect ratio for histogram equilization
                    options.a_fixed = QS.a_fixed ;
                end
                if ~isfield(options, 'ntiles')
                    % Number of tiles per circumference and per unit ap
                    % length, so that if aspect ratio is two, there will be
                    % 2*ntiles samplings for histogram equilization along
                    % the ap axis
                    options.ntiles = 50 ;
                end
                if ~isfield(options, 'overwrite')
                    options.overwrite = false ;
                end
            else
                % Default options
                % Default value of coordsys = 'sp' ;
                options = struct() ;
                imDir = QS.dir.im_sp ;
                imDir_e = QS.dir.im_spe ;
                fn0 = QS.fileBase.im_sp ;
                ofn = QS.fileBase.im_spe ;
                options.histeq = false ;
                options.a_fixed = QS.a_fixed ;
                options.ntiles = ntiles ;
            end
            options.outFnBase = ofn ;
            extendImages(imDir, imDir_e, fn0, QS.xp.fileMeta.timePoints, options)
            disp(['done ensuring extended tiffs for ' imDir ' in ' imDir_e])
        end
        
        % measure mean hoop centerline, Length, & Writhe of curve
        function [mss, mcline, avgpts] = getCurrentClineDVhoop(QS)
            clineDVhoopBase = QS.fullFileBase.clineDVhoop ;
            % load the centerline for this timepoint
            if isempty(QS.currentTime)
                error('Please first set currentTime')
            end
            fn = sprintf(clineDVhoopBase, QS.currentTime) ;
            disp(['Loading DVhoop centerline from ' fn])
            load(fn, 'mss', 'mcline', 'avgpts')
            
            % Is this needed here? Yes!
            if QS.flipy
                mcline(:, 2) = -mcline(:, 2) ;
                avgpts(:, 2) = -avgpts(:, 2) ;
            end
            QS.currentCline.mss = mss ;
            QS.currentCline.mcline = mcline ;
            QS.currentCline.avgpts = avgpts ;
        end

        [Wr, Wr_density, dWr, Length_t, clines_resampled] = ...
            measureWrithe(QS, options)
        function Length_t = measureLength(QS, options)
            if nargin < 2
                options = struct() ;
            end
            [~,~,~, Length_t, ~] = QS.measureWrithe(options) ;
        end
        
        % folds & lobes
        identifyFolds(QS, options)
        measureFoldRadiiVariance(QS, options)
        [lengths, areas, volumes] = measureLobeDynamics(QS, options)
        plotLobes(QS, options) 
        function plotConstrictionDynamics(QS, overwrite)
            % Plot the location of the constrictions over time along with
            % centerlines over time
            QS.getXYZLims() ;
            QS.getRotTrans() ;
            QS.t0set() ;
            QS.loadFeatures() ;
            tidxMap = nan(max(QS.xp.fileMeta.timePoints), 1) ;
            tidxMap(QS.xp.fileMeta.timePoints) = 1:length(QS.xp.fileMeta.timePoints); 
            
            % Plot motion of avgpts at folds in yz plane over time
            aux_plot_avgptcline_lobes(QS.features.folds, ...
                QS.features.fold_onset, QS.dir.lobe, ...
                QS.uvexten, QS.plotting.save_ims, ...
                overwrite, QS.xp.fileMeta.timePoints,...
                QS.fullFileBase.spcutMesh, QS.fullFileBase.clineDVhoop, ...
                QS.t0, QS.timeInterval, QS.timeUnits, QS.spaceUnits, ...
                tidxMap)
            
            % Plot motion of DVhoop at folds in yz plane over time
            aux_plot_constriction_DVhoops(QS.features.folds, ...
                QS.features.fold_onset, QS.dir.foldHoopIm,...
                QS.uvexten, QS.plotting.save_ims, ...
                overwrite, ...
                QS.xp.fileMeta.timePoints, QS.fullFileBase.spcutMesh, ...
                QS.fullFileBase.alignedMesh, ...
                QS.normalShift, QS.APDV.rot, QS.APDV.trans, QS.APDV.resolution, ...
                QS.plotting.colors, QS.plotting.xyzlim_um_buff, QS.flipy, ...
                QS.t0, QS.timeInterval, QS.timeUnits, QS.spaceUnits)
            
            % Plot shape of DVhoop at folds in locally transverse plane 
            % over time
            if ~exist(fullfile(QS.dir.lobe, 'constriction_hoops_transverse'), 'dir')
                mkdir(fullfile(QS.dir.lobe, 'constriction_hoops_transverse'))
            end
            aux_plot_constriction_DVhoops_transverse(QS.features.folds, ...
                QS.xp.tIdx(QS.t0), ...
                fullfile(QS.dir.lobe, 'constriction_hoops_transverse'),...
                QS.uvexten, QS.plotting.save_ims, ...
                overwrite, ...
                QS.xp.fileMeta.timePoints, QS.fullFileBase.spcutMesh, ...
                QS.fullFileBase.alignedMesh, ...
                QS.normalShift, QS.APDV.rot, QS.APDV.trans, QS.APDV.resolution, ...
                QS.plotting.colors, QS.plotting.xyzlim_um_buff, QS.flipy, ...
                QS.t0, QS.timeInterval, QS.timeUnits, QS.spaceUnits)
        end
        generateFoldCrossSections(QS, options)
        measureFoldCrossSectionDynamics(QS, options)
        
        % Smooth meshes in time
        [v3dsmM, nsmM] = smoothDynamicSPhiMeshes(QS, options) ;
        function plotSPCutMeshSm(QS, options) 
            plotSPCutMeshSmSeriesUtility(QS, 'spcutMeshSm', options)
        end
        function plotSPCutMeshSmRS(QS, options) 
            plotSPCutMeshSmSeriesUtility(QS, 'spcutMeshSmRS', options)
        end
        function plotSPCutMeshSmRSC(QS, options) 
            plotSPCutMeshSmSeriesUtility(QS, 'spcutMeshSmRSC', options)
        end
        function [v3dsmM, nsmM] = loadSPCutMeshSm(QS) 
            timePoints = QS.xp.fileMeta.timePoints ;
            v3dsmM = zeros(length(timePoints), QS.nU*QS.nV, 3);
            nsmM = zeros(length(timePoints), QS.nU*QS.nV, 3) ;
            % Load each mesh into v3dsmM and nsmM    
            for qq = 1:length(timePoints)
                load(sprintf(QS.fullFileBase.spcutMeshSm, ...
                    timePoints(qq)), 'spcutMeshSm') ;
                v3dsmM(qq, :, :) = spcutMeshSm.v ;
                nsmM(qq, :, :) = spcutMeshSm.vn ;
            end
        end
        
        % spcutMeshSm at DoubleRes (2x resolution)
        function generateSPCutMeshSm2x(QS, overwrite)
            % Double the resolution of spcutMeshSm meshes
            %
            if nargin < 2
                overwrite = false ;
            end
            for tp = QS.xp.fileMeta.timePoints 
                sp2xfn = sprintf(QS.fullFileBase.spcutMeshSm2x, tp) ;
                if overwrite || ~exist(sp2xfn, 'file')
                    mesh1x = load(sprintf(QS.fullFileBase.spcutMeshSm, tp),...
                        'spcutMeshSm') ;
                    mesh1x = mesh1x.spcutMeshSm ;
                    spcutMeshSm2x = QS.doubleResolution(mesh1x) ;
                    disp(['saving ' sp2xfn])
                    save(sp2xfn, 'spcutMeshSm2x')
                end
                sp2xfn = sprintf(QS.fullFileBase.spcutMeshSmRS2x, tp) ;
                spC2fn = sprintf(QS.fullFileBase.spcutMeshSmRSC2x, tp) ;
                if overwrite || ~exist(sp2xfn, 'file') || ...
                        ~exist(spC2fn, 'file')
                    mesh1x = load(...
                        sprintf(QS.fullFileBase.spcutMeshSmRS, tp), ...
                        'spcutMeshSmRS') ;
                    mesh1x = mesh1x.spcutMeshSmRS ;
                    [spcutMeshSmRS2x, spcutMeshSmRSC2x] = ...
                        QS.doubleResolution(mesh1x) ;
                    disp(['saving ' sp2xfn])
                    save(sp2xfn, 'spcutMeshSmRS2x')
                    disp(['saving ' spC2fn])
                    save(spC2fn, 'spcutMeshSmRSC2x')
                end
            end
        end
        
        % Mean & Gaussian curvature videos
        measureCurvatures(QS, options)
        
        % Cell segmentation
        generateCellSegmentation2D(QS, options)
        processCorrectedCellSegmentation2D(QS, options)
        generateCellSegmentation3D(QS, options)
        generateCellSegmentationPathlines3D(QS, options)
        function seg2d = getCurrentSegmentation2D(QS, options)
            % Obtain the cell segmentation in 3D pullback space
            if isempty(QS.currentSegmentation.seg2d)
                try
                    QS.currentSegmentation.seg2d = load(sprintf(QS.fullFileBase.segmentation2d, QS.currentTime)) ;
                catch
                    options.timePoints = [QS.currentTime] ;
                    QS.generateCellSegmentation2D(options) ;
                    QS.currentSegmentation.seg2d = load(sprintf(QS.fullFileBase.segmentation2d, QS.currentTime)) ;
                end
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg2d = QS.currentSegmentation.seg2d ;
            end
        end
        function seg2dCorr = getCurrentSegmentation2DCorrected(QS, options)
            % Decide on coordinate system for corrected pullback binary 
            coordSys = 'spsme';
            if nargin > 1
                if isfield(options, 'coordSys')
                    coordSys = options.coordSys ;
                end
            end
                
            % Obtain the cell segmentation in 3D pullback space
            if isempty(QS.currentSegmentation.seg2dCorrected)
                try
                    QS.currentSegmentation.seg2dCorrected = ...
                        load(sprintf(QS.fullFileBase.segmentation2dCorrected, coordSys, QS.currentTime)) ;
                catch
                    options.timePoints = [QS.currentTime] ;
                    QS.processCorrectedCellSegmentation2D(options) ;
                    QS.currentSegmentation.seg2dCorrected = ...
                        load(sprintf(QS.fullFileBase.segmentation2dCorrected, coordSys, QS.currentTime)) ;
                end
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg2dCorr = QS.currentSegmentation.seg2dCorrected ;
                assert(strcmpi(coordSys, seg2dCorr.coordSys))
            end
        end
        function seg3d = getCurrentSegmentation3D(QS, options)
            % Obtain the cell segmentation in 3D pushforward space
            if isempty(QS.currentSegmentation.seg3d)
                try
                    QS.loadCurrentSegmentation3D() ;
                catch
                    options.timePoints = [QS.currentTime] ;
                    QS.generateCellSegmentation3D(options) ;
                    QS.loadCurrentSegmentation3D() ;
                end    
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg3d = QS.currentSegmentation.seg3d ;
            end
        end
        function seg3dCorr = getCurrentSegmentation3DCorrected(QS, options)
            % Obtain the cell segmentation in 3D pullback space
            if isempty(QS.currentSegmentation.seg3dCorrected)
                try
                    QS.currentSegmentation.seg3dCorrected = ...
                        load(sprintf(QS.fullFileBase.segmentation3dCorrected, QS.currentTime)) ;
                catch
                    options.timePoints = [QS.currentTime] ;
                    options.corrected = true ;
                    QS.processCorrectedCellSegmentation3D(options) ;
                    QS.currentSegmentation.seg3dCorrected = ...
                        load(sprintf(QS.fullFileBase.segmentation3dCorrected, QS.currentTime)) ;
                end
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg3dCorr = QS.currentSegmentation.seg3dCorrected ;
            end
        end
        function seg3d = loadCurrentSegmentation3D(QS) 
            QS.currentSegmentation.seg3d = ...
                load(sprintf(QS.fullFileBase.segmentation3d, QS.currentTime)) ;
            if nargout > 0
                seg3d = QS.currentSegmentation.seg3d ;
            end
        end
        function seg3dCorr = loadCurrentSegmentation3DCorrected(QS) 
            QS.currentSegmentation.seg3dCorrected = ...
                load(sprintf(QS.fullFileBase.segmentation3dCorrected, QS.currentTime)) ;
            if nargout > 0
                seg3dCorr = QS.currentSegmentation.seg3dCorrected ;
            end
        end
        % Note: anisotropy is stored in seg3d.quality
        % no need for: measureCellAnisotropy(QS, options)
        plotSegmentationStatisticsLobes(QS, options)
        estimateIntercalationRate(QS, options)
        
        % density of cells -- nuclei or membrane based
        measureCellDensity(QS, nuclei_or_membrane, options)
        function loadCurrentCellDensity(QS)
            if QS.currentData
                disp('Loading from self')
            else
                disp('Loading from disk')
            end
        end
        plotCellDensity(QS, options)
        plotCellDensityKymograph(QS, options)
        
        % spcutMeshStack 
        %   --> similar to spcutMeshSmStack, but no smoothing
        %   --> this is useful for fixed data
        generateSPCutMeshStack(QS, spcutMeshStackOptions)
        
        % spcutMeshSmStack
        generateSPCutMeshSmStack(QS, spcutMeshSmStackOptions)
        measureThickness(QS, thicknessOptions)
        phi0_fit = fitPhiOffsetsViaTexture(QS, uspace_ds_umax, vspace,...
            phi0_init, phi0TextureOpts)
       
        % uvprime cutMeshSm
        generateUVPrimeCutMeshes(QS, options)
        
        % ricci Mesh (truly conformal mapped mesh)
        [ricciMesh, ricciMu] = generateRicciMeshTimePoint(QS, tp, options) 
        
        % spcutMeshSm coordinate system demo
        coordinateSystemDemo(QS)
        
        % flow measurements
        function getPIV(QS, options)
            % Load PIV results and store in QS.piv if not already loaded
            if isempty(fieldnames(QS.piv.raw)) || isempty(QS.piv.Lx) ...
                    || isempty(QS.piv.Lx) 
                % Load raw PIV results
                if nargin > 1
                    QS.loadPIV(options)
                else
                    QS.loadPIV() 
                end
            end
            
            if isempty(fieldnames(QS.piv.smoothed)) 
                % Additionally smooth the piv output by sigma
                if QS.piv.smoothing_sigma > 0
                    piv = QS.piv.raw ;
                    disp(['Smoothing piv output with sigma=' num2str(QS.piv.smoothing_sigma)])
                    for tidx = 1:length(QS.xp.fileMeta.timePoints)-1
                        velx = piv.u_filtered{tidx} ;
                        vely = piv.v_filtered{tidx} ;
                        piv.u_filtered{tidx} = imgaussfilt(velx, QS.piv.smoothing_sigma) ;
                        piv.v_filtered{tidx} = imgaussfilt(vely, QS.piv.smoothing_sigma) ;
                    end
                    disp('done smoothing')
                    QS.piv.smoothed = piv ;
                end
            end
            
        end
        function piv = loadPIV(QS, options)
            % Load PIV results from disk and store in QS.piv
            if ~isempty(fieldnames(QS.piv.raw)) && isempty(QS.piv.Lx) ...
                    && isempty(QS.piv.Lx) 
                disp("WARNING: Overwriting QS.piv with piv from disk")
            end
            QS.piv.raw = load(QS.fileName.pivRaw.raw) ;  
            timePoints = QS.xp.fileMeta.timePoints ;
            if strcmp(QS.piv.imCoords, 'sp_sme')
                im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, ...
                    timePoints(1))) ;
                % for now assume all images are the same size
                QS.piv.Lx = size(im0, 1) * ones(length(timePoints), 1) ;
                QS.piv.Ly = size(im0, 2) * ones(length(timePoints), 1) ;
            else
                error(['Unrecognized imCoords: ' QS.piv.imCoords])
            end
            if nargout > 0
                piv = QS.piv ;
            end
        end
        measurePIV3d(QS, options)
        
        % Multi-channel PIV --> one struct for each channel
        function piv = getPIVMultiChannel(QS, options)
            % Load PIV results and store in QS.piv if not already loaded
            do_load = false ;
            for qq = 1:length(QS.xp.expMeta.channelsUsed)
                ch = QS.xp.expMeta.channelsUsed(qq) ;
                if isempty(QS.pivMultiChannel.raw) || ...
                        isempty(fieldnames(QS.pivMultiChannel.raw{qq})) || ...
                        isempty(QS.pivMultiChannel.Lx) || ...
                        isempty(QS.pivMultiChannel.Ly) 
                    do_load = do_load || true ;
                    disp(['ch = ' num2str(ch)])
                end
            end
            if do_load
                % Load raw PIV results
                if nargin > 1
                    QS.loadPIVMultiChannel(options)
                else
                    QS.loadPIVMultiChannel() 
                end
            end
            
            if isempty(QS.pivMultiChannel.smoothed) || ...
                    isempty(fieldnames(QS.pivMultiChannel.smoothed{1})) 
                for qq = 1:length(QS.xp.expMeta.channelsUsed)
                    % Additionally smooth the piv output by sigma
                    if QS.pivMultiChannel.smoothing_sigma > 0
                        piv = QS.pivMultiChannel.raw{qq} ;
                        disp(['Smoothing piv output with sigma=' num2str(QS.piv.smoothing_sigma)])
                        for tidx = 1:length(QS.xp.fileMeta.timePoints)-1
                            velx = piv.u_filtered{tidx} ;
                            vely = piv.v_filtered{tidx} ;
                            piv.u_filtered{tidx} = imgaussfilt(velx, QS.piv.smoothing_sigma) ;
                            piv.v_filtered{tidx} = imgaussfilt(vely, QS.piv.smoothing_sigma) ;
                        end
                        disp('done smoothing')
                        QS.pivMultiChannel.smoothed{qq} = piv ;
                    end
                end
            end
            if nargout > 0 
                piv = QS.pivMultiChannel ;
            end
        end
        function piv = loadPIVMultiChannel(QS, options)
            % Load PIV results from disk and store in QS.piv
            if ~isempty(QS.pivMultiChannel.raw)
                if ~isempty(fieldnames(QS.pivMultiChannel.raw{1})) && ...
                    isempty(QS.pivMultiChannel.Lx) ...
                    && isempty(QS.pivMultiChannel.Ly) 
                    
                    disp("WARNING: Overwriting QS.pivMultiChannel with results from disk")
                end
            end
            for ch = QS.xp.expMeta.channelsUsed
                QS.pivMultiChannel.raw{ch} = ...
                    load(sprintf(QS.fileName.pivRawMultiChannel.raw, ch)) ;  
            end
            timePoints = QS.xp.fileMeta.timePoints ;
            if strcmp(QS.pivMultiChannel.imCoords, 'sp_sme')
                im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, ...
                    timePoints(1))) ;
                % for now assume all images are the same size
                QS.pivMultiChannel.Lx = size(im0, 1) * ones(length(timePoints), 1) ;
                QS.pivMultiChannel.Ly = size(im0, 2) * ones(length(timePoints), 1) ;
            else
                error(['Unrecognized imCoords: ' QS.pivMultiChannel.imCoords])
            end
            if nargout > 0
                piv = QS.pivMultiChannel ;
            end
        end
        measurePIV3dMultiChannel(QS, options)
        measurePIV3dDoubleResolution(QS, options)
        % Note: To timeAverage Velocities at Double resolution, pass
        % options.doubleResolution == true
        
        %% Metric
        plotMetric(QS, options) 
        
        %% Pathlines
        function [p2d, p3d] = samplePullbackPathlines(QS, XY0, options)
            %[p2d, p3d] = samplePullbackPathlines(QS, XY0, options)
            % start at XY0, folow flow using barycentric coordinates of PIV
            % pullback pathlines. This is the same as advecting along a
            % fixed Langrangian coordinate. 
            
            % default options
            t0 = QS.t0set() ;
            preview = false ;
            forceShift = NaN ;
            
            % unpack options
            if nargin < 3
                options = struct() ;
            end
            if isfield(options, 't0')
                t0 = options.t0 ;
            elseif isfield(options, 't0Pathlines')
                t0 = options.t0Pathlines ;
            end
            if isfield(options, 'preview')
                preview = options.preview ;
            end
            
            tidx0 = QS.xp.tIdx(t0) ;
            timePoints = QS.xp.fileMeta.timePoints ;
            
            % get image size if we push to 3d or forceShift
            if strcmp(QS.piv.imCoords, 'sp_sme')
                im = imread(sprintf(QS.fullFileBase.im_sp_sme, t0)) ;
                forceShift = -size(im, 1) * 0.5 ;
            else
                error('Handle coordSys here')
            end
            
            pathlines = QS.getPullbackPathlines(t0, 'vertexPathlines', ...
                'vertexPathlines3D') ;
            XX = squeeze(pathlines.vertices.vX(tidx0, :, :)) ;
            YY = squeeze(pathlines.vertices.vY(tidx0, :, :)) ;
            nX = length(unique(XX(:))) ;
            nY = length(unique(YY(:))) ;
            % pivfaces = defineFacesRectilinearGrid(XYgrid, nX, nY) ;
            faces = pathlines.refMesh.f ;
            XY = [XX(:), YY(:)] ;
            cmesh = struct('f', faces, 'u', XY, ...
                'pathPairs', pathlines.refMesh.pathPairs) ;
            tileOpts = struct() ;
            tileOpts.forceShift = forceShift ;
            [tiledFaces, tiledXY] = tileAnnularCutMesh2D(cmesh, ...
                [1, 1], tileOpts) ;
            
            % Note: the following is an exerpt from barycentricMap2d(faces, v2d, vmap, uv)
            % tr0 = triangulation(pivfaces, XY) ;
            % [fieldfaces, baryc0] = pointLocation(tr0, XY0) ; 
            
            p2d = zeros(length(timePoints), size(XY0, 1), 2) ;
            if nargout > 1
                p3d = zeros(length(timePoints), size(XY0, 1), 3) ;
            end
            for tidx = 1:length(timePoints)
                tp = timePoints(tidx) ;
                if mod(tidx, 10) == 0
                    disp(['t=' num2str(tp)])
                end
                % X1 = squeeze(pathlines.piv.XX(tidx, :, :))' ;
                % Y1 = squeeze(pathlines.piv.YY(tidx, :, :))' ;
                X1 = squeeze(pathlines.vertices.vX(tidx, :, :)) ;
                Y1 = squeeze(pathlines.vertices.vY(tidx, :, :)) ;
                XY1 = [X1(:), Y1(:)] ;
                
                
                if nargout > 1
                    % tile current mesh in 2d and 3d
                    v3d1x = squeeze(pathlines.vertices3d.vXrs(tidx, :, :)) ;
                    v3d1y = squeeze(pathlines.vertices3d.vYrs(tidx, :, :)) ;
                    v3d1z = squeeze(pathlines.vertices3d.vZrs(tidx, :, :)) ;
                    
                    cmesh = struct('f', faces, 'u', XY1, ...
                        'v', [v3d1x(:), v3d1y(:), v3d1z(:)], ...
                        'pathPairs', pathlines.refMesh.pathPairs) ;
                    [~, XY1, XYZ1] = tileAnnularCutMesh(cmesh, [1, 1], tileOpts) ;
                    p2d(tidx, :, :) = barycentricMap2d(tiledFaces, tiledXY, XY1, XY0) ;
                    p3d(tidx, :, :) = barycentricMap2d(tiledFaces, tiledXY, XYZ1, XY0) ;
                    assert(~any(isnan(p2d(:))))
                                        
                    if preview && mod(tidx, 10) == 0
                        clf
                        subplot(1, 2, 1)
                        plot(p2d(tidx, :, 1), p2d(tidx, :, 2), '.')
                        axis equal
                        title(['t=' num2str(tp)])
                        subplot(1, 2, 2)
                        plot3(p3d(tidx, :, 1), p3d(tidx, :, 2), ...
                            p3d(tidx, :, 3), '.')
                        axis equal
                        pause(0.1)
                    end
                else
                    % tile current mesh
                    cmesh = struct('f', faces, 'u', XY1, ...
                        'pathPairs', pathlines.refMesh.pathPairs) ;
                    [~, XY1] = tileAnnularCutMesh2D(cmesh, [1, 1]) ;
                    p2d(tidx, :, :) = barycentricMap2d(...
                        tiledFaces, tiledXY, XY1, XY0) ;
                    assert(~any(isnan(p2d(:))))
                end
                
                % if nargout > 1
                %     QS.setTime(tp)
                %     if strcmp(QS.piv.imCoords, 'sp_sme')
                %         mesh = QS.getCurrentSPCutMeshSmRS() ;
                %         umax = max(mesh.u(:, 1)) ;
                %         vmax = 1.0 ;
                %         doubleCovered = true ;
                %         uv = QS.XY2uv(im, squeeze(p2d(tidx, :, :)), ...
                %             doubleCovered, umax, vmax) ;
                %     else
                %         error('Handle this 3d mapping coordSys here')
                %     end
                %     p3d(tidx, :, :) = interpolate2Dpts_3Dmesh(mesh.f, mesh.u, ...
                %         mesh.v, uv) ;
                % 
                %     if preview && mod(tidx, 10) == 0
                %         clf
                %         subplot(1, 2, 1)
                %         plot(p2d(tidx, :, 1), p2d(tidx, :, 2), '.')
                %         axis equal
                %         title(['t=' num2str(tp)])
                %         subplot(1, 2, 2)
                %         plot3(p3d(tidx, :, 1), p3d(tidx, :, 2), ...
                %             p3d(tidx, :, 3), '.')
                %         axis equal
                %         pause(0.1)
                %     end
                % else
                %     if preview && mod(tidx, 10) == 0
                %         clf
                %         plot(p2d(tidx, :, 1), p2d(tidx, :, 2), '.')
                %         axis equal
                %         title(['t=' num2str(tp)])
                %         pause(0.1)
                %     end
                % end
            end
            
            if nargout > 1
                QS.clearTime() ;
            end
        end
        
        measurePullbackPathlines(QS, options)
        function pathlines = getPullbackPathlines(QS, t0, varargin)
            % Discern if we must load pathlines or if already loaded
            doneLoading = false ;
            if nargin > 1 
                if QS.pathlines.t0 ~= t0
                    % The timestamp at which pathlines form grid that is 
                    % requested is different than the one that is loaded,
                    % if any are indeed already loaded. Therefore, we load
                    % anew
                    if nargin > 2
                        % pass varargin along to load method
                        QS.loadPullbackPathlines(t0, varargin)
                    else
                        QS.loadPullbackPathlines(t0)
                    end
                    doneLoading = true ;
                end
            end
            if ~doneLoading
                % No t0 supplied or t0 is the same as what is stored
                % in QS.pathlines.t0, if any is already stored (ie if any
                % pathlines are already loaded)
                if isempty(QS.pathlines.t0)
                    % no pathlines loaded. Load here
                    if nargin > 2
                        QS.loadPullbackPathlines(t0, varargin)
                    elseif nargin > 1
                        QS.loadPullbackPathlines(t0)
                    else
                        QS.loadPullbackPathlines()
                    end
                else
                    % There are pathlines loaded already. Which are
                    % requested here in varargin? First check if varargin 
                    % is empty or not  
                    if nargin > 2            
                        if any(contains(lower(varargin), 'pivpathlines')) || ...
                                any(contains(lower(varargin), 'piv'))
                            if isempty(QS.pathlines.piv)
                                disp('Loading pivPathlines') 
                                QS.loadPullbackPathlines(t0, 'pivPathlines')
                            end          
                        end
                        if any(contains(lower(varargin), 'vertexpathlines')) || ...
                                any(contains(lower(varargin), 'vertex')) || ...
                                any(contains(lower(varargin), 'vertices'))
                            if isempty(QS.pathlines.vertices)
                                disp('Loading vertexPathlines') 
                                QS.loadPullbackPathlines(t0, 'vertexPathlines')
                            end            
                        end            
                        if any(contains(lower(varargin), 'facepathlines')) || ...
                                any(contains(lower(varargin), 'face')) || ...
                                any(contains(lower(varargin), 'faces'))
                            disp('Loading facePathlines') 
                            if isempty(QS.pathlines.faces)
                                QS.loadPullbackPathlines(t0, 'facePathlines')
                            end
                        end
                        if any(contains(lower(varargin), 'vertexpathlines3d')) || ...
                                any(contains(lower(varargin), 'vertex3d')) || ...
                                any(contains(lower(varargin), 'vertices3d'))
                            disp('Loading vertexPathlines3d') 
                            if isempty(QS.pathlines.vertices3d)
                                QS.loadPullbackPathlines(t0, 'vertexPathlines3d')
                            end
                        end
                    else
                        % varargin is not supplied, so load all three if
                        % not already loaded
                        % First grab t0
                        if nargin < 2
                            t0 = QS.t0set() ;
                        end
                        if isempty(QS.pathlines.piv)
                            disp('Loading pivPathlines') 
                            QS.loadPullbackPathlines(t0, 'pivPathlines')
                        end            
                        if isempty(QS.pathlines.vertices)
                            disp('Loading vertexPathlines') 
                            QS.loadPullbackPathlines(t0, 'vertexPathlines')
                        end            
                        if isempty(QS.pathlines.faces)
                            disp('Loading facePathlines') 
                            QS.loadPullbackPathlines(t0, 'facePathlines')
                        end  
                        if isempty(QS.pathlines.vertices3d)
                            disp('Loading vertexPathlines3d') 
                            QS.loadPullbackPathlines(t0, 'vertexPathlines3d')
                        end
                    end
                end
            end
            if nargout > 0
                pathlines = QS.pathlines ;
            end
        end
        function loadPullbackPathlines(QS, t0, varargin)
            if nargin < 2
                t0 = QS.t0set() ;
            elseif isempty(t0)
                t0 = QS.t0set() ;
            else
                try
                    assert(isnumeric(t0))
                catch
                    error('t0 supplied must be numeric')
                end
            end
            % assign t0 as the pathline t0
            QS.pathlines.t0 = t0 ;
            tmp = load(sprintf(QS.fileName.pathlines.refMesh, t0), ...
                'refMesh') ;
            QS.pathlines.refMesh = tmp.refMesh ;
            
            if nargin < 3
                varargin = {'pivPathlines', 'vertexPathlines', ...
                    'facePathlines'} ;
            end
            if any(contains(varargin, 'pivPathlines'))
                load(sprintf(QS.fileName.pathlines.XY, t0), 'pivPathlines')
                QS.pathlines.piv = pivPathlines ;
            end
            if any(contains(varargin, 'vertexPathlines'))
                load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines')
                QS.pathlines.vertices = vertexPathlines ;
            end
            if any(contains(varargin, 'facePathlines'))
                load(sprintf(QS.fileName.pathlines.fXY, t0), 'facePathlines')
                QS.pathlines.faces = facePathlines ;
            end
            if any(contains(varargin, 'vertexPathlines3d'))
                tmp = load(sprintf(QS.fileName.pathlines.v3d, t0)) ;
                QS.pathlines.vertices3d = tmp.v3dPathlines ;
                QS.pathlines.vertices3d.smoothing_sigma = tmp.smoothing_sigma ;
            end
        end
        featureIDs = measurePathlineFeatureIDs(QS, pathlineType, options)
        function featureIDs = getPathlineFeatureIDs(QS, pathlineType, options)
            % featureIDs = GETPATHLINEFEATUREIDS(QS, pathlineType, options)
            %   recall, load, or interactively identify feature locations 
            %   as positions in zeta, the longitudinal pullback coordinate 
            %
            if nargin < 2
                pathlineType = 'vertices' ;
            end
            if nargin < 3
                options = struct() ;
            end
            if strcmpi(pathlineType, 'vertices')
                if isempty(QS.pathlines.featureIDs.vertices)
                    featureIDs = measurePathlineFeatureIDs(QS, ...
                        pathlineType, options) ;
                    QS.pathlines.featureIDs.vertices = featureIDs ;
                else
                    featureIDs = QS.pathlines.featureIDs.vertices ;
                end
            else
                error('Code for this pathlineType here')
            end
        end
        
        function featureIDs = getUVPrimePathlineFeatureIDs(QS, pathlineType, options)
            % featureIDs = getUVPrimePathlineFeatureIDs(QS, pathlineType, options)
            %   recall, load, or interactively identify feature locations 
            %   as positions in zeta, the longitudinal pullback coordinate 
            %
            if nargin < 2
                pathlineType = 'vertices' ;
            end
            if nargin < 3
                options = struct() ;
            end
            if strcmpi(pathlineType, 'vertices')
                if isempty(QS.pathlines_uvprime.featureIDs.vertices)
                    options.field2 = 'radius' ;
                    featureIDs = ...
                        QS.measureUVPrimePathlineFeatureIDs( ...
                        pathlineType, options) ;
                    QS.pathlines.featureIDs.vertices = featureIDs ;
                else
                    featureIDs = QS.pathlines_uvprime.featureIDs.vertices ;
                end
            else
                error('Code for this pathlineType here')
            end            
        end
        
        %% Velocities -- loading Raw / noAveraging
        function loadVelocityRaw(QS, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivRaw.v3d, 'vsmM') ;
                QS.velocityAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivRaw.v2dum, 'v2dsmMum') ;
                QS.velocityAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivRaw.vn, 'vnsmM') ;
                QS.velocityAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivRaw.vf, 'vfsmM') ;
                QS.velocityAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(QS.fileName.pivRaw.vv, 'vvsmM') ;
                QS.velocityAverage.vv = vvsmM ;
            end
        end
        function getVelocityRaw(QS, varargin)
            % todo: check if all varargin are already loaded
            loadVelocityRaw(QS, varargin{:})
        end
        
        %% Velocities -- Lagrangian Averaging
        timeAverageVelocities(QS, samplingResolution, options)
        function loadVelocityAverage(QS, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivAvg.v3d, 'vsmM') ;
                QS.velocityAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivAvg.v2dum, 'v2dsmMum') ;
                QS.velocityAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivAvg.vn, 'vnsmM') ;
                QS.velocityAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivAvg.vf, 'vfsmM') ;
                QS.velocityAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(QS.fileName.pivAvg.vv, 'vvsmM') ;
                QS.velocityAverage.vv = vvsmM ;
            end
        end
        function getVelocityAverage(QS, varargin)
            % todo: check if all varargin are already loaded
            loadVelocityAverage(QS, varargin{:})
        end
        function loadVelocityAverage2x(QS, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivSimAvg2x.v3d, 'vsmM') ;
                QS.velocityAverage2x.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocityAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'v2d'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocityAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivSimAvg2x.vn, 'vnsmM') ;
                QS.velocityAverage2x.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivSimAvg2x.vf, 'vfsmM') ;
                QS.velocityAverage2x.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(QS.fileName.pivSimAvg2x.vv, 'vvsmM') ;
                QS.velocityAverage2x.vv = vvsmM ;
            end
        end
        function getVelocityAverage2x(QS, varargin)
            if isempty(QS.velocityAverage2x.v3d)
                loadVelocityAverage2x(QS, varargin{:})
            end
        end
        plotTimeAvgVelocities(QS, options)
        helmholtzHodge(QS, options)
        measurePathlineVelocities(QS, options)
        plotPathlineVelocities(QS, options)
        
        %% Velocities -- simple/surface-Lagrangian averaging
        timeAverageVelocitiesSimple(QS, samplingResolution, options)
        function loadVelocitySimpleAverage(QS, varargin)
            % Load and pack into struct
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivSimAvg.v3d, 'vsmM') ;
                QS.velocitySimpleAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivSimAvg.v2dum, 'v2dsmMum') ;
                QS.velocitySimpleAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivSimAvg.vn, 'vnsmM') ;
                QS.velocitySimpleAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivSimAvg.vf, 'vfsmM') ;
                QS.velocitySimpleAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'v2v'))
                load(QS.fileName.pivSimAvg.vf, 'vvsmM') ;
                QS.velocitySimpleAverage.vv = vvsmM ;
            end
        end
        function getVelocitySimpleAverage(QS, varargin)
            % todo: check if all varargin are already loaded
            loadVelocitySimpleAverage(QS, varargin{:})
        end
        function loadVelocitySimpleAverage2x(QS, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivSimAvg2x.v3d, 'vsmM') ;
                QS.velocitySimpleAverage2x.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocitySimpleAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'v2d'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocitySimpleAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivSimAvg2x.vn, 'vnsmM') ;
                QS.velocitySimpleAverage2x.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivSimAvg2x.vf, 'vfsmM') ;
                QS.velocitySimpleAverage2x.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(QS.fileName.pivSimAvg2x.vv, 'vvsmM') ;
                QS.velocitySimpleAverage2x.vv = vvsmM ;
            end
        end
        function getVelocitySimpleAverage2x(QS, varargin)
            if isempty(QS.velocitySimpleAverage2x.v3d)
                loadVelocitySimpleAverage2x(QS, varargin{:})
            end
        end
        % NOTE: the following have a simple option for averagingStyle
        % plotTimeAvgVelocities(QS, options)
        % helmholtzHodge(QS, options)
        
        %% compressible/incompressible flow on evolving surface
        measureMetricKinematics(QS, options)
        plotMetricKinematics(QS, options)
        measurePathlineMetricKinematics(QS, options)
        plotPathlineMetricKinematics(QS, options)
        
        %% infer stokes forces
        measureStokesForces(QS, options)
        
        %% Strain RATE
        measureMetricStrainRate(QS, options)
        measureStrainRate(QS, options)
        plotStrainRate(QS, options)
        plotStrainRate3DFiltered(QS, options)
        measurePathlineStrainRate(QS, options)
        measureDxDyStrainFiltered(QS, options)
        % Also makes fund forms in regularlized (zeta, phi) t0 Lagrangian frame
        
        % strain RATE along pathlines
        plotPathlineStrainRate(QS, options)
        
        % Strain by following flow lines
        function strain = getCurrentPathlineStrain(QS, t0Pathlines, varargin)
            %
            % Parameters
            % ----------
            % QS : current class instance
            % t0Pathlines : int or empty
            %   reference timepoint. If empty, set to QS.t0set()
            % varargin : strings ('strain', 'beltrami')
            %   which strain measures to load from disk and attribute to
            %   self
            % 
            % Returns
            % -------
            % strain : QS.currentStrain.pathline
            
            if isempty(QS.currentTime)
                error('Must first set time with QS.setTime()')
            end
            if nargin < 2 
                t0Pathlines = QS.t0set() ;
            else
                if isempty(t0Pathlines)
                    t0Pathlines = QS.t0set() ;
                end
                if nargin < 3
                    varargin = {'strain', 'beltrami'} ;
                end
            end
            
            QS.currentStrain.pathline.t0Pathlines = t0Pathlines ;
            if any(contains(varargin, 'strain'))
                ffn = sprintf(QS.fullFileBase.pathlines.strain, t0Pathlines, QS.currentTime) ;
                QS.currentStrain.pathline.strain = load(ffn) ;
            end

            % Grab all beltramis and assign to self since stored as single
            % file on disk
            if any(contains(varargin, 'beltrami'))
                tidx = QS.xp.tIdx(QS.currentTime) ;
                % stores beltramis 
                QS.getBeltramiCoefficient() ;
                QS.currentStrain.pathline.beltrami.mu_material = ...
                    QS.pathlines.beltrami.mu_material(tidx, :);
                QS.currentStrain.pathline.beltrami.mu_material_filtered = ...
                    QS.pathlines.beltrami.mu_material_filtered(tidx, :);
                QS.currentStrain.pathline.beltrami.mu_material_vertices = ...
                    QS.pathlines.beltrami.mu_material_vertices(tidx, :);
                QS.currentStrain.pathline.beltrami.filterOptions = ...
                    QS.pathlines.beltrami.filterOptions ;
            end
            if nargout > 0
                strain = QS.currentStrain.pathline ;
            end
        end
        measurePathlineStrain(QS, options)
        function getPathlineStrain(QS, options)
            for tidx = 1:length(QS.xp.fileMeta.timePoints)
                tp = QS.xp.fileMeta.timePoints(tidx) ;
                ffn = sprintf(QS.fullFileBase.pathlines.strain, t0Pathlines, tp) ;
                strains(tidx, :) = load(ffn) ;
            end
            QS.pathlines.strains = strains ;
        end
        % Note: plotPathlineStrain also plots pathline beltrami coeffs
        plotPathlineStrain(QS, options)
        function plotPathlineBeltramiKymograph(QS, t0Pathlines, options)
            % Example usage for 2021 gut paper:
            % options = struct('ylim', [0, 2] )
            % QS.plotPathlineBeltramiKymograph([], options)
            % 
            if nargin < 2 
                t0Pathlines = QS.t0set() ;
            else
                if isempty(t0Pathlines)
                    t0Pathlines = QS.t0set() ;
                end
            end 
            
            convert2hrs = false ;
            timeunits = QS.timeUnits ;
            timepoints = QS.xp.fileMeta.timePoints - QS.t0 ;
            if nargin < 3
                options = struct() ;
                if contains(lower(QS.timeUnits), 'min')
                    convert2hrs = length(QS.xp.fileMeta.timePoints) > 60 / QS.timeInterval ;
                    timeunits = 'hr' ;
                    timepoints = timepoints/ 60 ;
                end
            else
                if isfield(options, 'convert2hrs')
                    convert2hrs = options.convert2hrs ;
                elseif contains(lower(QS.timeUnits), 'min')
                    convert2hrs = length(QS.xp.fileMeta.timePoints) > 60 / QS.timeInterval ;
                    timeunits = 'hr' ;
                    timepoints = timepoints/ 60 ;
                end
            end
            
            muAPfn = sprintf(...
                QS.fileName.pathlines.kymographs.mu, t0Pathlines) ;
            tmp = load(muAPfn) ;
            
            re = real(tmp.mu_apM) ;
            im = imag(tmp.mu_apM) ;
            magAP = sqrt(re.^2 + im .^2) ;
            phase = atan2(im, re) ;
            % Create pullback pathline grid as mesh
            mesh = struct() ;
            [uu, vv] = meshgrid(linspace(0,1,QS.nU), timepoints) ;
            mesh.v = [uu(:), vv(:), 0*uu(:)] ;
            pOptions = struct('cbarlabel', '$\mu$', ...
                'mesh', mesh, 'xlabel', 'ap position, $\zeta/L$', ...
                'ylabel', ['time [' timeunits ']']) ;
            if isfield(options, 'ylim') 
                pOptions.ylim = options.ylim ; 
            end
            
            close all
            hf = figure('Position', [100 100 320 320], 'units', 'centimeters');
            [meshHandle, cbs] = plotNematicField(magAP, phase, pOptions)
            axis square 
            set(gcf, 'color','w')
            % set(gcf, 'renderer', 'painters')
            
            saveas(gcf, fullfile(sprintf(...
                QS.dir.pathlines.kymographs, t0Pathlines), ...
                sprintf(...
                'mu_apM_kumograph_pullbackPathlines_%06dt0.png', ...
                t0Pathlines)))
            
            % Save axis content as png image with accompanying pdf of
            % figure frame
            FF = getframe(gca) ;
            imwrite(FF.cdata, fullfile(sprintf(...
                QS.dir.pathlines.kymographs, t0Pathlines), ...
                sprintf(...
                'mu_apM_kumograph_pullbackPathlines_%06dt0_image.png', ...
                t0Pathlines)))
            
            cla ; grid off
            axes(cbs{1})
            % Save figure frame as PDF without image content (for fusion in
            % vector graphics software)
            saveas(gcf, fullfile(sprintf(...
                QS.dir.pathlines.kymographs, t0Pathlines), ...
                sprintf(...
                'mu_apM_kumograph_pullbackPathlines_%06dt0_frame.pdf', ...
                t0Pathlines)))
            
            % Make colorwheel
            
        end
        
        % DEPRICATED
        % compareBeltramiToLinearizedConstriction.m
        
        
        % DEPRICATED -- could integrate rates
        measurePathlineIntegratedStrain(QS, options)
        plotPathlineIntegratedStrain(QS, options)
        
        %% Measurement of relative motion between 2 sets of objects  
        % (for example, multiple layers of the tissue)
        measureRelativeMotionTracks(QS, track1fn, track2fn, Options) 
        plotRelativeMotionTracks2D(QS, Options)
        plotRelativeMotionTracks3D(QS, Options)
        measureRelativeMotionTracksLagrangianFrame(QS, Options)
        
               
        %% Tracking -- manual workflow
        function tracks = manualTrackingAdd(QS, options)
            % Add to current tracks on disk
            trackOutfn = fullfile(QS.dir.tracking, 'manualTracks.mat') ;
            tracks = [] ;
            timePoints = QS.xp.fileMeta.timePoints ;
            imDir = fullfile(QS.dir.im_sp_sm, 'endoderm') ;
            imFileBase = QS.fullFileBase.im_sp_sm ;
            tracks2Add = length(tracks)+1:length(tracks) + 10 ;
                
            if isfield(options, 'trackOutfn')
                trackOutfn = options.trackOutfn ;
            end
            if exist(trackOutfn, 'file')
                load(trackOutfn, 'tracks')
            end
            if isfield(options, 'timePoints')
                timePoints = options.timePoints ;
            end
            if isfield(options, 'imFileBase')
                imFileBase = options.fileBase ;
            end
            if isfield(options, 'tracks2Add')
                tracks2Add = options.tracks2Add ;
            end
            % Decide on tidx0
            if isfield(options, 'tidx0')
                tidx0 = options.tidx0 ;
            elseif all(timePoints == QS.xp.fileMeta.timePoints)
                tidx0 = QS.xp.tIdx(QS.t0set()) ;
            else
                tidx0 = 1 ;
            end
            tracks = manualTrack2D(tracks, imFileBase, timePoints, trackOutfn, tracks2Add, tidx0) ;
        end
        
        function manualTrackingCorrect(options)
            disp('Review tracking results: manualCorrectTracks2D')
            %% Ensure all pairIDs are good tracks in muscle layer
            subdir = 'muscle' ;
            imDir = fullfile(QS.dir.im_sp_sm, subdir) ;
            timePoints = 1:60 ;
            fileBase = fullfile(imDir, QS.fileBase.im_sp_sm) ;
            trackOutfn = fullfile(QS.dir.tracking, 'muscle_tracks.mat') ;
            load(trackOutfn, 'tracks') ;
            tracks2Correct = 1:length(tracks) ;
            if all(timePoints == QS.xp.fileMeta.timePoints)
                tidx = QS.xp.tidx(QS.t0set()) ;
            else
                tidx = 1 ;
            end

            [newTracks, newG] = manualCorrectTracks2D(tracks, fileBase, timePoints, trackOutfn, tracks2Correct) ;
        end
        
        %% Visualize tracked segmentation
        % for visualizing T1 transitions
        visualizeDemoTracks(QS, Options)
        visualizeTracking3D(QS, Options)
        visualizeSegmentationPatch(QS, Options)
        
        %% timepoint-specific coordinate transformations
        sf = interpolateOntoPullbackXY(QS, XY, scalar_field, options)
        
        %% Reconstruction of experiment via NES simulation 
        simulateNES(QS, options)
        
        %% coordSysDemo
        function coordSystemDemo(options)
            % Image for publication/presentation on method & coordinate system
            % Create coordinate system charts visualization using smoothed meshes

            QS.setTime(QS.t0 + 90)
            mesh = QS.loadCurrentSPCutMeshSmRS() ;
            fig = figure('units', 'centimeters', 'position', [0, 0, 13, 13]) ;
            uu = mesh.u(:, 1) ;
            vv = mesh.u(:, 2) ;
            xx = (mesh.v(:, 1)) ;
            yy = (mesh.v(:, 2)) ;
            zz = (mesh.v(:, 3)) ;
            colorsV = viridis(mesh.nV) ;
            colorsU = viridis(mesh.nU) ;

            % LONGITUDE 3D
            subplot(2, 3, [1,2])
            hold off
            for qq = 1:mesh.nV  
                plot3(xx(qq:nU:end), yy(qq:nU:end), zz(qq:nU:end), '-', ...
                    'color', colorsV(qq, :));
                hold on ;
            end
            axis equal
            xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('lateral position [$\mu$m]', 'interpreter', 'latex')
            zlabel('dv position [$\mu$m]', 'interpreter', 'latex')

            % AZIMUTH 3D
            subplot(2, 3, [4,5])
            hold off
            for qq = 1:2:mesh.nU
                inds = (qq-1)*nU+1:qq*nU ;
                plot3(xx(inds), yy(inds), zz(inds), '-', ...
                    'color', colorsV(qq, :));
                hold on ;
            end
            axis equal
            xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('lateral position [$\mu$m]', 'interpreter', 'latex')
            zlabel('dv position [$\mu$m]', 'interpreter', 'latex')

            % LONGITUDE
            subplot(2, 3, 3)
            hold off
            for qq = 1:mesh.nV  
                plot(uu(qq:nU:end), vv(qq:nU:end), '-', ...
                    'color', colorsV(qq, :));
                hold on ;
            end
            axis square
            xlabel('$s$ [$\mu$m]', 'interpreter', 'latex')
            ylabel('$\phi$ [1/$2\pi$]', 'interpreter', 'latex')

            % AZIMUTH
            subplot(2, 3, 6)
            hold off
            for qq = 1:mesh.nV  
                inds = (qq-1)*nU+1:qq*nU ;
                plot(uu(inds), vv(inds), '-', ...
                    'color', colorsV(qq, :));
                hold on ;
            end
            axis square
            xlabel('$s$ [$\mu$m]', 'interpreter', 'latex')
            ylabel('$\phi$ [1/$2\pi$]', 'interpreter', 'latex')

            set(gcf, 'color', 'w')
            export_fig(fullfile(QS.dir.uvCoord, 'coordSystemDemo.png'), '-nocrop','-r600')
        end
    end
    
    methods (Static)
        function uv = XY2uv(im, XY, doubleCovered, umax, vmax)
            %XY2uv(im, XY, doubleCovered, umax, vmax)
            %   Map pixel positions (1, sizeImX) and (1, sizeImY) to
            %   (0, umax) and (0, vmax) of pullback space if singleCover,
            %   or y coords are mapped to (-0.5, 1.5)*vmax if doubleCover
            % 
            % NOTE THAT MAP IS
            % [xesz, yesz] = [size(im, 1), size(im, 2)]
            % uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
            %
            % NOTE THAT INVERSE MAP IS
            % x--> (xy(:, 1) * (Xsz-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (Ysz-1)) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            %
            % Parameters
            % ----------
            % im : NxM numeric array or 2 ints for size
            %   image in which pixel coordinates are defined or dimensions
            %   of the image (pullback image in pixels)
            % XY : Qx2 numeric array
            %   pixel coordinates to convert to pullback space
            % doubleCovered : bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            %
            % Returns
            % -------
            % uv : Qx2 numeric array
            %   pullback coordinates of input pixel coordinates
            
            % Input defaults
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            
            % Input checking
            if size(XY, 2) ~= 2
                if size(XY, 1) == 2
                    XY = XY' ;
                else
                    error('XY must be passed as #pts x 2 numeric array')
                end
            end
            % size of extended image
            if any(size(im) > 2) 
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                % Interpret im as imsize
                Xsz = im(1) ;
                Ysz = im(2) ;
            end
            % map extended image size to (0, 1), (-0.5, 1.5) if double
            % covered. 
            % subtract 1 since pixel positions range from (1, sizeIm)
            xesz = double(Xsz - 1) ;
            yesz = double(Ysz - 1) ;
            % map from pixel y to network y (sphi)
            uv = zeros(size(XY)) ;
            % convert x axis
            uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % convert y axis
            % subtract 1 since pixel positions range from (1, sizeIm)
            if doubleCovered
                uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
            else
                uv(:, 2) = vmax * (XY(:, 2) - 1) / yesz ;
            end
        end
        
        function duv = dXY2duv(im, dXY, doubleCovered, umax, vmax)
            %XY2uv(im, XY, doubleCovered, umax, vmax)
            %   Map difference in pixel positions defined on 
            %   (1, sizeImX) and (1, sizeImY) to
            %   (0, umax) and (0, vmax) of pullback space if singleCover,
            %   or y coords are mapped to (-0.5, 1.5)*vmax if doubleCover
            % 
            % NOTE THAT FULL COORDINATE MAP IS
            % [xesz, yesz] = [size(im, 1), size(im, 2)]
            % uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
            %
            % NOTE THAT INVERSE MAP IS
            % x--> (xy(:, 1) * (Xsz-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (Ysz-1)) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            %
            % Parameters
            % ----------
            % im : NxM numeric array or 2 ints for size
            %   image in which pixel coordinates are defined or dimensions
            %   of the image (pullback image in pixels)
            % dXY : Qx2 numeric array
            %   difference in pixel coordinates as vector, to convert to 
            %   pullback space
            % doubleCovered : bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            %
            % Returns
            % -------
            % uv : Qx2 numeric array
            %   pullback coordinates of input pixel coordinates
            
            % Input defaults
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            
            % Input checking
            if size(dXY, 2) ~= 2
                if size(dXY, 1) == 2
                    dXY = dXY' ;
                elsed
                    error('XY must be passed as #pts x 2 numeric array')
                end
            end
            % size of extended image
            if any(size(im) > 2) 
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                % Interpret im as imsize
                Xsz = im(1) ;
                Ysz = im(2) ;
            end
            % map extended image size to (0, 1), (-0.5, 1.5) if double
            % covered. 
            % subtract 1 since pixel positions range from (1, sizeIm)
            xesz = double(Xsz - 1) ;
            yesz = double(Ysz - 1) ;
            % map from pixel y to network y (sphi)
            duv = zeros(size(dXY)) ;
            % convert x axis
            duv(:, 1) = umax * (dXY(:, 1)) / xesz ;
            % convert y axis
            % WOULD subtract 1 since pixel positions range from (1, sizeIm)
            if doubleCovered
                duv(:, 2) = vmax * 2.0 * (dXY(:, 2)) / yesz - 0.5 ;
            else
                duv(:, 2) = vmax * (dXY(:, 2)) / yesz ;
            end
        end
        
        function XY = uv2XY(im, uv, doubleCovered, umax, vmax) 
            % XY = uv2XY(im, uv, doubleCovered, umax, vmax) 
            %   Map from pullback uv u=(0,1), v=(0,1) to pixel XY
            % x--> (xy(:, 1) * (size(im, 2)-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (size(im, 1)-1)) / (2*vmax) + 0.75 + (size(im,1)-1)*0.25
            %
            % Parameters
            % ----------
            % im : NxM numeric array or length(2) int array
            %   2D image into whose pixel space to map or size(im)
            % uv : Q*2 numeric array
            %   mesh coordinates to convert to pullback pixel space (XY)
            % doubleCovered: bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float 
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            %
            % Returns
            % -------
            % XY : N x 2 float array
            %   positions of uv coordinates in pullback pixel space
            %
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            
            if any(size(im) > 2) 
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                % Interpret im as imsize
                Xsz = im(1) ;
                Ysz = im(2) ;
            end
            XY = 0*uv ;
            XY(:, 1) = uv(:, 1) * (Xsz-1) / (1*umax) + 1 ;
            if doubleCovered
                % image is double cover of physical object (periodic
                % cylinder)
                XY(:, 2) = uv(:, 2) * (Ysz-1) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            else
                % singleCover image of physical cylindrical object
                XY(:, 2) = uv(:, 2) * (Ysz-1) / (1*vmax) + 1  ;
            end        
        end
        
        function [xx, yy] = clipXY(xx, yy, Lx, Ly)
            % Clip x at (1, Lx) and clip Y as periodic (1=Ly, Ly=1), for
            % image that is periodic in Y. Consider Y in [1, Ly]. If the
            % pullback is a doubleCover, Ly = 2*mesh width in pixels
            %
            % Parameters
            % ----------
            % xx : 
            % yy : 
            % Lx : int
            %   number of pixels along x dimension
            % Ly : int
            %   number of pixels along y dimension
            % 
            % Returns
            % -------
            % [xx, yy] : x and y coordinates clipped to [1, Lx] and [1, Ly]
            
            % Note we use minimum values of 1 (in pixels)
            minX = 1 ;
            minY = 1 ;
            
            % Clip in X
            xx(xx > Lx) = Lx ;
            xx(xx < minX ) = 1 ;
            % modulo in Y
            yy(yy > Ly) = yy(yy > Ly) - Ly + minY;
            yy(yy < minY) = yy(yy < minY) + Ly ;
        end
        
        function XY = doubleToSingleCover(XY, Ly)
            % detect if XY is passed as a pair of grids
            if length(size(XY)) > 2 && size(XY, 2) > 2
                % XY is a pair of position grids each as 2D arrays. Clip Y
                tmp = XY(:, :, 2) ;
                tmp(tmp < Ly * .25) = tmp(tmp < Ly * .25) + Ly * 0.5 ;
                tmp(tmp > Ly * .75) = tmp(tmp > Ly * .75) - Ly * 0.5 ;
                XY(:, :, 2) = tmp ;
            elseif size(XY, 2) == 2
                % XY is input as Nx2 array
                XY(XY(:, 2) < Ly * .25, 2) = XY(XY(:, 2) < Ly * .25, 2) + Ly * 0.5 ;
                XY(XY(:, 2) > Ly * .75, 2) = XY(XY(:, 2) > Ly * .75, 2) - Ly * 0.5 ;
            else
                error('XY must be passed as Nx2 or QxRx2 array')
            end
        end
        
        % function uv2pix_old(im, aspect)
        %     % map from network xy to pixel xy
        %     % Note that we need to flip the image (Yscale - stuff) since saved ims had
        %     % natural ydirection.
        %     % Assume here a coord system xy ranging from (0, xscale) and (0, 1) 
        %     % maps to a coord system XY ranging from (0, Yscale * 0.5) and (0, Yscale)
        %     x2Xpix = @(x, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
        %     % y2Ypix = @(y, h, Yscale) Yscale - (Yscale*0.5)*(y+0.5) + h ;
        %     y2Ypix = @(y, h, Yscale) (Yscale*0.5)*(y+0.5) + h ;
        % 
        %     dx2dX = @ (y, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
        %     dy2dY = @ (y, Yscale) (Yscale*0.5)*y ;
        % 
        %     % Now map the coordinates
        % end
        
        [cutMesh, cutMeshC] = doubleResolution(cutMesh, preview)
        
        function [mag_ap, theta_ap] = dvAverageNematic(magnitude, theta)
            %[mag_ap, theta_ap] = DVAVERAGENEMATIC(magnitude, theta)
            % Given a nematic field defined on a rectilinear grid, with 
            % Q(i,j) being in the ith ap position and jth dv position,
            % average the nematic field over the dv positions (dimension 2)
            % 
            % Parameters
            % ----------
            % magnitude : nU x nV numeric array
            %   magnitude of nematic field in 2D rectilinear grid
            % theta : nU x nV float array, with values as angles in radians
            %   angle (mod pi) of nematic field in 2D rectilinear grid
            %
            % Returns
            % -------
            % mag_ap : nU x 1 float array
            %   average magnitude of nematic field along the ap dimension
            % theta_ap : nU x 1 float array
            %   average angle (mod pi) of nematic field along the ap 
            %   dimension
            ap_x = magnitude .* cos(2*theta) ;
            ap_y = magnitude .* sin(2*theta) ;
            ap_xy = [mean(ap_x, 2) , mean(ap_y, 2)] ;
            mag_ap = vecnorm(ap_xy, 2, 2) ;
            theta_averages = atan2(ap_xy(:, 2), ap_xy(:, 1)) ;
            theta_ap = 0.5 * mod(theta_averages, 2*pi) ;
        end
       
        function invRot = invertRotation(rot)
            % Obtain rotation matrix that undoes the APDV rotation 
            rotM = [rot(1, :), 0; rot(2,:), 0; rot(3,:), 0; 0,0,0,1] ;
            tform = affine3d(rotM) ;
            invtform = invert(tform) ;
            invRot = invtform.T(1:3,1:3) ;
        end
        
        function [dorsal, ventral, left, right] = quarterIndicesDV(nV)
            %[dorsal, ventral, left, right] = quarterIndicesDV(nV)
            % indices for each quarter of a DV section in grid coordinates
            if nargin < 1
                nV = QS.nV ;
            end            
            q0 = round(nV * 0.125) ;
            q1 = round(nV * 0.375) ;
            q2 = round(nV * 0.625) ;
            q3 = round(nV * 0.875) ;
            left = q0:q1 ;
            ventral = q1:q2 ;
            right = q2:q3 ;
            dorsal = [q3:nV, 1:q1] ;
        end
    end
    
end