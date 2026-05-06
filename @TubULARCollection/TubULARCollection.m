classdef TubULARCollection < handle
    % TubULARCollection  Container for same-genotype TubULAR instances with intra-group timing.
    %
    %   A TubULARCollection aggregates several TubULAR datasets that belong to
    %   the SAME genotype (or the same experimental group).  Relative time shifts
    %   are defined only within a collection — i.e., all datasets are aligned to
    %   one chosen reference specimen inside the collection.
    %
    %   To compare two genotypes, build a separate two-dataset collection whose
    %   members are the reference specimen of each genotype, then call
    %   measureRelativeTimeShifts on that cross-genotype collection.
    %
    % Properties
    %   dataDirs    - (N×1 cell) Absolute paths to per-dataset root directories.
    %   genotypes   - (N×1 cell) Genotype label strings matching each dataDir.
    %   outdir      - (char)     Directory for shared output files.
    %   tubis       - (N×1 cell) TubULAR objects, one per dataset.
    %   foldsNorm   - (N×1 cell) Normalised fold positions for each dataset.
    %   radiiInterp - (N×1 cell) Spatially-interpolated radius kymographs.
    %
    % Typical usage
    %   % Within-genotype alignment
    %   wtCol  = TubULARCollection(wtDirs,  repmat({'wt'},3,1),  wtOutdir);
    %   bynCol = TubULARCollection(bynDirs, repmat({'byn'},3,1), bynOutdir);
    %   wtShifts  = wtCol.measureRelativeTimeShifts();
    %   bynShifts = bynCol.measureRelativeTimeShifts();
    %
    %   % Cross-genotype alignment using one reference per genotype
    %   crossCol = TubULARCollection({wtRefDir; bynRefDir}, {'wt';'byn'}, crossOutdir);
    %   crossShifts = crossCol.measureRelativeTimeShifts(struct('refID', 1));
    %
    % See also TubULAR, compute_radius_timestamps,
    %          interpolate_spatially_between_folds

    properties
        dataDirs
        genotypes
        outdir
        tubis
        foldsNorm
        radiiInterp
    end

    methods
        function tc = TubULARCollection(dataDirs, genotypes, outdir)
            % TubULARCollection  Construct a collection from a list of dataset directories.
            %
            %   tc = TubULARCollection(dataDirs) builds the collection using the
            %   directories in the cell array dataDirs.  Each directory must contain
            %   an 'xp.mat' file loadable by TubULAR and a pre-computed strain
            %   output file produced by segment_strains_updated.
            %
            %   tc = TubULARCollection(dataDirs, genotypes) additionally tags each
            %   dataset with a genotype string.  genotypes must be a cell array of
            %   the same length as dataDirs.  Defaults to empty strings.
            %
            %   tc = TubULARCollection(dataDirs, genotypes, outdir) sets the shared
            %   output directory used when saving time-shift results.  When outdir
            %   is omitted, defaults to ./TubULARCollection_<id1>_<id2>_.../ where
            %   each id is the first 8 hex characters of a dataset signature.
            %
            % Parameters
            %   dataDirs  - (N×1 cell of char) Dataset root directories.
            %   genotypes - (N×1 cell of char, optional) Genotype label per dataset.
            %   outdir    - (char, optional) Shared output directory path.
            %
            % Returns
            %   tc - Initialised TubULARCollection handle object.
            if nargin < 1
                error('TubULARCollection:MissingDataDirs', 'dataDirs is required');
            end
            if nargin < 2 || isempty(genotypes)
                genotypes = repmat({''}, numel(dataDirs), 1);
            end
            if nargin < 3
                outdir = TubULARCollection.defaultOutdir(dataDirs);
            end

            tc.dataDirs = dataDirs(:);
            tc.genotypes = genotypes(:);
            tc.outdir = outdir;

            nData = numel(tc.dataDirs);
            tc.tubis = cell(nData, 1);
            tc.foldsNorm = cell(nData, 1);
            tc.radiiInterp = cell(nData, 1);

            tc.initializeTubis();
        end

        function [relativeShifts, timeshifts, optimalShifts, results] = measureRelativeTimeShifts(tc, options)
            % measureRelativeTimeShifts  Compute temporal offsets among datasets in the collection.
            %
            %   All datasets in the collection are treated as a single group.
            %   Every dataset is aligned to one reference specimen (options.refID)
            %   via cross-correlation of spatially-interpolated radius kymographs.
            %
            %   For within-genotype alignment, build a collection from all datasets
            %   of that genotype and call this method.
            %
            %   For cross-genotype alignment, build a two-dataset collection whose
            %   members are the chosen reference specimen of each genotype, set
            %   options.refID = 1, and interpret the returned shift of dataset 2
            %   as the cross-genotype offset.
            %
            % Options struct fields (all optional)
            %   method       - (char) Alignment method.  Currently only
            %                  'radiusKymos' is supported.  Default: 'radiusKymos'.
            %   overwrite    - (logical) Recompute even if a saved result file
            %                  exists.  Default: false.
            %   refID        - (int) Index of the reference dataset within this
            %                  collection.  All other datasets are shifted relative
            %                  to this one.  Default: 1.
            %   saveResults  - (logical) Write outputs to relShiftFile.
            %                  Default: true (false when outdir is empty).
            %   relShiftFile - (char) Full path for the .mat results file.
            %                  Auto-generated from tc.outdir when empty.
            %
            % Returns
            %   relativeShifts - (N×1) Array of time shifts (in timepoints) for
            %                    each dataset relative to the reference.  Entry
            %                    options.refID is always 0.
            %   timeshifts     - Same as relativeShifts (kept for API symmetry).
            %   optimalShifts  - (N×1 cell) Raw cross-correlation peak offsets.
            %   results        - Struct with all of the above plus metadata
            %                    (correlations, shiftValues, dataDirs, genotypes,
            %                    refID, readme).
            %
            % Errors
            %   TubULARCollection:UnknownMethod   - Unsupported options.method.
            %   TubULARCollection:BadRefID        - options.refID out of range.
            %   TubULARCollection:MissingFunction - Helper not on MATLAB path.
            if nargin < 2 || isempty(options)
                options = struct();
            end
            options = tc.applyDefaultOptions(options);

            if options.refID < 1 || options.refID > numel(tc.dataDirs)
                error('TubULARCollection:BadRefID', ...
                    'options.refID (%d) is out of range [1, %d].', ...
                    options.refID, numel(tc.dataDirs));
            end

            if isempty(options.relShiftFile)
                if isempty(tc.outdir)
                    options.saveResults = false;
                else
                    options.relShiftFile = fullfile(tc.outdir, 'timeshifts_fromRadiusKymos.mat');
                end
            end

            if ~isempty(options.relShiftFile) && exist(options.relShiftFile, 'file') && ~options.overwrite
                data = load(options.relShiftFile);
                relativeShifts = data.relativeShifts;
                timeshifts     = data.timeshifts;
                optimalShifts  = data.optimalShifts;
                results = data;
                return
            end

            switch options.method
                case 'radiusKymos'
                    tc.computeRadiiInterp();
                otherwise
                    error('TubULARCollection:UnknownMethod', ...
                        'Unknown options.method: %s', options.method);
            end

            if exist('compute_radius_timestamps', 'file') ~= 2
                error('TubULARCollection:MissingFunction', ...
                    ['compute_radius_timestamps.m is not on the MATLAB path. ', ...
                    'Add it before calling measureRelativeTimeShifts().']);
            end

            [timeshifts, optimalShifts, correlations, shiftValues] = ...
                compute_radius_timestamps(tc.radiiInterp, options.refID, tc.outdir);

            relativeShifts = timeshifts;

            readme = struct();
            readme.description = 'Radius-based intra-collection time-shift alignment';
            readme.relativeShifts = 'Array of shifts (timepoints) per dataset relative to refID; refID entry is 0';
            readme.optimalShifts  = 'Cell array of raw cross-correlation peak offsets, one per dataset';
            readme.dataDirs       = 'Cell array of data directory paths in this collection';
            readme.genotypes      = 'Cell array of genotype labels for each dataset';
            readme.refID          = 'Index of reference dataset within this collection';

            results = struct();
            results.relativeShifts = relativeShifts;
            results.timeshifts     = timeshifts;
            results.optimalShifts  = optimalShifts;
            results.correlations   = correlations;
            results.shiftValues    = shiftValues;
            results.dataDirs       = tc.dataDirs;
            results.genotypes      = tc.genotypes;
            results.refID          = options.refID;
            results.readme         = readme;

            if options.saveResults
                if ~isempty(tc.outdir) && exist(tc.outdir, 'dir') ~= 7
                    mkdir(tc.outdir);
                end
                save(options.relShiftFile, 'relativeShifts', 'timeshifts', 'optimalShifts', 'results');
            end
        end
    end

    methods (Access = private)
        function initializeTubis(tc)
            % initializeTubis  Construct TubULAR objects and load fold data for each dataset.
            %
            %   Called automatically by the constructor.  For every entry in
            %   tc.dataDirs the method:
            %     1. Loads 'xp.mat' and creates a TubULAR instance.
            %     2. Patches file-path fields inside xp and opts to match the
            %        dataset's own directory.
            %     3. Loads normalised fold positions from the DVLRstrain.mat file
            %        produced by segment_strains_updated and stores them in
            %        tc.foldsNorm.
            %
            % Errors
            %   TubULARCollection:MissingFoldData - DVLRstrain.mat not found for a
            %                                       dataset; run segment_strains_updated
            %                                       first.
            for dsetID = 1:numel(tc.dataDirs)
                dataDir = tc.dataDirs{dsetID};
                xpfn = fullfile(dataDir, 'xp.mat');
                load(xpfn, 'xp', 'opts');

                xp.detectOptions.mslsDir = fullfile(dataDir, 'tubular_output');
                xp.fileMeta.dataDir = dataDir;
                opts.meshDir = fullfile(dataDir, 'tubular_output');

                tc.tubis{dsetID} = TubULAR(xp, opts);

                stripesDir = fullfile(dataDir, 'tubular_output', 'cellTracking', 'muscleStripes');
                strainOutfn = fullfile(stripesDir, 'DVLRstrain.mat');
                if exist(strainOutfn, 'file') ~= 2
                    error('TubULARCollection:MissingFoldData', ...
                        'Run segment_strains_updated first for dataset %d', dsetID);
                end
                foldData = load(strainOutfn, 'folds_norm');
                tc.foldsNorm{dsetID} = foldData.folds_norm;
            end
        end

        function computeRadiiInterp(tc)
            % computeRadiiInterp  Build spatially-interpolated radius kymographs.
            %
            %   For each dataset, iterates over all time points, loads the
            %   current SP-cut mesh, and extracts the mean radius profile along
            %   the body-axis coordinate s/L.  The resulting time×s/L matrix is
            %   then passed to interpolate_spatially_between_folds, which
            %   re-samples the profile between the three fold positions so that
            %   the spatial axis is comparable across specimens.
            %
            %   Results are stored in tc.radiiInterp as a cell array, one entry
            %   per dataset.  The method is idempotent; calling it multiple times
            %   simply overwrites existing entries.
            %
            % Errors
            %   TubULARCollection:MissingFunction - interpolate_spatially_between_folds
            %                                       is not on the MATLAB path.
            if exist('interpolate_spatially_between_folds', 'file') ~= 2
                error('TubULARCollection:MissingFunction', ...
                    ['interpolate_spatially_between_folds.m is not on the MATLAB path. ', ...
                    'Add it before calling measureRelativeTimeShifts().']);
            end

            for dsetID = 1:numel(tc.tubis)
                tubi = tc.tubis{dsetID};
                foldsNorm = tc.foldsNorm{dsetID};

                fold1Norm = foldsNorm(:, 1);
                fold2Norm = foldsNorm(:, 2);
                fold3Norm = foldsNorm(:, 3);

                timePoints = tubi.xp.fileMeta.timePoints;
                radii = zeros(length(timePoints), tubi.nU);
                for tidx = 1:length(timePoints)
                    tubi.setTime(timePoints(tidx));
                    mesh = tubi.getCurrentSPCutMesh();
                    radii(tidx, :) = mean(mesh.radii_from_mean_uniform_rs, 2);
                end

                soverL = linspace(0, 1, tubi.nU);
                tc.radiiInterp{dsetID} = interpolate_spatially_between_folds( ...
                    radii, fold1Norm, fold2Norm, fold3Norm, ...
                    'soverL', soverL, 'chamber_sizes', [100, 99, 99, 99]);
            end
        end
    end

    methods (Static, Access = private)
        function outdir = defaultOutdir(dataDirs)
            datasetIDs = cell(numel(dataDirs), 1);
            for qq = 1:numel(dataDirs)
                datasetIDs{qq} = TubULARCollection.datasetSignature(dataDirs{qq});
            end
            outdir = fullfile('.', ['TubULARCollection_' strjoin(datasetIDs, '_')]);
        end

        function signature = datasetSignature(dataDir)
            md = java.security.MessageDigest.getInstance('MD5');
            md.update(uint8(char(dataDir(:).')));
            hash = typecast(md.digest(), 'uint8');
            signature = lower(reshape(dec2hex(hash, 2).', 1, []));
            signature = signature(1:8);
        end

        function options = applyDefaultOptions(options)
            % applyDefaultOptions  Fill missing fields in an options struct with defaults.
            %
            %   options = TubULARCollection.applyDefaultOptions(options) inspects
            %   each recognised field and sets it to its default value when the
            %   field is absent or empty.  Unknown fields are left untouched.
            %
            % Parameters
            %   options - Scalar struct (may be empty or partially populated).
            %
            % Returns
            %   options - Struct with all recognised fields guaranteed to be set.
            %
            % Default values
            %   method       'radiusKymos'
            %   overwrite    false
            %   refID        1
            %   saveResults  true
            %   relShiftFile ''
            if ~isfield(options, 'method') || isempty(options.method)
                options.method = 'radiusKymos';
            end
            if ~isfield(options, 'overwrite') || isempty(options.overwrite)
                options.overwrite = false;
            end
            if ~isfield(options, 'refID') || isempty(options.refID)
                options.refID = 1;
            end
            if ~isfield(options, 'saveResults') || isempty(options.saveResults)
                options.saveResults = true;
            end
            if ~isfield(options, 'relShiftFile') || isempty(options.relShiftFile)
                options.relShiftFile = '';
            end
        end
    end
end