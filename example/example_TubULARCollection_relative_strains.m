%% example_TubULARCollection_relative_strains.m
%
% Example showing how to use TubULARCollection to compute relative time
% shifts among specimens of the same genotype, and then compare the two
% genotype reference specimens against each other.
%
% Design
%   - One TubULARCollection per genotype: relative shifts are defined only
%     within a collection (all specimens vs. one chosen reference).
%   - A separate two-dataset cross-genotype collection holds one reference
%     per genotype; its single time-shift value is the genotype-level offset.
%
% Prerequisites
%   - segment_strains_updated.m must have been run for every dataset so
%     that tubular_output/cellTracking/muscleStripes/DVLRstrain.mat exists.
%   - compute_radius_timestamps.m and
%     interpolate_spatially_between_folds.m must be on the MATLAB path.

addpath(genpath('/mnt/data/code/gptoolbox'))
addpath(genpath('/mnt/data/code/canto/tubular'))
addpath(genpath('/mnt/data/code/tubular/tubular/'))

%% -------------------------------------------------------------------------
%  Dataset paths
% --------------------------------------------------------------------------
wtDirs = { ...
    '/mnt/data/Chris/HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/20240527_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/Deconvolved/All_timepoint_channels/deconvolved_16bit_0810/', ...
    '/mnt/data/Chris/HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/deconvolved/deconvolved_16bit_0810', ...
    '/mnt/data/Chris/HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/2024-06-07_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/deconvolved32bit/deconvolved_16bit_0810/' ...
    };

bynDirs = { ...
    '/mnt/data/Chris/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/20240528_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/Deconvolved_2024-05-28_150855/All_timepoint_channels/deconvolved_16bit_0810/', ...
    '/mnt/data/Chris/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-06-26_175116_crisp_120s_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/deconvolved32bit/deconvolved_16bit_0810', ...
    '/mnt/data/Chris/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-06-28_180134_good_120s_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/unpacked_flipped_rotated_TIFFs_new/deconvolved/deconvolved_16bit_0810/' ...
    };

% Reference specimen index within each genotype collection (1-based).
% All other specimens are shifted relative to this one.
wtRefID  = 2;   % second wt specimen is the reference
bynRefID = 1;   % first bynMyo1c specimen is the reference

outdir     = '/mnt/data/midgut_chirality/analysis/DVLR_strain_comparisons/';
wtOutdir   = fullfile(outdir, 'wt');
bynOutdir  = fullfile(outdir, 'bynMyo1c');
crossOutdir= fullfile(outdir, 'cross_genotype');

overwrite = false;

%% -------------------------------------------------------------------------
%  Step 1: Within-WT relative time shifts
% --------------------------------------------------------------------------
fprintf('========== Step 1: within-WT alignment ==========\n');

wtCol = TubULARCollection(wtDirs, repmat({'wt'}, numel(wtDirs), 1), wtOutdir);

opts_wt.overwrite = overwrite;
opts_wt.refID     = wtRefID;

[wtShifts, ~, ~, wtResults] = wtCol.measureRelativeTimeShifts(opts_wt);

fprintf('WT relative shifts (timepoints, ref = dataset %d):\n', wtRefID);
disp(wtShifts)

%% -------------------------------------------------------------------------
%  Step 2: Within-bynMyo1c relative time shifts
% --------------------------------------------------------------------------
fprintf('========== Step 2: within-bynMyo1c alignment ==========\n');

bynCol = TubULARCollection(bynDirs, repmat({'bynMyo1c'}, numel(bynDirs), 1), bynOutdir);

opts_byn.overwrite = overwrite;
opts_byn.refID     = bynRefID;

[bynShifts, ~, ~, bynResults] = bynCol.measureRelativeTimeShifts(opts_byn);

fprintf('bynMyo1c relative shifts (timepoints, ref = dataset %d):\n', bynRefID);
disp(bynShifts)

%% -------------------------------------------------------------------------
%  Step 3: Cross-genotype shift between the two reference specimens
%
%  Build a two-dataset collection: [wt reference, bynMyo1c reference].
%  The returned shift of dataset 2 relative to dataset 1 is the
%  cross-genotype developmental timing offset.
% --------------------------------------------------------------------------
fprintf('========== Step 3: cross-genotype alignment ==========\n');

crossDirs = { wtDirs{wtRefID}; bynDirs{bynRefID} };
crossCol  = TubULARCollection(crossDirs, {'wt'; 'bynMyo1c'}, crossOutdir);

opts_cross.overwrite = overwrite;
opts_cross.refID     = 1;   % wt reference is the anchor

[crossShifts, ~, ~, crossResults] = crossCol.measureRelativeTimeShifts(opts_cross);

fprintf('Cross-genotype shift: bynMyo1c ref is %+d timepoints relative to wt ref.\n', ...
    crossShifts(2));
