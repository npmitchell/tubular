function [outGrid, modeData] = modeFilterQuasi1D(gridData, options)
%[ynew, fft_output] = modeFilter(yy, options)
% FFT reconstruction of lowest nmodes of 1d periodic signal with evenly
% spaced sampling in "time" or other cyclic coordinate, evaluated on a grid
% or field with another dimension (for ex space).
% This is appropriate to bandpass filter the circumferential direction of a
% field living on cylinder with evenly spaced sampling along phi and 
% lines of longitude connected in rows of yy.
%
% Method is a brute force low-pass filter without any fancy Nyquist 
% handling, or passband ripples, or etc. Just select lowest n modes 
% and rebuilt signal from those alone. 
%
% For reference, FFT scalings are:
%     Scale by dt for the FFT, and by Fs for the IFFT
%     Scale by 1/M for the FFT, and by M for the IFFT
%     Scale by 1 for the FFT, and by 1 for the IFFT
%     Scale by 1/sqrt(M) for the FFT, and by sqrt(M) for the IFFT.
%
% Parameters
% ----------
% yy : numeric 2d array
%   input signal to filter spectrally along second dimension (keeping 
%   nmodes) and with tripulse filter for mode coefficients along first 
%   dimension with width widthX
% options : struct with fields
%   nmodesY or nmodes : int (default=5)
%       number of modes (including DC offset/average) to use in reconstruction 
%   widthX or Xwidth or Zwidth: int (default=3)
%       half width of pulse filter along first dimension
%       for ex, if widthX=3, then filter=[0,0.11,0.22,0.33,0.22,0.11,0]
%   extrapolationMethod : 'nearest', 'periodic', 'reflect'
%       (default='nearest')
%       How to handle
%   preview : bool (default=false)
%       preview the results in figure form
%
% Returns
% -------
% ynew : numeric 1d array
%   lowpass-filtered output signal reconstructed from lowest nmodes of DFT
%   of input data yy
% modeData : optional output struct with fields
%   amplitudes : nU x nModes float array
%   phases : nU x nModes float array
%   filter1D : struct with fields
%       Yin : 
%           two-sided DFT of input data, with phases
%       Yout : 
%           two-sided DFT of filtered output data, with phases
%       SSBand : 
%           single-sided DFT band of filtered output data
%       magnitudes :
%           magnitudes of DFT of filtered output data
%       readme : struct with same fields as above
%           descriptions of output
%
% See also
% --------
% https://medium.com/analytics-vidhya/breaking-down-confusions-over-fast-fourier-transform-fft-1561a029b1ab
%
%
% Example usage
% -------------
% nV = 100 ; % note that period T=1
% tt = linspace(0, (nV-1) / nV, nV-1) ;
% yy = sin(6 * pi * tt) + 0.5 * cos(2 * pi * tt + pi/4) + 0.3 * rand(1, nV-1) ;
% opts = struct('nModes', 5, 'preview', true) ;
% modeFilter(yy, opts)
%
% NPMitchell 2020

%% Default options
preview = false ;
nmodesY = 5 ;
widthX = 3 ;
extrapolationMethod = 'nearest' ;

%% Unpack options
if isfield(options, 'nmodesY')
    nmodesY = options.nmodesY ;
elseif isfield(options, 'nModesY')
    nmodesY = options.nModesY ;
elseif isfield(options, 'nmodesy')
    nmodesY = options.nmodesy ;
elseif isfield(options, 'nmodes')
    nmodesY = options.nmodes ;
end
if isfield(options, 'widthX')
    widthX = options.widthX ;
elseif isfield(options, 'widthx')
    widthX = options.widthx ;
elseif isfield(options, 'Xwidth')
    widthX = options.Xwidth ;
elseif isfield(options, 'xwidth')
    widthX = options.xwidth ;
elseif isfield(options, 'Zwidth')
    widthX = options.Zwidth ;
elseif isfield(options, 'zwidth')
    widthX = options.zwidth ;
end
if isfield(options, 'extrapolationMethod')
    extrapolationMethod = options.extrapolationMethod ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end

%% Create theta values for circumferential dimension (dim 2)
tt = linspace(0, 2*pi, size(gridData,2)+1) ;
tt = tt(1:end-1) ;

%% Cycle over each y strip
nU = size(gridData, 1) ;
for qq = 1:nU 
    ystrip = squeeze(gridData(qq, :)) ;
    [ynew, data] = modeFilter(ystrip, nmodesY) ;
    dat(qq, :) = ynew ;
    amps(qq, :) = data.amplitudes ;
    thetas(qq, :) = data.theta ;
    
    % Visualize the constitution of this reconstructed strip
    if preview
        clf
        plot(tt, ystrip, 'o') ;
        hold on;
        plot(tt, ynew, '.') ;
        for pp = 1:length(amps(qq, :))
            if pp == 1
                terms = amps(qq, pp) * cos((pp-1) * tt + thetas(qq, pp)) ;
            else
                terms = terms + amps(qq, pp) * cos((pp-1) * tt + thetas(qq, pp)) ;    
            end
            plot(tt, terms) ;
            pause(1)
        end
    end
end

%% Average in pulse filter along first dimension axis
pulse = [ 1:widthX, fliplr(1:widthX-1) ] ;
pulse = pulse ./ sum(pulse) ;
pmIDs = -(widthX-1):(widthX-1)  ;
assert(length(pulse) == length(pmIDs))

% Consider each strip along first dimension, average modes with phase
% interference
% Initialize output 
outGrid = 0 * gridData ;
if widthX > 0
    for qq = 1:nU
        % for this site along first dim, average each mode with interference
        for modeID = 1:nmodesY
            % for this mode at this site, average with pulse filter

            % First initialize amplitudes of this mode
            ampX = 0 ;
            ampY = 0 ;  
            for pp = 1:length(pulse) 
                % Consider each mode

                % Handle extrapolation if pulse filter extends out of bounds
                % of the range [1,nU]
                if strcmpi(extrapolationMethod, 'nearest')
                    % clip IDs at range bounds [1, nU] 
                    pID = min(max(1, pmIDs(pp)+qq), nU) ;
                elseif strcmpi(extrapolationMethod, 'reflect')
                    pID = pmIDs(pp)+qq ;
                    if pID < 1
                        pID = 1 + abs(1-pID) ;
                    elseif pID > nU
                        pID = nU - abs(pID(pID > nU) - nU) ;
                    end
                elseif strcmpi(extrapolationMethod, 'periodic')
                    pID = mod(pID, nU) ;
                    if pID == 0
                        pID = nU ;
                    end
                end

                theta = thetas(pID, modeID) ;
                amp = amps(pID, modeID) ;
                % DC component does not have interference
                if modeID == 0
                    ampX = pulse ;
                    ampY = 0 ;
                else
                    % higher modes do have phase interference
                    ampX = ampX + amp * cos(theta) * pulse(pp) ;
                    ampY = ampY + amp * sin(theta) * pulse(pp) ;
                end
            end
            ampFinal(qq, modeID) = sqrt(ampX^2 + ampY^2) ;
            thetaFinal(qq, modeID) = atan2(ampY, ampX) ;
        end

        % Add up these modes to reconstruct spectrally-filtered signal
        for modeID = 1:nmodesY
            outGrid(qq, :) = outGrid(qq, :) + ...
                ampFinal(qq, modeID) * cos((modeID-1) * tt + thetaFinal(qq, modeID)) ;
        end

        if preview 
            imagesc(outGrid)
            pause(0.001)
            if qq == 1
                colorbar()
            end
        end
    end
    
    % Supply spectral output if desired
    if nargout > 1
        modeData.amplitudes = ampFinal ;
        modeData.phases = thetaFinal ;
        modeData.filter1D = dat ;
    end
else
    outGrid = dat ;
    % Supply spectral output if desired
    if nargout > 1
        modeData.amplitudes = amps ;
        modeData.phases = thetas ;
        modeData.filter1D = dat ;
    end
    
end

