function [ynew, data_output, fft_output] = modeFilter(yy, nmodes, preview)
%[ynew, fft_output] = modeFilter(yy, nmodes, preview)
% FFT reconstruction of lowest nmodes of 1d periodic signal with evenly
% spaced sampling in "time" or other presumably cyclic coordinate.
% Brute force low-pass filter without any fancy Nyquist handling, passband
% ripples, etc. Just select lowest n modes and rebuilt signal from those 
% alone. 
%
% For reference, FFT scalings are:
%     Scale by dt for the FFT, and by Fs for the IFFT
%     Scale by 1/M for the FFT, and by M for the IFFT
%     Scale by 1 for the FFT, and by 1 for the IFFT
%     Scale by 1/sqrt(M) for the FFT, and by sqrt(M) for the IFFT.
%
% Parameters
% ----------
% yy : numeric 1d array
%   input signal to filter
% nmodes : int (default=5)
%   number of modes (including DC offset/average) to use in reconstruction 
% preview : bool (default=false)
%   preview the results in figure form
%
% Returns
% -------
% ynew : numeric 1d array
%   lowpass-filtered output signal reconstructed from lowest nmodes of DFT
%   of input data yy
% data_output : optional output struct with fields
%   Yin : 
%       two-sided DFT of input data, with phases
%   amplitudes :
%       amplitudes used to reconstruct data via 
%       ynew(j) = sum(amp .* cos(2*pi*wn*tt(j)+theta))
%   theta : 
%       phases used to reconstruct data via 
%       ynew(j) = sum(amp .* cos(2*pi*wn*tt(j)+theta)) 
% fft_output : optional output struct with fields
%   Yout : 
%       two-sided DFT of filtered output data, with phases
%   SSBand : 
%       single-sided DFT band of filtered output data
%   magnitudes :
%       magnitudes of DFT of filtered output data
%   readme : struct with same fields as above
%       descriptions of output
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
% modeFilter(yy, 5, true)
%
% NPMitchell 2020

%% Options unpacking
if nargin < 2
    nmodes = 5 ;
    preview = false ;
else
    if nargin < 3
        preview = false ;
    end
end
nV = length(yy) + 1 ;
tt = linspace(0, (nV-1) / nV, nV-1) ;

% if mod(length(yy), 2) == 0
%     error('Input signal might need to be odd in length. Might actually work with even?')
% end

%% Take FFT
YY = fft(yy, nV);  % take nV samples so that we can cleanly div by two
% Create single-sided band
SSB = YY(1:round(nV/2));
% Fold over the second side onto the first side
SSB(2:end) = 2*SSB(2:end);
f = (0:round(nV/2)-1) ; % NOTE: we would multiply this by *(Fs/nV) if sampling frequency Fs were not nV;
% Amplitude -- note normalize the amplitude by L, the length of the signal
magnitudes = abs(SSB/length(tt)) ;

try
    assert(length(magnitudes) > nmodes)
catch
    error('Not enough data for requested number of modes in filter')
end

%% RECONSTRUCTION 
wn = 0:(nmodes-1) ;   %Vector of desired frequencies
amp = zeros(size(wn));
for k = 1:length(wn)
    p = find(f>=wn(k),1);
    
    % magnitude of the mode
    amp(k) = magnitudes(p);
    
    % get angle of this mode component
    teta(k) = atan2(imag(YY(p)), real(YY(p))) ;
end

%% Reconstruct timepoint by timpoint
ynew = zeros(size(tt)) ;
for j = 1:length(tt)
    ynew(j) = sum(amp .* cos(2*pi*wn*tt(j)+teta) );
end

%% FORWARD DISCRET FOURIER TRANSFORM of the RECONSTRUCTED SIGNAL
if nargout > 1 || preview
    Ynew = fft(ynew,nV);
    % Create single-sided band
    SSBnew = Ynew(1:round(nV/2));
    % Fold over the second side onto the first side
    SSBnew(2:end) = 2*SSBnew(2:end);
    % Amplitude -- note normalize the amplitude by L, the length of the signal
    magnitudes_new = abs(SSBnew/length(tt)) ;
end

%% Plot signal
if preview
    clf
    npanels = 3 ;
    subplot(npanels, 1, 1); 
    plot([tt, tt+1], [yy, yy]) ;
    hold on;
    plot([tt, tt+1], [ynew, ynew]); 
    plot([1, 1], [min(yy), max(yy)], 'k--', 'HandleVisibility','off')
    xlabel('time, $t$ [$T$]', 'interpreter', 'latex')
    title(['Lowest ' num2str(nmodes) ' modes of signal'], 'interpreter', 'latex')
    legend({'$y$', '$\tilde{y}$'}, 'interpreter', 'latex')
    ylabel('raw signal $y$ and filtered $\tilde{y}$', 'interpreter', 'latex')
    % Plot fft raw
    subplot(npanels, 1, 2) ;
    plot(0:length(YY)-1, real(YY)) ;
    hold on;
    plot(0:length(YY)-1, imag(YY), '--') ;
    plot(0:length(YY)-1, real(Ynew)) ;
    plot(0:length(YY)-1, imag(Ynew), '--') ;
    legend({'$\Re{Y}$', '$\Im{Y}$', ...
        '$\Re{\tilde{Y}}$', '$\Im{\tilde{Y}}$'}, ...
        'interpreter', 'latex')
    ylabel('FFTs', 'interpreter', 'latex')
    xlabel('mode number', 'interpreter', 'latex')
    subplot(npanels, 1, 3)
    % Plot fft magnitudes
    plot(f, magnitudes)
    hold on;
    plot(f, magnitudes_new)
    xlabel('frequency $f$ [$1/T$]', 'interpreter', 'latex')
    legend({'$Y$', '$\tilde{Y}$'}, 'interpreter', 'latex')
    ylabel('$|$FFT$(y)|$, $|$FFT($\tilde{y})|$',...
        'interpreter', 'latex')
    set(gcf, 'visible', 'on')
end
    
%% Second output
if nargout > 1
    data_output = struct() ;
    data_output.readme = struct() ;
    data_output.readme.Yin = 'two-sided DFT of input data, with phases' ;
    data_output.readme.amplitudes = 'amplitudes used to reconstruct data via ynew(j) = sum(amp .* cos(2*pi*wn*tt(j)+theta))' ;
    data_output.readme.theta = 'phases used to reconstruct data via ynew(j) = sum(amp .* cos(2*pi*wn*tt(j)+theta))' ;
    data_output.Yin = YY ;
    data_output.amplitudes = amp ;
    data_output.theta = teta ;
end

%% Third output
if nargout > 2
    fft_output.Yout = Ynew ;   
    fft_output.SSBand = SSBnew ;
    fft_output.magnitudes = magnitudes_new ;
    fft_output.readme.Yout = 'two-sided DFT of filtered output data, with phases' ;
    fft_output.readme.output_SSBand = 'single-sided DFT band of filtered output data' ;
    fft_output.readme.output_magnitudes = 'magnitudes of DFT of filtered output data' ;
end

    
    
    