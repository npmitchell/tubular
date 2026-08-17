function ySmoothed = rsavgol(y, order, windowSize, varargin)
% ySmoothed = RSAVGOL(y, order, windowSize)
% ySmoothed = RSAVGOL(y, order, windowSize, 'Name',Value, ...)
%
% Robust Savitzky-Golay filter with median spike suppression. Drop-in
% replacement for SAVGOL. Two stages:
%
%   Stage 1 (median filter): every sample is replaced by the median of a local
%   window of width MedianWin. This is an ordinary sliding-median filter and it
%   bridges ISOLATED spikes -- e.g. the phi0(u) values pinned at 0 where
%   fminbnd failed to match a hoop -- because the real backbone is the window
%   majority. Wide, genuinely-sharp features (the real dips) survive because
%   they span more than half the median window.
%
%   Stage 2 (robust smoothing): a polynomial of the given order is fit to each
%   sliding window of the de-spiked signal by IRLS with Tukey bisquare, giving
%   a smooth backbone. With spikes already gone the IRLS just reduces to a
%   clean least-squares SG fit; it stays as a second line of defence.
%
% OUTPUT ORIENTATION: like SAVGOL, RSAVGOL ALWAYS returns a 1 x N ROW vector,
% independent of input orientation. The pipeline calls it as
% phi0_fit = rsavgol(phi0s, order, win)' -- the transpose expects a row so the
% result is a column; returning a column here would broadcast against an
% nU x 1 accumulator into an nU x nU matrix.
%
% Parameters
% ----------
% y          : 1d array (Nx1 or 1xN)
% order      : polynomial order per smoothing window
% windowSize : odd integer, smoothing window width
%
% Optional name/value
% -------------------
% 'MedianWin' : width of the sliding-median spike filter (default 7). This is
%               the knob that removes spikes. Set larger if spikes come in
%               short CLUSTERS (a cluster wider than MedianWin/2 cannot be
%               bridged by any median); set to 1 to disable stage 1 entirely.
%               NOTE: reachable only if you add it to the rsavgol call inside
%               fitPhiOffsetsFromPrevMesh, which currently passes no varargin,
%               so the default of 7 is what the pipeline uses as-is.
% 'Tune'      : bisquare tuning constant (default 4.685)
% 'MaxIter'   : max IRLS iterations per window (default 10)
% 'Tol'       : relative convergence tol on coefficients (default 1e-4)
%
% Returns
% -------
% ySmoothed : 1 x N row vector
%
% edited NPMitchell 2019; robust IRLS variant 2026; median spike filter + shape fix 2026

%% ---- parameters ----
p = inputParser;
p.addParameter('Tune',      4.685, @(v) isscalar(v) && v > 0);
p.addParameter('MaxIter',   10,    @(v) isscalar(v) && v >= 1);
p.addParameter('Tol',       1e-4,  @(v) isscalar(v) && v > 0);
p.addParameter('MedianWin', 7,     @(v) isscalar(v) && v >= 1);
p.parse(varargin{:});
tune    = p.Results.Tune;
maxIter = p.Results.MaxIter;
tol     = p.Results.Tol;
medwin  = round(p.Results.MedianWin);

y = reshape(y, [1, length(y)]);
N = length(y);

if mod(windowSize,2) == 0 || windowSize < 1
    error('windowSize must be a positive odd integer.');
end
if windowSize <= order
    error('windowSize must be greater than order.');
end

x = y(:);                        % Nx1

%% ---- Stage 1: sliding-median spike filter ----
mHalf = floor((medwin-1)/2);
xmed  = x;
if mHalf > 0
    for i = 1:N
        lo = max(1, i-mHalf);
        hi = min(N, i+mHalf);
        xmed(i) = median(x(lo:hi));
    end
end

%% ---- Stage 2: robust polynomial (SG) fit on the de-spiked signal ----
halfWindow = (windowSize-1)/2;
F          = windowSize;
expo       = 0:order;            % basis t^0 ... t^order
ySmoothed  = zeros(N,1);

for i = 1:N
    % Window placement: centred in the interior, clamped to the first/last F
    % samples near the edges (matching SAVGOL's transient handling).
    if i <= halfWindow
        lo = 1;        hi = F;
    elseif i > N - halfWindow
        lo = N-F+1;    hi = N;
    else
        lo = i-halfWindow;  hi = i+halfWindow;
    end

    % Local coordinate measured from the OUTPUT point i, so the fitted value at
    % i is the constant coefficient c(1) whether or not the window is centred.
    t    = (lo:hi).' - i;        % F x 1
    A    = t .^ expo;            % F x (order+1) Vandermonde
    ywin = xmed(lo:hi);         % F x 1  (de-spiked)

    c = irlsFit(A, ywin, tune, maxIter, tol);
    ySmoothed(i) = c(1);
end

%% ---- return a ROW, matching SAVGOL's contract ----
ySmoothed = reshape(ySmoothed, [1, N]);

end

% ========================================================================
function c = irlsFit(A, y, tune, maxIter, tol)
% IRLS polynomial fit with Tukey bisquare weights.
n = size(A,1);
c = A \ y;                                   % OLS seed
if n <= size(A,2)                            % not over-determined: return OLS
    return;
end
for iter = 1:maxIter
    r = y - A*c;
    s = median(abs(r - median(r))) / 0.6745; % robust scale (MAD)
    if s < eps
        break;                               % residuals ~0
    end
    u = r / (tune * s);
    w = (abs(u) < 1) .* (1 - u.^2).^2;       % bisquare; 0 for gross outliers
    sw   = sqrt(w);
    cNew = (sw .* A) \ (sw .* y);            % weighted LS via row scaling
    if norm(cNew - c) <= tol * max(1, norm(c))
        c = cNew; break;
    end
    c = cNew;
end
end