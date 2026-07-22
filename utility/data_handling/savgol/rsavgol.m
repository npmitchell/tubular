function ySmoothed = rsavgol(y,order,windowSize,varargin)
% ySmoothed = RSAVGOL(y, order, windowSize)
% ySmoothed = RSAVGOL(y, order, windowSize, 'Tune',T, 'MaxIter',M, 'Tol',t)
%
% Robust Savitzky-Golay filter. Like SAVGOL, a polynomial of the given order
% is fit to each sliding window, but the ordinary least-squares fit is
% replaced by iteratively reweighted least squares (IRLS) with Tukey's
% bisquare loss. Samples whose residual is large relative to a robust scale
% estimate (the MAD) are down-weighted, and gross outliers (residual beyond
% Tune*sigma) are ignored entirely, so a few bad points no longer drag the
% local fit.
%
% Note: the plain SAVGOL precomputes one filter matrix B and applies it with
% FILTER because every window shares the same (identity) weight matrix. Here
% each window gets its own data-dependent weights, so shift-invariance is
% lost and each window is fit explicitly in a loop.
%
% Parameters
% ----------
% y : 1d array of data (Nx1 or 1xN)
% order : order of the polynomial fit to each window
% windowSize : odd integer
%
% Optional name/value parameters
% ------------------------------
% 'Tune' : bisquare tuning constant (default 4.685, ~95% Gaussian efficiency)
% 'MaxIter' : max IRLS iterations per window (default 10)
% 'Tol' : relative convergence tolerance on the coefficients (default 1e-4)
%
% Returns
% -------
% ySmoothed : same shape as y
%
% edited NPMitchell 2019; robust IRLS variant 2026

% ---- parse optional args ----
p = inputParser;
p.addParameter('Tune',    4.685, @(v) isscalar(v) && v > 0);
p.addParameter('MaxIter', 10,    @(v) isscalar(v) && v >= 1);
p.addParameter('Tol',     1e-4,  @(v) isscalar(v) && v > 0);
p.parse(varargin{:});
tune    = p.Results.Tune;
maxIter = p.Results.MaxIter;
tol     = p.Results.Tol;

origShape = size(y);
y = reshape(y, [1,length(y)]);
N = length(y);

if mod(windowSize,2) == 0 || windowSize < 1
    error('windowSize must be a positive odd integer.');
end
if windowSize <= order
    error('windowSize must be greater than order.');
end

halfWindow = (windowSize-1)/2;
F = windowSize;

% Local polynomial basis exponents: columns t^0, t^1, ..., t^order.
% (Same idea as the original S = temp2.^temp3, but built per window below so
% the abscissa is centred on each output point.)
expo = 0:order;

x = y(:);                 % Nx1
ySmoothed = zeros(N,1);

for i = 1:N
    % Window placement: centred in the interior, and clamped to the first or
    % last F samples near the edges (matching the original transient handling).
    if i <= halfWindow
        lo = 1;        hi = F;
    elseif i > N - halfWindow
        lo = N-F+1;    hi = N;
    else
        lo = i-halfWindow;  hi = i+halfWindow;
    end

    % Local coordinate is measured from the OUTPUT point i, so the fitted
    % value at i is simply the constant coefficient c(1) (t^0 term), whether
    % or not the window is centred.
    t     = (lo:hi).' - i;       % F x 1
    A     = t .^ expo;           % F x (order+1) Vandermonde
    ywin  = x(lo:hi);            % F x 1

    c = irlsFit(A, ywin, tune, maxIter, tol);
    ySmoothed(i) = c(1);
end

ySmoothed = reshape(ySmoothed, origShape);

end

% ------------------------------------------------------------------------
function c = irlsFit(A, y, tune, maxIter, tol)
% IRLS polynomial fit with Tukey bisquare weights.

n = size(A,1);

% Seed with the ordinary least-squares fit.
c = A \ y;

% If the system is not over-determined (can happen only for tiny windows),
% robust reweighting is meaningless -- return the OLS solution.
if n <= size(A,2)
    return;
end

for iter = 1:maxIter
    r = y - A*c;

    % Robust scale: MAD about the median, scaled to estimate sigma for
    % Gaussian residuals.
    s = median(abs(r - median(r))) / 0.6745;
    if s < eps
        break;   % residuals ~0: fit is essentially exact, nothing to reweight
    end

    u = r / (tune * s);

    % Tukey bisquare weights: 0 for |u| >= 1, i.e. gross outliers are ignored.
    w = (abs(u) < 1) .* (1 - u.^2).^2;

    % Weighted least-squares update via row scaling (numerically stable).
    sw   = sqrt(w);
    cNew = (sw .* A) \ (sw .* y);

    if norm(cNew - c) <= tol * max(1, norm(c))
        c = cNew;
        break;
    end
    c = cNew;
end

end
