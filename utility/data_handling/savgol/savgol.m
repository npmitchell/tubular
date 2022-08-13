function ySmoothed = savgol(y,order,windowSize)
% ySmoothed = SAVGOL(y, order, windowSize)
%
% Parameters
% ----------
% y : 1d array of data (Nx1 or 1xN)
% order : order of the polynomial fit to each window
% windowSize : odd integer
% 
% Returns
% -------
% ySmoothed : N x 1 array of data
% 
% edited NPMitchell 2019

y = reshape(y, [1,length(y)]);

halfWindow = (windowSize-1)/2;

omat = ones(2*halfWindow+1,2*halfWindow+1);
temp1 = tril(omat)*omat;
temp2 = temp1-halfWindow-1;
temp3 = temp1'-1;

% coefficients
s = temp2.^temp3;
S = s(:,1:order+1);

W = eye(windowSize,'double');
[~,R] = qr(sqrt(W)*S,0);

G = S*inv(R)*inv(R)'; 

B = G*S'*W; 

% Preallocate output
x = y'; % x should be Nx1
ySmoothed = zeros(size(x));
F = windowSize;

% Compute the transient on
ySmoothed(1:(F+1)/2-1,:) = flipud(B((F-1)/2+2:end,:))*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B((F-1)./2+1,:),1,x);
ySmoothed((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
ySmoothed(end-(F+1)/2+2:end,:) = flipud(B(1:(F-1)/2,:))*flipud(x(end-(F-1):end,:));

ySmoothed = reshape(ySmoothed,size(y));

end
