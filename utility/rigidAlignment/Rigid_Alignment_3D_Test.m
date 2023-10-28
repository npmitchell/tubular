%% 3D Rigid Alignment Tests ===============================================
clear; close all; clc;

% Read in target point set
P = readOBJ('horse.obj');

noTrans = false;
noRot = false;

if noTrans
    
    trueT = zeros(3,1);
else
    
    % Generate a random translation
    trueT = 30*rand(3,1)-15;
    
end

if noRot
    
    trueR = eye(3); trueK = [1 0 0];
    trueAA = 0; trueTT = pi/2; truePP = 0;
    
else
    
    % Generate a random rotation
    trueAA = 2*pi*rand(1) - pi;
    trueTT = pi*rand(1);
    truePP = 2*pi*rand(1) - pi;
    
    % Construct the rotation vector
    trueK = [ sin(trueTT) * cos(truePP); ...
        sin(trueTT) * sin(truePP); cos(trueTT) ];
    
    % The cross product matrix of the rotation vector
    CKQ = [ 0, -trueK(3), trueK(2); trueK(3), 0, -trueK(1); ...
        -trueK(2), trueK(1), 0 ];
    
    % Construct the rotation matrix from the axis-angle parameters
    trueR = cos(trueAA) * eye(3) + sin(trueAA) * CKQ + ...
        (1-cos(trueAA)) * (trueK * trueK.');

    clear CKQ
    
end

% Moving point set
% Q = ((trueR.' * P.') - trueT).';
Q = (trueR.' * (P.' - trueT)).';


% Find Alignment Parameters
[R, T, K, aa, tt, pp] = rigidAlignment3D(P, Q, 'DisplayStyle', 'none', ...
    'MaxIterations', 10000, 'NoTranslation', noTrans, 'NoRotation', noRot );
regQ = ((R * Q.') + T).';

fprintf('Rotation Error = %0.5e\n', norm((R-trueR), 'fro'));
fprintf('Translation Error = %0.5e\n', sqrt(sum((T-trueT).^2, 1)));
fprintf('Distance Error = %0.5e\n', max(sqrt(sum((P-regQ).^2, 2))));

scatter3(P(:,1), P(:,2), P(:,3), 'ks');
hold on
scatter3(regQ(:,1), regQ(:,2), regQ(:,3), 'rx');
hold off
axis equal
