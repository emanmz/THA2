% pt. j test func 
% clear; clc; close all;
addpath("Functions");

% FR3 Parameters
L = [0.333 0.316 0.384 0.107];
A = 0.088;
M = [1 0 0 A; 0 -1 0 0; 0 0 -1 L(1)+L(2)+L(3)-L(4); 0 0 0 1];

% S_space axes (from the franka file)
ws = {[0;0;1], [0;-1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qs = {[0;0;0], [0;0;L(1)], [0;0;L(1)], [A;0;L(1)+L(2)], [0;0;L(1)+L(2)+L(3)], ...
      [0;0;L(1)+L(2)+L(3)], [A;0;L(1)+L(2)+L(3)-L(4)]};
S_space = zeros(6,7);
for i=1:7
    S_space(:, i) = [ws{i}; cross(-ws{i}, qs{i})];
end

% random reachable poser
theta_test = [0, -pi/4, 0, -3*pi/4, 0, pi/2, pi/4];
T_desired = FK_space(M, S_space, theta_test);

% guess
theta0 = zeros(1, 7);

% redundancy (this is the functiontst)
[theta_sol, success] = redundancy_resolution(S_space, M, T_desired, theta0);

if success
    disp('Target reached while maximizing manipulability!');
    fprintf('Final Joint Angles: '); disp(theta_sol); % for the report too?
    
    % Compare Manipulability (might be nice to see for presentation /
    % report? ) 
    w_start = J_ellipsoid_volume(J_space(S_space, theta0));
    w_final = J_ellipsoid_volume(J_space(S_space, theta_sol));
    fprintf('Manipulability improved from %.4f to %.4f\n', w_start, w_final);
end