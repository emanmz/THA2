%% Test Function for Space Jacobian & Singularity
addpath("Functions");
clear; clc; close all;

%% SCARA Parameters (W7-L1 Slide 9)
L1 = 2; L2 = 2;
theta = [pi/2, pi/2, pi/2, 0.5]; % th1, th2, th3 (revolute), d4 (prismatic)

% Space Jacobian Verification ---
% Define screws in Home Position {s}
% Joints 1, 2, 3 are revolute Z-axes. Joint 4 is prismatic Z-axis.
ws = {[0;0;1], [0;0;1], [0;0;1], [0;0;0]}; 
qs = {[0;0;0], [L1;0;0], [L1+L2;0;0], [0;0;1]}; % v4 for prismatic is the direction [0;0;1]

S_space = zeros(6,4);
for i=1:4
    if norm(ws{i}) == 0 % Prismatic
        S_space(:, i) = [0; 0; 0; qs{i}]; 
    else % Revolute
        S_space(:, i) = [ws{i}; cross(-ws{i}, qs{i})];
    end
end

% Manual Js from Slide 9 (Calculated at current theta)
Js_manual = [ 0, 0, 0, 0;
              0, 0, 0, 0;
              1, 1, 1, 0;
              0, L1*sin(theta(1)), L1*sin(theta(1))+L2*sin(theta(1)+theta(2)), 0;
              0, -L1*cos(theta(1)), -L1*cos(theta(1))-L2*cos(theta(1)+theta(2)), 0;
              0, 0, 0, 1];

% Function Js
Js_func = J_space(S_space, theta);

disp("--- Jacobian Comparison ---")
disp("Is Js_func equal to Js_manual?")
disp(all(abs(Js_func - Js_manual) < 1e-6, 'all')); % Should return 1

% Body Jacobian & Transformation ---
% Home Position M
M = [1 0 0 L1+L2; 0 1 0 0; 0 0 1 0; 0 0 0 1];

% Define Body Screws (Screws relative to {b} at home)
wb = {[0;0;1], [0;0;1], [0;0;1], [0;0;0]};
qb_pos = {[-L1-L2;0;0], [-L2;0;0], [0;0;0], [0;0;1]}; 

B_list = zeros(6,4);
for i=1:4
    if norm(wb{i}) == 0
        B_list(:, i) = [0; 0; 0; qb_pos{i}];
    else
        B_list(:, i) = [wb{i}; cross(-wb{i}, qb_pos{i})];
    end
end

Jb_func = J_body(B_list, theta);
T_sb = FK_space_no_plot(M, S_space, theta);

% THE CHECK: Js = Ad(Tsb) * Jb ONLY if Jb is defined in the body frame
Js_from_Jb = Adjoint(T_sb) * Jb_func;

disp("--- Body to Space Adjoint Check ---")
disp("Does Ad(Tsb)*Jb match J_space?")
disp(all(abs(Js_from_Jb - Js_func) < 1e-6, 'all'));

% Singularity Test 
fprintf('\n--- Singularity Analysis (Analytical) ---\n');
syms q1 q2 q3 q4 real
thetaSym = [q1 q2 q3 q4];
Sings = singularity(S_space, thetaSym);
disp("Singular configurations found for:");
disp(Sings);