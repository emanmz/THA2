% test script for pt. h & i
%% Test Script for 2-Link Planar Arm
clear; clc; close all;
addpath("Functions"); 

L1 = 1; L2 = 1;
% Home position: Arm extended along the x-axis
M = [1 0 0 L1+L2; 
     0 1 0 0; 
     0 0 1 0; 
     0 0 0 1];

S1 = [0; 0; 1; 0; 0; 0];      % Joint 1 at origin
S2 = [0; 0; 1; 0; -L1; 0];   % Joint 2 at (L1, 0, 0)
Slist = [S1, S2];

% Desired Tas
Tsd = [1 0 0 1.0; 
       0 1 0 1.0; 
       0 0 1 0.0; 
       0 0 0 1];

% starting point
theta0 = [0.1, 0.1]; 

% Pseudoinverse
fprintf('J_inverse_kinematics\n');
tic;
[theta_nr, success_nr] = J_inverse_kinematics(Slist, M, Tsd, theta0, 1e-3, 1e-3);
toc;

if success_nr
    fprintf('NR Solved Theta: [%.4f, %.4f]\n', theta_nr(1), theta_nr(2));
    T_final_nr = FK_space(M, Slist, theta_nr);
    fprintf('NR Final Position: [%.2f, %.2f]\n\n', T_final_nr(1,4), T_final_nr(2,4));
else
    fprintf('NR Failed to converge.\n\n');
end

% Jacobian Transpose IK
fprintf('J_transpose_kinematics \n');
tic;
[theta_tr, success_tr] = J_transpose_kinematics(Slist, M, Tsd, theta0);
toc;
% 
if success_tr
    fprintf('Transpose Solved Theta: [%.4f, %.4f]\n', theta_tr(1), theta_tr(2));
    T_final_tr = FK_space(M, Slist, theta_tr);
    fprintf('Transpose Final Position: [%.2f, %.2f]\n', T_final_tr(1,4), T_final_tr(2,4));
else
    fprintf('Transpose Failed (increase max_iter or alpha).\n');
end
