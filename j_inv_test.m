%% IK Test Case: Franka Research 3
addpath("Functions");

% Define FR3 Home M and Slist (Simplified placeholder dimensions)
L1 = 0.333; L2 = 0.316; % m
M = [1 0 0 0; 0 1 0 0; 0 0 1 L1+L2; 0 0 0 1];
Slist = rand(6, 7); % Replace with your actual FR3 Screw Axes

% Target Configuration (The "Goal")
theta_dest = [0, -pi/4, 0, -3*pi/4, 0, pi/2, pi/4];
T_desired = FK_space(M, Slist, theta_dest);

% Initial Guess (The "Start")
% IK is sensitive to this; if it's too far, it might hit a singularity
theta_start = theta_dest + 0.1*randn(1, 7); 

% Run IK
[theta_solved, success] = J_inverse_kinematics(Slist, M, T_desired, theta_start);

if success
    disp('IK Converged! Solved Theta:');
    disp(theta_solved);
    
    % Verify by checking FK of the result
    T_check = FK_space(M, Slist, theta_solved);
    fprintf('Final Positional Error: %e meters\n', norm(T_check(1:3,4) - T_desired(1:3,4)));
end