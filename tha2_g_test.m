addpath("Functions");

%% Test Script (most classic example i could find) 
clear; clc; close all;
%  3-Link 3 revolute joint easy stuff W6-L1 slide 2 - 7 (its the fwd
%  kinemetics of a 3R 
L1 = 1; L2 = 1; L3 = 1;

% Configuration (Change these to see ellipsoid shapes change!)
theta = [0, pi/4, pi/4]; 

% home position: end-effector at (L1+L2+L3, 0, 0)
M = [1 0 0 L1+L2+L3; 
     0 1 0 0; 
     0 0 1 0; 
     0 0 0 1];

% spatial screw axes (Z-axis)
S1 = [0; 0; 1; 0; 0; 0];          % Joint 1 at origin
S2 = [0; 0; 1; 0; -L1; 0];       % Joint 2 at (L1, 0, 0)
S3 = [0; 0; 1; 0; -(L1+L2); 0];  % Joint 3 at (L1+L2, 0, 0)
Slist = [S1, S2, S3];

% current Pose and Jacobian
T = FK_space(M, Slist, theta);
J = Jacobian_space(Slist, theta); % replace this jsut jacobian space function sarah makes later rn it random matlab forum one

% split jacob
Jw = J(1:3, :); % ang part
Jv = J(4:6, :); % lin part

% testing Metrics
disp('Current Theta:');
disp(theta);
fprintf('Manipulability Metrics \n');

%  it is a planar robot look at the 2D plane (X-Y) for non-zero metrics
%  (was getting 0.00 and 0.000 :P and inf 
Jv_2D = J(4:5, :); 
fprintf('Isotropy (Linear 2D):    %.4f\n', J_isotropy(Jv_2D));
fprintf('Condition Number (2D):   %.4f\n', J_condition(Jv_2D));
fprintf('Ellipsoid Volume (Area): %.4f\n', J_ellipsoid_volume(Jv_2D));

%% plotting 
% fixing the thing i did with points
n = length(theta);
points = zeros(3, n + 1); 

% Define where the joints are at home position
q_home = [0, 0, 0;      % Joint 1
          L1, 0, 0;     % Joint 2
          L1+L2, 0, 0]'; % Joint 3

points(:, 1) = [0; 0; 0]; % Base is fixed
T_curr = eye(4);

for i = 1:n
    % Move to the current joint frame
    T_curr = T_curr * screw_to_exp(Slist(:,i), theta(i));
    
    if i < n
        % Find where the next joint moved to
        p_next = T_curr * [q_home(:, i+1); 1];
        points(:, i+1) = p_next(1:3);
    end
end

% The last point is the end-effector (FK)
T_final = FK_space(M, Slist, theta);
points(:, n+1) = T_final(1:3, 4);

% Combine into the points matrix
figure('Color', 'w', 'Name', 'Robot Manipulability Analysis');
hold on; grid on; axis equal;
view(3);

%links
h_robot = plot3(points(1,:), points(2,:), points(3,:), 'k-o', 'LineWidth', 3, 'MarkerFaceColor', 'k');

% ellipsoids 
ellipsoid_plot_linear(J, T);  % Should be Magenta: [1 0 1]
ellipsoid_plot_angular(J, T); % Should be Green: [0 1 0]

% no legend from listing every individual axis line
h_lin = findobj(gca, 'Type', 'Surface', 'FaceColor', [1 0 1]);
h_ang = findobj(gca, 'Type', 'Surface', 'FaceColor', [0 1 0]);

xlabel('X-axis (m)'); ylabel('Y-axis (m)'); zlabel('Z-axis (m)');
title('Robot Manipulability Analysis');

% Only include the main structure and the two ellipsoid surfaces
if ~isempty(h_lin) && ~isempty(h_ang)
    legend([h_robot, h_lin(1), h_ang(1)], ...
           {'Robot Structure', 'Linear Manipulability', 'Angular Manipulability'}, ...
           'Location', 'northeast');
end

%% Jacobian_space (taken from a matlab forum just for this example lol idk if this right 
% (Calculates the space Jacobian for PoE) W6-L1 slide 7
function J = Jacobian_space(Slist, theta)
    n = length(theta);
    J = zeros(6, n);
    T = eye(4);
    J(:, 1) = Slist(:, 1);
    for i = 2:n
        % This is the correct Adjoint implementation for the Space Jacobian
        T = T * screw_to_exp(Slist(:, i-1), theta(i-1));
        % Adjoint transformation to shift screw to current frame
        R = T(1:3, 1:3);
        p = T(1:3, 4);
        p_sk = [0 -p(3) p(2); p(3) 0 -p(1); -p(2) p(1) 0];
        AdT = [R, zeros(3); p_sk*R, R];
        J(:, i) = AdT * Slist(:, i);
    end
end