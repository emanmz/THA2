addpath("Functions");
close all; clc;
%% Franka Emika

% link lengths, m
L = [0.333 0.316 0.384 0.107];

% flange offset, m
A = 0.088;


% home position, https://frankarobotics.github.io/docs/robot_specifications.html#kinematic-configuration

M = [1 0 0 A;
    0 -1 0 0;
    0 0 -1 L(1)+L(2)+L(3)-L(4);
    0 0 0 1];
S_space = zeros(6,7);

for i=1:7
    wi = ws{i};
    if norm(wi) == 0
        vi = qs{i};
    else
    vi = cross(-wi, qs{i});
    end
    
    S_space(:, i) = [wi; vi];
end
disp("S_space")
disp(S_space)
%% --- Configuration Setup ---
pose_home = zeros(1, 7);
pose_packing = [0, -32.08, 0, -170.17, 0, 0, 45] * (pi/180);
pose_ready = [0, -pi/4, 0, -3*pi/4, 0, pi/2, pi/4];
steps = 150;
pause_time = 0.02; % Speed up animation slightly

%% --- Compute SE(3) Targets ---
T_packing = FK_space_no_plot(M, S_space, pose_packing);
T_ready = FK_space_no_plot(M, S_space, pose_ready);

%% --- Execute Animation ---
fprintf('Animating IK: Packing to Ready using DLS...\n');
% Note: Calling with 7 arguments to match your script logic
animate_transition_IK(T_packing, T_ready, pose_packing, S_space, M, steps, pause_time);

%% --- Animation Function ---
function animate_transition_IK(T_start, T_end, theta_init, Slist, M, steps, pause_time)
    % Precompute Trajectory
    T_traj = cell(steps,1);
    for i = 1:steps
        s = (i-1)/(steps-1);
        % Standard screw linear interpolation
        T_traj{i} = T_start * expm(VecTose3(MatrixLog6(T_start \ T_end)) * s);
    end

    % Robot Plotting Constants
    L = [0.333 0.316 0.384 0.107]; A = 0.088;
    qs_fixed = [0, 0, 0; 0, 0, L(1); 0, 0, L(1); A, 0, L(1)+L(2); 
                0, 0, L(1)+L(2)+L(3); 0, 0, L(1)+L(2)+L(3); 
                A, 0, L(1)+L(2)+L(3)-L(4)]';
    
    theta_curr = theta_init;
    figure(1); set(gcf, 'Color', 'w');

    for t = 1:steps
        Tsd = T_traj{t};
        
        % CALL DLS INVERSE KINEMATICS HERE
        [theta_sol, success] = DLS_inverse_kinematics(Slist, M, Tsd, theta_curr);
        
        if success
            theta_curr = theta_sol;
        end

        clf; hold on; grid on; axis equal; view(3);
        axis([-0.8 0.8 -0.8 0.8 0 1.2]);
        
        draw_arm_fixed(Slist, M, theta_curr, qs_fixed);
        plot3(Tsd(1,4), Tsd(2,4), Tsd(3,4), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        
        title(sprintf('Step %d / %d | DLS Tracking', t, steps));
        drawnow; pause(pause_time);
    end
end

%% --- Damped Least Squares IK ---
function [theta, success] = DLS_inverse_kinematics(Slist, M, Tsd, theta0)
    max_iter = 100;
    lambda_max = 0.1; 
    epsilon = 0.05;   
    theta = theta0(:)'; % Ensure row vector
    success = false;
    
    for i = 1:max_iter
        Tsb = FK_space_no_plot(M, Slist, theta);
        
        % Calculate Space Twist Error Vs
        V_mat = MatrixLog6(Tsd / Tsb);
        % Check tolerances (eomg=1e-3, ev=1e-3)
        if (norm(V_mat(1:3)) < 1e-3) && (norm(V_mat(4:6)) < 1e-3)
            success = true;
            break;
        end
        
        Js = J_space(Slist, theta);
        Vs = V_mat; 
        
        % Singularity Detection via SVD
        s = svd(Js);
        s_min = min(s);
        
        % Engaged Damping logic
        if s_min < epsilon
            lambda = (1 - (s_min/epsilon)^2) * lambda_max^2;
        else
            lambda = 0;
        end
        
        % DLS Update: Delta_Theta = J' * inv(J*J' + lambda*I) * Vs
        % Using the backslash or explicit inverse for Damped pseudo-inverse
        d_theta = (Js' / (Js * Js' + lambda * eye(6))) * Vs;
        theta = theta + d_theta';
    end
end

function se3mat = VecTose3(V)
    % Converts a 6x1 twist vector into a 4x4 se(3) matrix
    w = V(1:3);
    v = V(4:6);
    
    w_skew = [   0   -w(3)  w(2);
               w(3)   0    -w(1);
              -w(2)  w(1)   0 ];
    
    se3mat = [w_skew, v;
              0 0 0 0];
end

function draw_arm_fixed(Slist, M, theta, qs_fixed)
    n = length(theta);
    joint_coords = zeros(3, n + 1);
    joint_coords(:, 1) = [0; 0; 0];
    
    T_accum = eye(4);
    for i = 1:n
        bracket_S = [0, -Slist(3,i), Slist(2,i), Slist(4,i);
                     Slist(3,i), 0, -Slist(1,i), Slist(5,i);
                    -Slist(2,i), Slist(1,i), 0, Slist(6,i);
                     0, 0, 0, 0];
        T_accum = T_accum * expm(bracket_S * theta(i));
        
        if i < n
            home_pos_next = [qs_fixed(:, i+1); 1];
            current_pos_next = T_accum * home_pos_next;
            joint_coords(:, i+1) = current_pos_next(1:3);
        end
    end
    
    T_ee = T_accum * M;
    joint_coords(:, end) = T_ee(1:3, 4);
    
    plot3(joint_coords(1,:), joint_coords(2,:), joint_coords(3,:), ...
          '-o', 'LineWidth', 4, 'MarkerSize', 6, ...
          'MarkerFaceColor', 'k', 'Color', [0.2 0.2 0.2]);

    plot3(joint_coords(1,end), joint_coords(2,end), joint_coords(3,end), ...
          'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end
