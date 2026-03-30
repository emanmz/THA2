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
%% Configuration Setup
pose_home = zeros(1, 7);
pose_packing = [0, -32.08, 0, -170.17, 0, 0, 45] * (pi/180);
pose_ready = [0, -pi/4, 0, -3*pi/4, 0, pi/2, pi/4];
pose_singularity = zeros(1, 7);

steps = 150;
pause_time = 0.05;

eomg = 1e-3;
ev = 1e-3;

%% Compute SE(3) Targets from FK
T_home = FK_space_no_plot(M, S_space, pose_home);
T_packing = FK_space_no_plot(M, S_space, pose_packing);
T_ready = FK_space_no_plot(M, S_space, pose_ready);
T_singularity = FK_space_no_plot(M, S_space, pose_singularity);

%% 1. Home to Packing (IK-driven)
% fprintf('Animating IK: Home to Packing...\n');
% animate_transition_IK(T_home, T_packing, pose_home, S_space, M, steps, pause_time, eomg, ev);
% 2. Packing to Ready
fprintf('Animating IK: Packing to Ready...\n');
animate_transition_IK(T_packing, T_ready, pose_packing, S_space, M, steps, pause_time);

% % 3. Ready to Singularity
% fprintf('Animating IK: Ready to Singularity...\n');
% animate_transition_IK(T_ready, T_singularity, pose_ready, S_space, M, steps, pause_time, eomg, ev);


%% --- 1. Robot Configuration (Definitions) ---
% Defining these so the code is self-contained and runnable
L = [0.333, 0.316, 0.384, 0.107];
A = 0.088;

% Space Jacobian Screws (S_space) - 7 DOF
S_space = [0, 0, 1, 0, 0, 0;
           0, 1, 0, -L(1), 0, 0;
           0, 0, 1, 0, 0, 0;
           0, 1, 0, -(L(1)+L(2)), 0, A;
           0, 0, 1, 0, 0, 0;
           0, 1, 0, -(L(1)+L(2)+L(3)), 0, 0;
           0, 0, 1, 0, 0, 0]';

% Home Position M
M = [1, 0, 0, A; 
     0, 1, 0, 0; 
     0, 0, 1, L(1)+L(2)+L(3)-L(4); 
     0, 0, 0, 1];

%% --- 2. Configuration Setup ---
pose_home = zeros(1, 7);
pose_packing = [0, -32.08, 0, -170.17, 0, 0, 45] * (pi/180);
pose_ready = [0, -pi/4, 0, -3*pi/4, 0, pi/2, pi/4];

steps = 100;
pause_time = 0.02;

%% --- 3. Compute SE(3) Targets ---
T_packing = FK_space_no_plot(M, S_space, pose_packing);
T_ready = FK_space_no_plot(M, S_space, pose_ready);

%% --- 4. Run Animation ---
fprintf('Animating IK: Packing to Ready using Redundancy Resolution...\n');
animate_transition_IK(T_packing, T_ready, pose_packing, S_space, M, steps, pause_time);

%% --- IK Animation Function ---
function animate_transition_IK(T_start, T_end, theta_init, Slist, M, steps, pause_time)
    % SE(3) Trajectory Interpolation
    T_traj = cell(steps,1);
    for i = 1:steps
        s = (i-1)/(steps-1);
        % Interpolate from T_start to T_end
        T_traj{i} = T_start * expm(VecTose3(MatrixLog6(T_start \ T_end)) * s);
    end

    L = [0.333 0.316 0.384 0.107]; A = 0.088;
    qs_fixed = [0, 0, 0; 0, 0, L(1); 0, 0, L(1); A, 0, L(1)+L(2); 
                0, 0, L(1)+L(2)+L(3); 0, 0, L(1)+L(2)+L(3); 
                A, 0, L(1)+L(2)+L(3)-L(4)]';

    theta_curr = theta_init;
    figure(1); set(gcf, 'Color', 'w');

    for t = 1:steps
        Tsd = T_traj{t};
        
        % Use your redundancy resolution logic
        [theta_sol, success] = redundancy_resolution(Slist, M, Tsd, theta_curr);
        
        if success
            theta_curr = theta_sol;
        end

        clf; hold on; grid on; axis equal; view(45, 30);
        axis([-0.8 0.8 -0.8 0.8 0 1.2]);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        
        draw_arm_fixed(Slist, M, theta_curr, qs_fixed);
        
        % Plot blue target dot
        plot3(Tsd(1,4), Tsd(2,4), Tsd(3,4), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
        
        title(sprintf('Step %d/%d | Tracking Target', t, steps));
        drawnow;
        pause(pause_time);
    end
end

%% --- Redundancy Resolution (Modified for stability) ---
function [theta, success] = redundancy_resolution(Slist, M, Tsd, theta0)
    max_iter = 50;
    k = 0.1; 
    theta = theta0(:)'; 
    success = false;
    
    for i = 1:max_iter
        Tsb = FK_space_no_plot(M, Slist, theta);
        
        % Calculate Space Twist Error
        Vs = MatrixLog6(Tsd / Tsb);
        
        if (norm(Vs(1:3)) < 1e-3) && (norm(Vs(4:6)) < 1e-3)
            success = true;
            break;
        end
        
        Js = J_space(Slist, theta);
        J_pinv = pinv(Js);
        
        % Secondary Task: Gradient of Manipulability
        eps = 1e-4;
        grad_w = zeros(length(theta), 1);
        for j = 1:length(theta)
            th_plus = theta; th_plus(j) = th_plus(j) + eps;
            w_plus = sqrt(det(J_space(Slist, th_plus) * J_space(Slist, th_plus)'));
            
            th_minus = theta; th_minus(j) = th_minus(j) - eps;
            w_minus = sqrt(det(J_space(Slist, th_minus) * J_space(Slist, th_minus)'));
            
            grad_w(j) = (w_plus - w_minus) / (2 * eps);
        end
        
        % Null Space Projection: Move toward goal + maximize manipulability
        null_projection = (eye(size(Js,2)) - J_pinv * Js) * (k * grad_w);
        theta = theta + (J_pinv * Vs + null_projection)';
    end
end

%% --- Core Helper Functions ---
function T = FK_space_no_plot(M, Slist, theta)
    T = eye(4);
    for i = 1:length(theta)
        T = T * expm(VecTose3(Slist(:,i)) * theta(i));
    end
    T = T * M;
end

function Js = J_space(Slist, theta)
    Js = Slist;
    T = eye(4);
    for i = 2:length(theta)
        T = T * expm(VecTose3(Slist(:,i-1)) * theta(i-1));
        Js(:,i) = Adjoint(T) * Slist(:,i);
    end
end

function Ad = Adjoint(T)
    R = T(1:3,1:3); p = T(1:3,4);
    p_skew = [0 -p(3) p(2); p(3) 0 -p(1); -p(2) p(1) 0];
    Ad = [R, zeros(3); p_skew*R, R];
end

function V = MatrixLog6(T)
    R = T(1:3,1:3); p = T(1:3,4);
    if abs(trace(R)-3) < 1e-6
        V = [0;0;0;p];
    else
        phi = acos(max(min((trace(R)-1)/2,1),-1));
        w_sk = (R-R')/(2*sin(phi));
        w = [w_sk(3,2); w_sk(1,3); w_sk(2,1)];
        G_inv = eye(3)/phi - 0.5*w_sk + (1/phi - 0.5*cot(phi/2))*(w_sk^2);
        V = [w; G_inv*p] * phi;
    end
end

function se3mat = VecTose3(V)
    w = V(1:3); v = V(4:6);
    se3mat = [0 -w(3) w(2) v(1); w(3) 0 -w(1) v(2); -w(2) w(1) 0 v(3); 0 0 0 0];
end

function draw_arm_fixed(Slist, M, theta, qs_fixed)
    n = length(theta);
    joint_coords = zeros(3, n + 1);
    joint_coords(:, 1) = [0; 0; 0];
    T_accum = eye(4);
    for i = 1:n
        T_accum = T_accum * expm(VecTose3(Slist(:,i)) * theta(i));
        if i < n
            pos = T_accum * [qs_fixed(:, i+1); 1];
            joint_coords(:, i+1) = pos(1:3);
        end
    end
    T_ee = T_accum * M;
    joint_coords(:, end) = T_ee(1:3, 4);
    plot3(joint_coords(1,:), joint_coords(2,:), joint_coords(3,:), '-ok', 'LineWidth', 3, 'MarkerFaceColor', 'y');
end