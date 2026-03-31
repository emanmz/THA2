close all; clc; clear all;

%% Link lengths (m)
L = [0.333 0.316 0.384 0.107];
%% Flange offset (m)
A = 0.088;

%% Screw axis definitions
ws = {[0;0;1], [0;-1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qs = {[0;0;0], [0;0;L(1)], [0;0;L(1)], [A;0;L(1)+L(2)], [0;0;L(1)+L(2)+L(3)], ...
      [0;0;L(1)+L(2)+L(3)], [A;0;L(1)+L(2)+L(3)-L(4)]};

%% Home configuration
M = [1 0 0 A; 0 -1 0 0; 0 0 -1 L(1)+L(2)+L(3)-L(4); 0 0 0 1];

S_space = zeros(6,7);
for i = 1:7
    wi = ws{i};
    vi = cross(-wi, qs{i});
    S_space(:, i) = [wi; vi];
end

%% Configuration Setup
pose_home    = zeros(1, 7);
% Standard "Ready" pose for this kinematic chain
pose_ready   = [0, -0.4, 0, -2.0, 0, 1.5, 0.7]; 
% A "Packing" pose that is physically tucked
pose_packing = [0, 0.5, 0, -2.5, 0, 1.0, 0];

% Use 'pose_ready' as your starting guess for the first animation
theta_start_guess = pose_ready;
steps = 300;
pause_time = 0.01;
eomg = 1e-3; ev = 1e-3;

T_home    = FK_space_no_plot(M, S_space, pose_home);
T_packing = FK_space_no_plot(M, S_space, pose_packing);
T_ready   = FK_space_no_plot(M, S_space, pose_ready);

% Define the IMPOSSIBLE target
T_extreme = T_home;
T_extreme(3,4) = T_extreme(3,4) + 0.8; % Push way past reach

%% Run Animations
fprintf('\n=== Running Success Trajectory (Ready -> Home) ===\n');
% Notice: T_ready matches pose_ready
animate_transition_IK_cartesian(T_ready, T_home, pose_ready, S_space, M, qs, steps, pause_time, eomg, ev);

fprintf('\n=== Running Failure Trajectory (Home -> Impossible) ===\n');
% Notice: T_home matches pose_home
% animate_transition_IK_cartesian(T_home, T_ready, pose_home, S_space, M, qs, steps, pause_time, eomg, ev);
%% --- Core Functions ---

function animate_transition_IK_cartesian(T_start, T_end, theta_init, Slist, M, qs, steps, pause_time, eomg, ev)
    % 1. Setup GIF filename
    gif_filename = 'robot_animation.gif';
    
    T_rel = inv(T_start) * T_end;
    V_rel = MatrixLog6(T_rel);
    se3_rel = VecTose3(V_rel);
    theta_guess = theta_init;
    
    h_fig = figure('Color','w'); % Get handle to the figure
    
    for t = 1:steps
        s = (t-1)/(steps-1);
        Tsd = T_start * MatrixExp6(se3_rel, s);
        
        [theta_sol, success, info] = J_inverse_kinematics_diagnostics(Slist, M, Tsd, theta_guess, eomg, ev);
        theta_guess = theta_sol;
        
        T_curr = FK_space_no_plot(M, Slist, theta_guess);
        pos_err = norm(T_curr(1:3,4) - Tsd(1:3,4));
        status_str = ternary(success, 'OK', 'NO CONV');
        
        clf; hold on; grid on; axis equal; view(3);
        axis([-0.8 0.8 -0.8 0.8 -0.1 1.2]);
        draw_arm_fixed(Slist, M, theta_guess, qs);
        plot3(Tsd(1,4), Tsd(2,4), Tsd(3,4), 'r*', 'MarkerSize', 10);
        
        title(sprintf('Step %d/%d | w err=%.1e | v err=%.1e | %s', t, steps, info.w_err, info.v_err, status_str));
        drawnow;

        % --- GIF CAPTURE LOGIC ---
        frame = getframe(h_fig);      % Capture the figure window
        im = frame2im(frame);         % Convert to image
        [imind, cm] = rgb2ind(im, 256); % Convert to indexed image (required for GIF)
        
        if t == 1
            % Create the file on the first step
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', pause_time);
        else
            % Append to the file for all other steps
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', pause_time);
        end
        % -------------------------
        
        pause(pause_time);
    end
    fprintf('GIF saved as: %s\n', gif_filename);
end

function [theta, success, info] = J_inverse_kinematics_diagnostics(Slist, M, Tsd, theta0, eomg, ev)
    % --- Your provided IK Logic (pt. h) ---
    max_iter = 50; 
    theta = theta0;
    i = 0;
    
    Tsb = FK_space_no_plot(M, Slist, theta);
    [Vb, ~] = calculate_twist_error(Tsb, Tsd);
    
    % Current error magnitudes
    w_err = norm(Vb(1:3));
    v_err = norm(Vb(4:6));
    
    % Continue until error is below tolerances or max_iter is hit
    while (w_err > eomg || v_err > ev) && i < max_iter
        % Space Jacobian at current theta
        Js = J_space(Slist, theta);
        
        % Body Twist (Vb) to Space Twist (Vs) 
        Vs = Adjoint(Tsb) * Vb; 
        
        % Newton-Raphson update
        theta = theta + (pinv(Js) * Vs)';
        
        % Update current state
        i = i + 1;
        Tsb = FK_space_no_plot(M, Slist, theta);
        [Vb, ~] = calculate_twist_error(Tsb, Tsd);
        w_err = norm(Vb(1:3));
        v_err = norm(Vb(4:6));
    end
    
    % --- Diagnostics for the GIF/Title ---
    success = (w_err <= eomg && v_err <= ev);
    Js_final = J_space(Slist, theta);
    S = svd(Js_final);
    
    info.w_err = w_err;
    info.v_err = v_err;
    info.iters = i;
    info.condJ = cond(Js_final);
    info.sigma_min = min(S);
    info.reason = ternary(success, 'None', 'Max Iterations / Singular');
end
%% --- Geometry & Math Helpers ---

function draw_arm_fixed(Slist, M_ee, theta, qs)
    n = length(theta);
    joint_coords = zeros(3, n+1);
    joint_coords(:,1) = [0;0;0];
    T_prev = eye(4);
    for i = 2:n
        T_prev = T_prev * MatrixExp6(VecTose3(Slist(:,i-1)), theta(i-1));
        p_world = T_prev * [qs{i}; 1];
        joint_coords(:,i) = p_world(1:3);
    end
    T_ee = FK_space_no_plot(M_ee, Slist, theta);
    joint_coords(:,n+1) = T_ee(1:3,4);
    plot3(joint_coords(1,:), joint_coords(2,:), joint_coords(3,:), '-o', 'LineWidth', 3, 'Color', [0.2 0.2 0.2]);
end

function T = FK_space_no_plot(M, Slist, theta)
    T = eye(4);
    for i = 1:length(theta)
        T = T * MatrixExp6(VecTose3(Slist(:,i)), theta(i));
    end
    T = T * M;
end

function Js = J_space(Slist, theta)
    n = size(Slist,2); Js = zeros(6,n); Js(:,1) = Slist(:,1); T = eye(4);
    for i = 2:n
        T = T * MatrixExp6(VecTose3(Slist(:,i-1)), theta(i-1));
        Js(:,i) = Adjoint(T) * Slist(:,i);
    end
end

function AdT = Adjoint(T)
    R = T(1:3,1:3); p_hat = VecToso3(T(1:3,4));
    AdT = [R, zeros(3); p_hat*R, R];
end

function [Vb, err] = calculate_twist_error(Tsb, Tsd)
    T_err = inv(Tsb) * Tsd; Vb = MatrixLog6(T_err); err = norm(Vb);
end

function V = MatrixLog6(T)
    R = T(1:3, 1:3); p = T(1:3, 4);
    if norm(R - eye(3), 'fro') < 1e-10, V = [0;0;0; p]; return; end
    theta = acos(max(min((trace(R) - 1) / 2, 1), -1));
    so3mat = theta/(2*sin(theta)) * (R - R');
    w_hat = so3mat / theta;
    G_inv = eye(3)/theta - 0.5*w_hat + (1/theta - 0.5*cot(theta/2))*(w_hat*w_hat);
    V = [so3ToVec(so3mat); G_inv * p * theta];
end

function T = MatrixExp6(se3mat, theta)
    omg_skew = se3mat(1:3,1:3); v = se3mat(1:3,4);
    omg_vec = so3ToVec(omg_skew);
    ang = norm(omg_vec);
    if ang < 1e-12
        T = [eye(3), v * theta; 0 0 0 1];
    else
        phi = ang * theta;
        w_hat = omg_skew / ang;
        R = eye(3) + sin(phi)*w_hat + (1-cos(phi))*(w_hat*w_hat);
        G = eye(3)*theta + ((1-cos(phi))/ang)*w_hat + ((phi-sin(phi))/(ang^2))*(w_hat*w_hat);
        T = [R, G*v; 0 0 0 1];
    end
end

function se3mat = VecTose3(V)
    se3mat = [VecToso3(V(1:3)), V(4:6); 0 0 0 0];
end

function so3mat = VecToso3(w)
    so3mat = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

function w = so3ToVec(so3mat)
    w = [so3mat(3,2); so3mat(1,3); so3mat(2,1)];
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end