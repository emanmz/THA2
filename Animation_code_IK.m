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

eomg = 1e-2;
ev = 1e-2;

%% Compute SE(3) Targets from FK
T_home = FK_space_no_plot(M, S_space, pose_home);
T_packing = FK_space_no_plot(M, S_space, pose_packing);
T_ready = FK_space_no_plot(M, S_space, pose_ready);
T_singularity = FK_space_no_plot(M, S_space, pose_singularity);

%% 1. Home to Packing (IK-driven)
fprintf('Animating IK: Home to Packing...\n');
animate_transition_IK(T_home, T_packing, pose_home, S_space, M, steps, pause_time, eomg, ev);

%% 2. Packing to Ready
fprintf('Animating IK: Packing to Ready...\n');
animate_transition_IK(T_packing, T_ready, pose_packing, S_space, M, steps, pause_time, eomg, ev);

%% 3. Ready to Singularity
fprintf('Animating IK: Ready to Singularity...\n');
animate_transition_IK(T_ready, T_singularity, pose_ready, S_space, M, steps, pause_time, eomg, ev);

%% --- IK Animation Function ---

function animate_transition_IK(T_start, T_end, theta_init, Slist, M, steps, pause_time, eomg, ev)

% Precompute trajectory in SE(3)
T_traj = cell(steps,1);
for i = 1:steps
    s = (i-1)/(steps-1);
    V = MatrixLog6(inv(T_start) * T_end);   % 6x1 twist
    se3mat = VecTose3(V);                  % convert to 4x4
    T_traj{i} = MatrixExp6(se3mat, s) * T_start;
end

% Robot geometry
L = [0.333 0.316 0.384 0.107];
A = 0.088;
qs_fixed = [0, 0, 0;
    0, 0, L(1);
    0, 0, L(1);
    A, 0, L(1)+L(2);
    0, 0, L(1)+L(2)+L(3);
    0, 0, L(1)+L(2)+L(3);
    A, 0, L(1)+L(2)+L(3)-L(4)]';

theta_guess = theta_init;

figure(1); set(gcf, 'Color', 'w');

for t = 1:steps
    clf; hold on; grid on; axis equal;
    view(3); axis([-0.8 0.8 -0.8 0.8 0 1.2]);
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % Solve IK for current desired pose
    Tsd = T_traj{t};
    [theta_sol, success, w_err, v_err] = J_inverse_kinematics(Slist, M, Tsd, theta_guess, eomg, ev);
    T_current = FK_space_no_plot(M, Slist, theta_guess);
    p_actual = T_current(1:3, 4);
    p_desired = Tsd(1:3, 4);
    pos_error = norm(p_actual - p_desired);
    % btw if we dont want to make the robot shake her ass like crazy we can
    % change it this
    % if success
    theta_guess = theta_sol;
    % end

    draw_arm_fixed(Slist, M, theta_guess, qs_fixed);


    title(sprintf('Step %d / %d | ang Error = %.4f m | lin Error = %.4f m', t, steps, w_err, v_err));
    drawnow;
    pause(pause_time);
end
end

%% drawwwwwing func

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

%% funcs 
function [theta, success, w_err, v_err] = J_inverse_kinematics(Slist, M, Tsd, theta0, eomg, ev)
    max_iter = 150;
    theta = theta0(:);   % force column vector for consistency
    i = 0;
    
    Tsb = FK_space_no_plot(M, Slist, theta');
    [Vb, ~] = calculate_twist_error(Tsb, Tsd);

    % Split twist error into angular + linear parts
    w_err = norm(Vb(1:3));
    v_err = norm(Vb(4:6));
    
    while ((w_err > eomg) || (v_err > ev)) && i < max_iter
        Js = J_space(Slist, theta');
        Vs = Adjoint(Tsb) * Vb; 
        
        % Pseudoinverse IK update
        theta = theta + pinv(Js) * Vs;
        
        i = i + 1;
        
        % Recompute FK + twist error
        Tsb = FK_space_no_plot(M, Slist, theta');
        [Vb, ~] = calculate_twist_error(Tsb, Tsd);

        % Update convergence checks
        w_err = norm(Vb(1:3));
        v_err = norm(Vb(4:6));
    end
    
    success = (w_err <= eomg) && (v_err <= ev);
    theta = theta';   % return row vector to match your other code

    if ~success
        warning('IK did not converge within %d iterations | w_err = %.4e, v_err = %.4e', ...
                max_iter, w_err, v_err);
    end
end
function V = MatrixLog6(T)
    % Extracts the twist V = [w; v] from a transformation matrix T
    R = T(1:3, 1:3);
    p = T(1:3, 4);
    
    if abs(trace(R) - 3) < 1e-6
        % Pure translation / no rotation
        w = [0; 0; 0];
        v = p;
        V = [w; v];
    else
        % Clamp acos argument for numerical stability
        c = (trace(R) - 1) / 2;
        c = max(min(c, 1), -1);

        theta = acos(c);
        w_sk = (1 / (2 * sin(theta))) * (R - R');
        w = [w_sk(3,2); w_sk(1,3); w_sk(2,1)];
        
        % Translation part of twist
        G_inv = eye(3)/theta - 0.5*w_sk + ...
                (1/theta - 0.5*cot(theta/2))*(w_sk^2);
        v = G_inv * p;

        V = [w; v] * theta;
    end
end

function T = MatrixExp6(Smat, theta)
    % NOTE: kept in case you use it elsewhere
    omega_skew = Smat(1:3, 1:3);
    v = Smat(1:3, 4);
    omg = [omega_skew(3,2); omega_skew(1,3); omega_skew(2,1)];
    
    if norm(omg) < 1e-6
        T = [eye(3), v * theta; 0 0 0 1];
    else
        R = axisangle2rot(omg, theta);
        G = eye(3)*theta + (1-cos(theta))*omega_skew + ...
            (theta-sin(theta))*(omega_skew^2);
        T = [R, G*v; 0 0 0 1];
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

function R = axisangle2rot(axis, angle)
    % Normalize axis so it's a unit vector
    if norm(axis) > 0
        axis = axis / norm(axis);
    end 

    ux = axis(1);
    uy = axis(2);
    uz = axis(3);

    % Skew symmetric matrix
    K = [0, -uz, uy;
         uz, 0, -ux;
        -uy, ux, 0];

    % Rodrigues formula
    R = cos(angle) * eye(3) + ...
        (1 - cos(angle)) * (axis * axis') + ...
        sin(angle) * K;
end