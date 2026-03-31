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


% 2. Packing to Ready
fprintf('Animating IK: Packing to Ready...\n');
animate_transition_IK(T_packing, T_ready, pose_packing, S_space, M, steps, pause_time, eomg, ev);

% 3. Ready to Singularity
fprintf('Animating IK: Ready to Singularity...\n');
animate_transition_IK(T_ready, T_singularity, pose_ready, S_space, M, steps, pause_time, eomg, ev);

%% --- IK Animation Function ---
function animate_transition_IK(T_start, T_end, theta_init, Slist, M, steps, pause_time, eomg, ev)
    % Correct trajectory interpolation
    T_traj = cell(steps,1);
    for i = 1:steps
        s = (i-1)/(steps-1);
        % Matrix power for interpolation: T_start * exp(log(inv(T_start)*T_end)*s)
        T_traj{i} = T_start * MatrixExp6(VecTose3(MatrixLog6(inv(T_start)*T_end)), s);
    end

    % Robot geometry for plotting
    L = [0.333 0.316 0.384 0.107]; A = 0.088;
    qs_fixed = [0, 0, 0; 0, 0, L(1); 0, 0, L(1); A, 0, L(1)+L(2); 
                0, 0, L(1)+L(2)+L(3); 0, 0, L(1)+L(2)+L(3); 
                A, 0, L(1)+L(2)+L(3)-L(4)]';

    theta_guess = theta_init;
    figure(1); set(gcf, 'Color', 'w');

    for t = 1:steps
        Tsd = T_traj{t};
        [theta_sol, success] = J_transpose_kinematics(Slist, M, Tsd, theta_guess, eomg, ev);
        
        if success
            theta_guess = theta_sol;
        end

        clf; hold on; grid on; axis equal; view(3);
        axis([-0.8 0.8 -0.8 0.8 0 1.2]);
        
        draw_arm_fixed(Slist, M, theta_guess, qs_fixed);
        plot3(Tsd(1,4), Tsd(2,4), Tsd(3,4), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        
        title(sprintf('Step %d / %d | IK Success: %s', t, steps, string(success)));
        drawnow; pause(pause_time);
    end
end

%% --- Drawing function ---
function draw_arm_fixed(Slist, M, theta, qs_fixed)
    n = length(theta);
    joint_coords = zeros(3, n + 1);
    joint_coords(:, 1) = [0; 0; 0];
    
    T_accum = eye(4);
    for i = 1:n
        bracket_S = [0, -Slist(3,i),  Slist(2,i), Slist(4,i);
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

%% --- IK + helper funcs ---

function [theta, success] = J_transpose_kinematics(Slist, M, Tsd, theta0, eomg, ev)
   max_iter = 1000; % Increased iterations
    alpha = 0.5;    % Increased step size for Jacobian Transpose
    theta = theta0(:);
    theta = theta0(:); % Ensure theta is a column vector
    success = false;
    
    for i = 1:max_iter
        % 1. Get current EE pose
        Tsb = FK_space_no_plot(M, Slist, theta');
        
        % 2. Calculate the matrix log of the transformation error
        % This returns a 6x1 vector: [omega_x; omega_y; omega_z; vx; vy; vz]
        Vs = MatrixLog6(Tsd / Tsb); 
        
        % 3. Check convergence
        if (norm(Vs(1:3)) < eomg) && (norm(Vs(4:6)) < ev)
            success = true;
            break;
        end
        
        % 4. Compute Space Jacobian
        Js = J_space(Slist, theta');
        
        % 5. Jacobian Transpose Update
        % theta_new = theta_old + alpha * J^T * Twist_Error
        theta = theta + alpha * (Js' * Vs);
    end
    theta = theta'; % Convert back to row vector for the rest of your script
end


%% 
function T = FK_space_no_plot(M, Slist, theta)
    T = eye(4);
    for i = 1:length(theta)
        T = T * MatrixExp6(VecTose3(Slist(:,i)), theta(i));
    end
    T = T * M;
end

function Js = J_space(Slist, theta)
    Js = Slist;
    T = eye(4);
    for i = 2:length(theta)
        T = T * MatrixExp6(VecTose3(Slist(:,i-1)), theta(i-1));
        Js(:,i) = Adjoint(T) * Slist(:,i);
    end
end

function Ad = Adjoint(T)
    R = T(1:3,1:3); p = T(1:3,4);
    p_skew = [0 -p(3) p(2); p(3) 0 -p(1); -p(2) p(1) 0];
    Ad = [R, zeros(3); p_skew*R, R];
end

function T = MatrixExp6(se3mat, theta)
    % Simplified for se3mat already being a 4x4 matrix
    T = expm(se3mat * theta);
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