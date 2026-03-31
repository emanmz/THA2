%% Pt M: Bonus - Robotic System Toolbox Simulation
close all; clear; clc;

% 1. Import the actual Franka URDF
try
    robot = loadrobot('frankaEmikaPanda', 'DataFormat', 'row');
catch
    error('Robotic System Toolbox Panda model not found.');
end

% Detect a reasonable URDF end-effector body name
eeName = get_end_effector_name(robot);
fprintf('Using URDF end-effector body: %s\n', eeName);

%% Franka Emika Math Parameters
L = [0.333 0.316 0.384 0.107];
A = 0.088;
M = [1 0 0 A; 
     0 -1 0 0; 
     0 0 -1 L(1)+L(2)+L(3)-L(4); 
     0 0 0 1];

ws = {[0;0;1], [0;-1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qs = {[0;0;0], ...
      [0;0;L(1)], ...
      [0;0;L(1)], ...
      [A;0;L(1)+L(2)], ...
      [0;0;L(1)+L(2)+L(3)], ...
      [0;0;L(1)+L(2)+L(3)], ...
      [A;0;L(1)+L(2)+L(3)-L(4)]};

S_space = zeros(6,7);
for i = 1:7
    wi = ws{i};
    vi = cross(-wi, qs{i});
    S_space(:, i) = [wi; vi];
end

pose_home = zeros(1, 7);
% pose_packing = [0, -32.08, 0, -170.17, 0, 0, 45] * (pi/180);
pose_packing = [0, -2, 0, -2, 0, 0, 2] * (pi/180);
pose_ready = [0, -pi/4, 0, -3*pi/4, 0, pi/2, pi/4];
pose_singularity = zeros(1, 7);

T_home = FK_space_no_plot(M, S_space, pose_home);
T_packing = FK_space_no_plot(M, S_space, pose_packing);
T_ready = FK_space_no_plot(M, S_space, pose_ready);
T_singularity = FK_space_no_plot(M, S_space, pose_singularity);

eomg = 1e-2;
ev = 1e-2;

% Quick singularity check at home
J0 = J_space(S_space, pose_home);
fprintf('Home pose Jacobian rank = %d, cond(J) = %.4e\n', rank(J0,1e-6), cond(J0));

%% Animation Setup
steps = 300;          % lower for easier debugging
pause_time = 0.02;

% %% 1. Home to Packing
% fprintf('\nAnimating IK: Home to Packing...\n');
% animate_transition_IK(robot, eeName, T_home, T_packing, pose_home, S_space, M, steps, pause_time, eomg, ev);

%% 2. Packing to Ready
fprintf('\nAnimating IK: Packing to Ready...\n');
animate_transition_IK(robot, eeName, T_home, T_ready, pose_home, S_space, M, steps, pause_time, eomg, ev);

%% 3. Ready to Singularity
fprintf('\nAnimating IK: Ready to Singularity...\n');
animate_transition_IK(robot, eeName, T_ready, T_singularity, pose_ready, S_space, M, steps, pause_time, eomg, ev);

%% --- IK Function (ALGORITHM UNCHANGED) ---
function [theta, success, w_err, v_err] = J_inverse_kinematics(Slist, M, Tsd, theta0, eomg, ev)
    max_iter = 200;
    alpha = 0.35;          % step size to avoid overshoot
    lambda = 1e-2;         % damping for near singularities

    theta = theta0(:);     % column vector
    i = 0;

    Tsb = FK_space_no_plot(M, Slist, theta');
    Vb = TwistErrorBody(Tsb, Tsd);

    w_err = norm(Vb(1:3));
    v_err = norm(Vb(4:6));

    while ((w_err > eomg) || (v_err > ev)) && i < max_iter
        Js = J_space(Slist, theta');
        Vs = Adjoint(Tsb) * Vb;   % convert body error to space frame

        % Damped least squares update (UNCHANGED)
        dtheta = Js' / (Js*Js' + (lambda^2)*eye(6)) * Vs;
        theta = theta + alpha * dtheta;

        i = i + 1;

        % Recompute pose and error
        Tsb = FK_space_no_plot(M, Slist, theta');
        Vb = TwistErrorBody(Tsb, Tsd);

        w_err = norm(Vb(1:3));
        v_err = norm(Vb(4:6));
    end

    success = (w_err <= eomg) && (v_err <= ev);
    theta = theta';   % return row

    if ~success
        warning('IK did not converge within %d iterations | w_err = %.4e, v_err = %.4e', ...
                max_iter, w_err, v_err);
    end
end

%% --- Animation Function with model mismatch diagnostics ---
function animate_transition_IK(robot, eeName, T_start, T_end, theta_init, Slist, M, steps, pause_time, eomg, ev)
    % Correct SE(3) interpolation:
    % T(s) = T_start * exp( log(T_start^{-1} T_end) * s )
    T_rel = inv(T_start) * T_end;
    se3Log = MatrixLog6_SE3(T_rel);

    T_traj = cell(steps,1);
    for i = 1:steps
        s = (i-1)/(steps-1);
        T_traj{i} = T_start * MatrixExp6_SE3(se3Log * s);
    end

    theta_guess = theta_init;

    hFig = figure('Color', 'w');
    ax = axes(hFig);

    success_count = 0;
    fail_count = 0;

    fprintf('------------------------------------------------------------\n');
    fprintf('Trajectory start\n');
    fprintf('Steps = %d | eomg = %.3e | ev = %.3e\n', steps, eomg, ev);
    fprintf('Legend: red*=target(custom), green o=achieved(custom), blue s=achieved(URDF)\n');
    fprintf('------------------------------------------------------------\n');

    for t = 1:steps
        Tsd = T_traj{t};

        % --- Solve IK using YOUR custom model ---
        [theta_sol, success, w_err, v_err] = J_inverse_kinematics(Slist, M, Tsd, theta_guess, eomg, ev);
        theta_guess = theta_sol;  % always carry best available solution forward

        if success
            success_count = success_count + 1;
        else
            fail_count = fail_count + 1;
        end

        % --- Custom model achieved pose ---
        T_custom = FK_space_no_plot(M, Slist, theta_guess);
        p_target = Tsd(1:3, 4);
        p_custom = T_custom(1:3, 4);

        % --- URDF robot configuration ---
        full_config = [theta_guess, 0, 0];  % 7 arm + 2 fingers

        % --- URDF end-effector pose (same joint values, different model) ---
        try
            T_urdf = getTransform(robot, full_config, eeName);
            p_urdf = tform2trvec(T_urdf)';
            R_urdf = T_urdf(1:3,1:3);
        catch
            % fallback if getTransform output handling differs
            T_urdf = getTransform(robot, full_config, eeName);
            p_urdf = T_urdf(1:3,4);
            R_urdf = T_urdf(1:3,1:3);
        end

        % --- Diagnostics ---
        Js = J_space(Slist, theta_guess);
        rankJ = rank(Js, 1e-6);
        condJ = cond(Js);

        % Custom model error (this is what YOUR IK is actually solving)
        p_err_custom_vec = p_target - p_custom;
        p_err_custom = norm(p_err_custom_vec);

        R_err_custom = T_custom(1:3,1:3)' * Tsd(1:3,1:3);
        rot_trace_custom = max(min((trace(R_err_custom)-1)/2,1),-1);
        ang_err_custom_deg = rad2deg(acos(rot_trace_custom));

        % URDF-vs-custom mismatch (this shows the two-model problem)
        urdf_custom_mismatch_vec = p_urdf - p_custom;
        urdf_custom_mismatch = norm(urdf_custom_mismatch_vec);

        % URDF-vs-target visual mismatch
        urdf_target_err_vec = p_target - p_urdf;
        urdf_target_err = norm(urdf_target_err_vec);

        % Console diagnostics every 25 steps, on failure, and final step
        if mod(t,25)==1 || ~success || t==steps
            fprintf('\nStep %4d / %4d | %s\n', t, steps, ternary(success,'SUCCESS','FAIL'));
            fprintf('  IK errors (custom model):\n');
            fprintf('    Angular twist err  = %.4e rad\n', w_err);
            fprintf('    Linear twist err   = %.4e m\n', v_err);
            fprintf('    Position err       = %.4e m\n', p_err_custom);
            fprintf('    Rotation err       = %.3f deg\n', ang_err_custom_deg);

            fprintf('  Jacobian:\n');
            fprintf('    rank(J)            = %d\n', rankJ);
            fprintf('    cond(J)            = %.4e\n', condJ);

            fprintf('  Positions:\n');
            fprintf('    Target (custom)    = [%.4f %.4f %.4f]\n', p_target(1), p_target(2), p_target(3));
            fprintf('    Achieved (custom)  = [%.4f %.4f %.4f]\n', p_custom(1), p_custom(2), p_custom(3));
            fprintf('    Achieved (URDF)    = [%.4f %.4f %.4f]\n', p_urdf(1), p_urdf(2), p_urdf(3));

            fprintf('  Model mismatch:\n');
            fprintf('    URDF - Custom EE   = %.4e m\n', urdf_custom_mismatch);
            fprintf('    URDF - Target      = %.4e m\n', urdf_target_err);

            if urdf_custom_mismatch > 2e-2
                fprintf('    >>> WARNING: large URDF/custom mismatch (different kinematic models)\n');
            end

            if rankJ < 6
                fprintf('    >>> WARNING: Jacobian rank deficient (singular / near-singular)\n');
            elseif condJ > 1e4
                fprintf('    >>> WARNING: Jacobian ill-conditioned (near singularity)\n');
            end
        end

        % --- Plot ---
        cla(ax);
        show(robot, full_config, 'Parent', ax, 'Frames', 'off', 'PreservePlot', false);
        hold(ax, 'on');

        % Red star = target from custom model
        plot3(ax, p_target(1), p_target(2), p_target(3), ...
            'r*', 'MarkerSize', 10, 'LineWidth', 2);

        % Green circle = achieved EE from custom model
        plot3(ax, p_custom(1), p_custom(2), p_custom(3), ...
            'go', 'MarkerSize', 8, 'LineWidth', 2);

        % Blue square = achieved EE from URDF robot model
        plot3(ax, p_urdf(1), p_urdf(2), p_urdf(3), ...
            'bs', 'MarkerSize', 8, 'LineWidth', 2);

        % Line from custom achieved to target (actual IK error)
        plot3(ax, [p_custom(1), p_target(1)], ...
                  [p_custom(2), p_target(2)], ...
                  [p_custom(3), p_target(3)], ...
                  'k--', 'LineWidth', 1.5);

        % Line from URDF achieved to target (visual mismatch)
        plot3(ax, [p_urdf(1), p_target(1)], ...
                  [p_urdf(2), p_target(2)], ...
                  [p_urdf(3), p_target(3)], ...
                  'm:', 'LineWidth', 1.5);

        title(ax, sprintf(['Step %d/%d | %s | w=%.2e rad | v=%.2e m | ', ...
                           'rankJ=%d | condJ=%.1e | URDF/custom=%.2e m'], ...
                           t, steps, ternary(success,'OK','NO CONV'), ...
                           w_err, v_err, rankJ, condJ, urdf_custom_mismatch));

        axis(ax, [-0.8 0.8 -0.8 0.8 0 1.2]);
        view(ax, 3);
        grid(ax, 'on');
        legend(ax, {'Target (custom)', 'Achieved (custom)', 'Achieved (URDF)', ...
                    'Custom error', 'URDF visual mismatch'}, ...
                    'Location', 'northeast');

        drawnow limitrate;
        pause(pause_time);
    end

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Trajectory complete\n');
    fprintf('Successful IK solves : %d / %d\n', success_count, steps);
    fprintf('Failed IK solves     : %d / %d\n', fail_count, steps);
    fprintf('------------------------------------------------------------\n');
end

%% --- Body Twist Error ---
function Vb = TwistErrorBody(Tsb, Tsd)
    T_err = inv(Tsb) * Tsd;          % body-frame error transform
    se3mat = MatrixLog6_SE3(T_err);  % 4x4 se(3)
    Vb = se3ToVec(se3mat);           % 6x1 twist
end

%% --- Matrix Log of SE(3): returns 4x4 se(3) matrix ---
function se3mat = MatrixLog6_SE3(T)
    R = T(1:3,1:3);
    p = T(1:3,4);

    if norm(R - eye(3), 'fro') < 1e-8
        % Pure translation
        se3mat = [zeros(3), p;
                  0 0 0 0];
        return;
    end

    % SO(3) log
    so3mat = MatrixLog3(R);
    omgtheta = so3ToVec(so3mat);
    theta = norm(omgtheta);

    omgmat = so3mat / theta;

    G_inv = eye(3)/theta - 0.5*omgmat + ...
        (1/theta - 0.5*cot(theta/2))*(omgmat*omgmat);

    v = G_inv * p;
    se3mat = [so3mat, v*theta;
              0 0 0 0];
end

%% --- Matrix Exp of SE(3): takes full 4x4 se(3) matrix ---
function T = MatrixExp6_SE3(se3mat)
    so3mat = se3mat(1:3,1:3);
    vtheta = se3mat(1:3,4);

    if norm(so3mat, 'fro') < 1e-10
        T = [eye(3), vtheta;
             0 0 0 1];
        return;
    end

    omgtheta = so3ToVec(so3mat);
    theta = norm(omgtheta);

    omgmat = so3mat / theta;
    v = vtheta / theta;

    R = MatrixExp3(so3mat);
    G = eye(3)*theta + (1 - cos(theta))*omgmat + ...
        (theta - sin(theta))*(omgmat*omgmat);

    T = [R, G*v;
         0 0 0 1];
end

%% --- SO(3) helpers ---
function so3mat = MatrixLog3(R)
    acosinput = (trace(R) - 1)/2;
    acosinput = max(min(acosinput,1),-1);

    if abs(acosinput - 1) < 1e-10
        so3mat = zeros(3);
    elseif abs(acosinput + 1) < 1e-10
        theta = pi;
        if ~NearZero(1 + R(3,3))
            omg = (1/sqrt(2*(1 + R(3,3)))) * [R(1,3); R(2,3); 1 + R(3,3)];
        elseif ~NearZero(1 + R(2,2))
            omg = (1/sqrt(2*(1 + R(2,2)))) * [R(1,2); 1 + R(2,2); R(3,2)];
        else
            omg = (1/sqrt(2*(1 + R(1,1)))) * [1 + R(1,1); R(2,1); R(3,1)];
        end
        so3mat = VecToso3(omg * theta);
    else
        theta = acos(acosinput);
        so3mat = theta/(2*sin(theta)) * (R - R');
    end
end

function R = MatrixExp3(so3mat)
    omgtheta = so3ToVec(so3mat);
    theta = norm(omgtheta);

    if NearZero(theta)
        R = eye(3);
        return;
    end

    omgmat = so3mat / theta;
    R = eye(3) + sin(theta)*omgmat + (1 - cos(theta))*(omgmat*omgmat);
end

function so3mat = VecToso3(omg)
    so3mat = [    0   -omg(3)  omg(2);
               omg(3)    0    -omg(1);
              -omg(2) omg(1)    0   ];
end

function omg = so3ToVec(so3mat)
    omg = [so3mat(3,2); so3mat(1,3); so3mat(2,1)];
end

function V = se3ToVec(se3mat)
    V = [so3ToVec(se3mat(1:3,1:3)); se3mat(1:3,4)];
end

function flag = NearZero(z)
    flag = abs(z) < 1e-8;
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

%% --- URDF end-effector body detection ---
function eeName = get_end_effector_name(robot)
    % Try common Panda end-effector names first
    candidateNames = {'panda_hand_tcp', 'panda_hand', 'panda_link8'};

    bodyNames = robot.BodyNames;

    for i = 1:numel(candidateNames)
        if any(strcmp(bodyNames, candidateNames{i}))
            eeName = candidateNames{i};
            return;
        end
    end

    % Fallback: use last body in the chain
    eeName = bodyNames{end};
    warning('Could not find common Panda EE body; using last body: %s', eeName);
end