close all; clc; clear;

%% =========================
%  COMPARABLE IK vs JT vs RR vs DLS TEST
%  Same robot, same path, same tolerances, same steps
%  =========================

%% Link lengths (m)
L = [0.333 0.316 0.384 0.107];
A = 0.088;

%% Screw axis definitions
ws = {[0;0;1], [0;-1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qs = {[0;0;0], [0;0;L(1)], [0;0;L(1)], [A;0;L(1)+L(2)], ...
      [0;0;L(1)+L(2)+L(3)], [0;0;L(1)+L(2)+L(3)], ...
      [A;0;L(1)+L(2)+L(3)-L(4)]};

%% Home configuration
M = [1 0 0 A;
     0 -1 0 0;
     0 0 -1 L(1)+L(2)+L(3)-L(4);
     0 0 0 1];

%% Space screw axes
S_space = zeros(6,7);
for i = 1:7
    wi = ws{i};
    vi = cross(-wi, qs{i});
    S_space(:, i) = [wi; vi];
end

%% Configurations
pose_home    = zeros(1,7);
pose_ready   = [0, -0.4, 0, -2.0, 0, 1.5, 0.7];
pose_packing = [0, 0.5, 0, -2.5, 0, 1.0, 0];

%% Choose ONE identical motion for all methods
% Example: Packing -> Ready
theta_start = pose_packing;
theta_goal  = pose_ready;

T_start = FK_space_no_plot(M, S_space, theta_start);
T_end   = FK_space_no_plot(M, S_space, theta_goal);

%% Common parameters (IDENTICAL for all)
steps      = 200;
pause_time = 0.01;
eomg       = 1e-3;
ev         = 1e-3;
eomg_T       = 1e-2;
ev_T         = 1e-2;
max_iter   = 100;

%% =========================
%  RUN ALL METHODS
%  =========================
fprintf('\n=== Running Pseudoinverse IK (Packing -> Ready) ===\n');
results_IK = animate_transition_common( ...
    'IK', T_start, T_end, theta_start, S_space, M, qs, ...
    steps, pause_time, eomg, ev, max_iter, 'ik_comparison.gif');

fprintf('\n=== Running Jacobian Transpose (Packing -> Ready) ===\n');
results_JT = animate_transition_common( ...
    'JT', T_start, T_end, theta_start, S_space, M, qs, ...
    steps, pause_time, eomg_T, ev_T, max_iter, 'jt_comparison.gif');

fprintf('\n=== Running Redundancy Resolution (Packing -> Ready) ===\n');
results_RR = animate_transition_common( ...
    'RR', T_start, T_end, theta_start, S_space, M, qs, ...
    steps, pause_time, eomg, ev, max_iter, 'rr_comparison.gif');

fprintf('\n=== Running Damped Least Squares (Packing -> Ready) ===\n');
results_DLS = animate_transition_common( ...
    'DLS', T_start, T_end, theta_start, S_space, M, qs, ...
    steps, pause_time, eomg, ev, max_iter, 'dls_comparison.gif');

%% Summary comparison
print_comparison_all({results_IK, results_JT, results_RR, results_DLS});

%% =========================================================
%  COMMON ANIMATION / DRIVER
%  =========================================================
function results = animate_transition_common(method, T_start, T_end, theta_init, ...
    Slist, M, qs, steps, pause_time, eomg, ev, max_iter, gif_filename)

    % Cartesian interpolation (same for all)
    T_rel   = inv(T_start) * T_end;
    V_rel   = MatrixLog6(T_rel);
    se3_rel = VecTose3(V_rel);

    theta_guess = theta_init;

    % -------- Figure 1: Animation figure --------
    h_fig = figure('Color','w','WindowState','maximized', ...
                   'Name', [upper(method) ' Animation']);

    % -------- Figure 2: Metrics figure (separate live plot) --------
    h_metrics = figure('Color','w','Name', [upper(method) ' Metrics']);
    tiledlayout(2,1);

    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on');
    title(ax1, [upper(method) ' Condition Number']);
    xlabel(ax1, 'Path Step'); ylabel(ax1, 'cond(J)');
    h_cond = plot(ax1, nan, nan, 'LineWidth', 1.8);

    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on');
    title(ax2, [upper(method) ' Isotropy']);
    xlabel(ax2, 'Path Step'); ylabel(ax2, 'Isotropy');
    h_iso = plot(ax2, nan, nan, '--', 'LineWidth', 1.8);

    % Preallocate histories
    cond_hist    = nan(steps,1);
    iso_hist     = nan(steps,1);
    mu_lin_hist  = nan(steps,1);
    mu_ang_hist  = nan(steps,1);
    iter_hist    = nan(steps,1);
    time_hist    = nan(steps,1);
    pos_err_hist = nan(steps,1);
    ang_err_hist = nan(steps,1);
    lin_err_hist = nan(steps,1);
    succ_hist    = false(steps,1);

    for t = 1:steps
        s   = (t-1)/(steps-1);
        Tsd = T_start * MatrixExp6(se3_rel, s);

        % Solve one waypoint with chosen method
        switch upper(method)
            case 'IK'
                [theta_sol, success, info, solve_time] = IK_solver_step( ...
                    Slist, M, Tsd, theta_guess, eomg, ev, max_iter);

            case 'JT'
                [theta_sol, success, info, solve_time] = JT_solver_step( ...
                    Slist, M, Tsd, theta_guess, eomg, ev, max_iter);

            case 'RR'
                [theta_sol, success, info, solve_time] = RR_solver_step( ...
                    Slist, M, Tsd, theta_guess, eomg, ev, max_iter);

            case 'DLS'
                [theta_sol, success, info, solve_time] = DLS_solver_step( ...
                    Slist, M, Tsd, theta_guess, eomg, ev, max_iter);

            otherwise
                error('Unknown method: %s', method);
        end

        theta_guess = theta_sol;

        % Current state
        T_curr = FK_space_no_plot(M, Slist, theta_guess);
        Js = J_space(Slist, theta_guess);

        Jw = Js(1:3,:);   % angular manipulability ellipsoid
        Jv = Js(4:6,:);   % linear manipulability ellipsoid

        % Metrics
        svals = svd(Js);
        sigma_max = max(svals);
        sigma_min = min(svals);

        if sigma_min < 1e-12
            condJ = inf;
            isoJ  = 0;
        else
            condJ = sigma_max / sigma_min;
            isoJ  = sigma_min / sigma_max;
        end

        % Yoshikawa measures
        A_lin = Jv * Jv.';
        mu_lin = sqrt(max(det(A_lin), 0));

        A_ang = Jw * Jw.';
        mu_ang = sqrt(max(det(A_ang), 0));

        % Store histories
        cond_hist(t)     = condJ;
        iso_hist(t)      = isoJ;
        mu_lin_hist(t)   = mu_lin;
        mu_ang_hist(t)   = mu_ang;
        iter_hist(t)     = info.iters;
        time_hist(t)     = solve_time;
        pos_err_hist(t)  = norm(T_curr(1:3,4) - Tsd(1:3,4));
        ang_err_hist(t)  = info.w_err;
        lin_err_hist(t)  = info.v_err;
        succ_hist(t)     = success;

        % -------- Animation Plot --------
        figure(h_fig);
        clf;

        subplot(2,2,[1 3]);
        hold on; grid on; axis equal; view(3);
        axis([-0.8 0.8 -0.8 0.8 -0.1 1.2]);

        draw_arm_fixed(Slist, M, theta_guess, qs);
        plot3(Tsd(1,4), Tsd(2,4), Tsd(3,4), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
        plot3(T_start(1,4), T_start(2,4), T_start(3,4), 'go', 'MarkerSize', 8, 'LineWidth', 1.5);
        plot3(T_end(1,4),   T_end(2,4),   T_end(3,4),   'bx', 'MarkerSize', 10, 'LineWidth', 2);

        % Better status labeling
        if success
            status_str = 'OK';
        elseif pos_err_hist(t) < 1e-2
            status_str = 'CLOSE (NO CONV)';
        else
            status_str = 'NO CONV';
        end

        title(sprintf(['%s | Step %d/%d | %s\n' ...
                       'iters=%d | step time=%.4fs\n' ...
                       'w err=%.2e | v err=%.2e | pos err=%.2e'], ...
                       upper(method), t, steps, status_str, ...
                       info.iters, solve_time, info.w_err, info.v_err, pos_err_hist(t)));

        % Linear manipulability ellipsoid
        subplot(2,2,2);
        hold on; grid on; axis equal; view(3);
        plot_ellipsoid(Jv);
        title(sprintf('Linear Manipulability Ellipsoid | \\mu_v = %.3e', mu_lin));
        xlabel('x'); ylabel('y'); zlabel('z');

        % Angular manipulability ellipsoid
        subplot(2,2,4);
        hold on; grid on; axis equal; view(3);
        plot_ellipsoid(Jw);
        title(sprintf('Angular Manipulability Ellipsoid | \\mu_\\omega = %.3e', mu_ang));
        xlabel('\omega_x'); ylabel('\omega_y'); zlabel('\omega_z');

        drawnow;

        % -------- Update separate metrics figure --------
        if isvalid(h_metrics)
            set(h_cond, 'XData', 1:t, 'YData', cond_hist(1:t));
            set(h_iso,  'XData', 1:t, 'YData', iso_hist(1:t));

            xlim(ax1, [1 steps]);
            xlim(ax2, [1 steps]);

            if any(isfinite(cond_hist(1:t)))
                ymax = max(cond_hist(isfinite(cond_hist(1:t))));
                ylim(ax1, [0, max(1, 1.05*ymax)]);
            end
            ylim(ax2, [0, 1]);

            drawnow limitrate;
        end

        % GIF
        frame = getframe(h_fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        if t == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', pause_time);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', pause_time);
        end

        pause(pause_time);
    end

    fprintf('%s GIF saved as: %s\n', upper(method), gif_filename);

    % Return results
    results.method        = upper(method);
    results.cond_hist     = cond_hist;
    results.iso_hist      = iso_hist;
    results.mu_lin_hist   = mu_lin_hist;
    results.mu_ang_hist   = mu_ang_hist;
    results.iter_hist     = iter_hist;
    results.time_hist     = time_hist;
    results.pos_err_hist  = pos_err_hist;
    results.ang_err_hist  = ang_err_hist;
    results.lin_err_hist  = lin_err_hist;
    results.succ_hist     = succ_hist;
    results.final_theta   = theta_guess;
    results.total_time    = sum(time_hist);
    results.avg_time      = mean(time_hist);
    results.avg_iters     = mean(iter_hist);
    results.max_cond      = max(cond_hist(isfinite(cond_hist)));
    results.min_iso       = min(iso_hist);
    results.success_rate  = mean(succ_hist) * 100;
    results.final_pos_err = pos_err_hist(end);
end

%% =========================================================
%  PSEUDOINVERSE IK SOLVER
%  =========================================================
function [theta, success, info, solve_time] = IK_solver_step(Slist, M, Tsd, theta0, eomg, ev, max_iter)
    tic;
    theta = theta0;
    i = 0;

    Tsb = FK_space_no_plot(M, Slist, theta);
    [Vb, ~] = calculate_twist_error(Tsb, Tsd);
    w_err = norm(Vb(1:3));
    v_err = norm(Vb(4:6));

    while (w_err > eomg || v_err > ev) && i < max_iter
        Js = J_space(Slist, theta);

        % body error -> space error
        Vs = Adjoint(Tsb) * Vb;

        % Slightly regularized pseudoinverse for stability
        theta = theta + (pinv(Js, 1e-6) * Vs).';

        i = i + 1;

        Tsb = FK_space_no_plot(M, Slist, theta);
        [Vb, ~] = calculate_twist_error(Tsb, Tsd);
        w_err = norm(Vb(1:3));
        v_err = norm(Vb(4:6));
    end

    success = (w_err <= eomg && v_err <= ev);

    info.w_err = w_err;
    info.v_err = v_err;
    info.iters = i;
    solve_time = toc;
end

%% =========================================================
%  IMPROVED JACOBIAN TRANSPOSE SOLVER
%  Adaptive step size (much better than constant alpha)
%  =========================================================
function [theta, success, info, solve_time] = JT_solver_step(Slist, M, Tsd, theta0, eomg, ev, max_iter)
    tic;
    theta = theta0;
    i = 0;

    Tsb = FK_space_no_plot(M, Slist, theta);
    [Vb, ~] = calculate_twist_error(Tsb, Tsd);
    w_err = norm(Vb(1:3));
    v_err = norm(Vb(4:6));

    while (w_err > eomg || v_err > ev) && i < max_iter
        Js = J_space(Slist, theta);

        % body error -> space error
        Vs = Adjoint(Tsb) * Vb;

        % Adaptive steepest descent gain:
        % alpha = (e^T J J^T e)/( (J J^T e)^T (J J^T e) )
        JJTe = Js * (Js.' * Vs);
        denom = dot(JJTe, JJTe);

        if denom < 1e-12
            alpha = 0.01;
        else
            alpha = dot(Vs, JJTe) / denom;
        end

        % Clamp alpha to prevent blow-up / overstepping
        alpha = max(min(alpha, 1.0), 1e-4);

        theta = theta + (alpha * Js.' * Vs).';

        i = i + 1;

        Tsb = FK_space_no_plot(M, Slist, theta);
        [Vb, ~] = calculate_twist_error(Tsb, Tsd);
        w_err = norm(Vb(1:3));
        v_err = norm(Vb(4:6));
    end

    success = (w_err <= eomg && v_err <= ev);

    info.w_err = w_err;
    info.v_err = v_err;
    info.iters = i;
    solve_time = toc;
end

%% =========================================================
%  REDUNDANCY RESOLUTION SOLVER
%  Pseudoinverse + Null-space manipulability maximization
%  =========================================================
function [theta, success, info, solve_time] = RR_solver_step(Slist, M, Tsd, theta0, eomg, ev, max_iter)
    tic;
    theta = theta0;
    i = 0;
    k = 0.05;     % secondary task gain
    eps_fd = 1e-4;

    Tsb = FK_space_no_plot(M, Slist, theta);
    [Vb, ~] = calculate_twist_error(Tsb, Tsd);
    w_err = norm(Vb(1:3));
    v_err = norm(Vb(4:6));

    while (w_err > eomg || v_err > ev) && i < max_iter
        Js = J_space(Slist, theta);
        Vs = Adjoint(Tsb) * Vb;

        J_pinv = pinv(Js, 1e-6);

        % Gradient of manipulability (full Jacobian volume)
        grad_w = zeros(length(theta), 1);
        for j = 1:length(theta)
            th_plus = theta;
            th_plus(j) = th_plus(j) + eps_fd;
            w_plus = J_ellipsoid_volume(J_space(Slist, th_plus));

            th_minus = theta;
            th_minus(j) = th_minus(j) - eps_fd;
            w_minus = J_ellipsoid_volume(J_space(Slist, th_minus));

            grad_w(j) = (w_plus - w_minus) / (2 * eps_fd);
        end

        % Normalize gradient so it doesn't overpower primary task
        if norm(grad_w) > 1e-10
            grad_w = grad_w / norm(grad_w);
        end

        N = eye(size(Js,2)) - J_pinv * Js;
        dtheta = J_pinv * Vs + k * N * grad_w;

        theta = theta + dtheta.';

        i = i + 1;

        Tsb = FK_space_no_plot(M, Slist, theta);
        [Vb, ~] = calculate_twist_error(Tsb, Tsd);
        w_err = norm(Vb(1:3));
        v_err = norm(Vb(4:6));
    end

    success = (w_err <= eomg && v_err <= ev);

    info.w_err = w_err;
    info.v_err = v_err;
    info.iters = i;
    solve_time = toc;
end

%% =========================================================
%  DAMPED LEAST SQUARES SOLVER
%  Better near singularities
%  =========================================================
function [theta, success, info, solve_time] = DLS_solver_step(Slist, M, Tsd, theta0, eomg, ev, max_iter)
    tic;
    theta = theta0;
    i = 0;

    lambda_max = 0.1;   % max damping
    epsilon    = 0.05;  % singularity neighborhood

    Tsb = FK_space_no_plot(M, Slist, theta);
    [Vb, ~] = calculate_twist_error(Tsb, Tsd);
    w_err = norm(Vb(1:3));
    v_err = norm(Vb(4:6));

    while (w_err > eomg || v_err > ev) && i < max_iter
        Js = J_space(Slist, theta);
        Vs = Adjoint(Tsb) * Vb;

        s = svd(Js);
        s_min = min(s);

        % Adaptive damping
        if s_min < epsilon
            lambda = (1 - (s_min/epsilon)^2) * lambda_max^2;
        else
            lambda = 0;
        end

        % DLS inverse
        J_dls = Js' / (Js * Js' + lambda * eye(6));

        theta = theta + (J_dls * Vs).';

        i = i + 1;

        Tsb = FK_space_no_plot(M, Slist, theta);
        [Vb, ~] = calculate_twist_error(Tsb, Tsd);
        w_err = norm(Vb(1:3));
        v_err = norm(Vb(4:6));
    end

    success = (w_err <= eomg && v_err <= ev);

    info.w_err = w_err;
    info.v_err = v_err;
    info.iters = i;
    solve_time = toc;
end

%% =========================================================
%  FINAL COMPARISON PRINT (ALL METHODS)
%  =========================================================
function print_comparison_all(results_cell)
    n = length(results_cell);

    fprintf('\n================================================================================================================\n');
    fprintf('FINAL COMPARISON: ALL METHODS\n');
    fprintf('================================================================================================================\n');
    fprintf('%-22s', 'Metric');
    for k = 1:n
        fprintf(' | %-12s', results_cell{k}.method);
    end
    fprintf('\n');
    fprintf(repmat('-',1,22 + 15*n));
    fprintf('\n');

    print_metric_row('Total Time (s)',      results_cell, 'total_time');
    print_metric_row('Avg Step Time (s)',   results_cell, 'avg_time');
    print_metric_row('Avg Internal Iters',  results_cell, 'avg_iters');
    print_metric_row('Success Rate (%)',    results_cell, 'success_rate');
    print_metric_row('Final Pos Error',     results_cell, 'final_pos_err', true);
    print_metric_row('Worst cond(J)',       results_cell, 'max_cond');
    print_metric_row('Worst Isotropy',      results_cell, 'min_iso');
    print_metric_mean_row('Avg Linear Mu',  results_cell, 'mu_lin_hist', true);
    print_metric_mean_row('Avg Angular Mu', results_cell, 'mu_ang_hist', true);

    fprintf('================================================================================================================\n');

    % Overlay plots
    figure('Color','w','Name','All Method Comparison');
    tiledlayout(3,1);

    nexttile; hold on; grid on;
    for k = 1:n
        plot(results_cell{k}.cond_hist, 'LineWidth', 1.6);
    end
    ylabel('cond(J)');
    title('Condition Number Comparison');
    legend(cellfun(@(r) r.method, results_cell, 'UniformOutput', false), 'Location', 'best');

    nexttile; hold on; grid on;
    for k = 1:n
        plot(results_cell{k}.iso_hist, 'LineWidth', 1.6);
    end
    ylabel('Isotropy');
    title('Isotropy Comparison');
    legend(cellfun(@(r) r.method, results_cell, 'UniformOutput', false), 'Location', 'best');

    nexttile; hold on; grid on;
    for k = 1:n
        plot(results_cell{k}.iter_hist, 'LineWidth', 1.6);
    end
    ylabel('Internal Iters');
    xlabel('Path Step');
    title('Solver Iterations per Path Step');
    legend(cellfun(@(r) r.method, results_cell, 'UniformOutput', false), 'Location', 'best');
end

function print_metric_row(label, results_cell, fieldname, scientific)
    if nargin < 4, scientific = false; end
    fprintf('%-22s', label);
    for k = 1:length(results_cell)
        val = results_cell{k}.(fieldname);
        if scientific
            fprintf(' | %-12.3e', val);
        else
            fprintf(' | %-12.4f', val);
        end
    end
    fprintf('\n');
end

function print_metric_mean_row(label, results_cell, fieldname, scientific)
    if nargin < 4, scientific = false; end
    fprintf('%-22s', label);
    for k = 1:length(results_cell)
        val = mean(results_cell{k}.(fieldname));
        if scientific
            fprintf(' | %-12.3e', val);
        else
            fprintf(' | %-12.4f', val);
        end
    end
    fprintf('\n');
end

%% =========================================================
%  GEOMETRY / KINEMATICS HELPERS
%  =========================================================
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

    plot3(joint_coords(1,:), joint_coords(2,:), joint_coords(3,:), '-o', ...
        'LineWidth', 3, 'Color', [0.2 0.2 0.2], 'MarkerFaceColor', 'k');
    plot3(joint_coords(1,end), joint_coords(2,end), joint_coords(3,end), ...
        'rs', 'MarkerFaceColor', 'r');
end

function T = FK_space_no_plot(M, Slist, theta)
    T = eye(4);
    for i = 1:length(theta)
        T = T * MatrixExp6(VecTose3(Slist(:,i)), theta(i));
    end
    T = T * M;
end

function Js = J_space(Slist, theta)
    n = size(Slist,2);
    Js = zeros(6,n);
    Js(:,1) = Slist(:,1);
    T = eye(4);
    for i = 2:n
        T = T * MatrixExp6(VecTose3(Slist(:,i-1)), theta(i-1));
        Js(:,i) = Adjoint(T) * Slist(:,i);
    end
end

function AdT = Adjoint(T)
    R = T(1:3,1:3);
    p_hat = VecToso3(T(1:3,4));
    AdT = [R, zeros(3); p_hat*R, R];
end

function [Vb, err] = calculate_twist_error(Tsb, Tsd)
    T_err = inv(Tsb) * Tsd;
    Vb = MatrixLog6(T_err);
    err = norm(Vb);
end

function V = MatrixLog6(T)
    R = T(1:3,1:3);
    p = T(1:3,4);

    if norm(R - eye(3), 'fro') < 1e-10
        V = [0;0;0; p];
        return;
    end

    c = max(min((trace(R)-1)/2, 1), -1);
    theta = acos(c);

    if abs(theta) < 1e-12
        V = [0;0;0; p];
        return;
    end

    so3mat = theta/(2*sin(theta)) * (R - R');
    w_hat = so3mat / theta;
    G_inv = eye(3)/theta - 0.5*w_hat + (1/theta - 0.5*cot(theta/2))*(w_hat*w_hat);

    V = [so3ToVec(so3mat); G_inv * p * theta];
end

function T = MatrixExp6(se3mat, theta)
    omg_skew = se3mat(1:3,1:3);
    v = se3mat(1:3,4);

    omg_vec = so3ToVec(omg_skew);
    ang = norm(omg_vec);

    if ang < 1e-12
        T = [eye(3), v*theta; 0 0 0 1];
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
    so3mat = [  0   -w(3)  w(2);
              w(3)   0    -w(1);
             -w(2)  w(1)   0   ];
end

function w = so3ToVec(so3mat)
    w = [so3mat(3,2); so3mat(1,3); so3mat(2,1)];
end

function plot_ellipsoid(J)
    A = J * J.';
    A = (A + A.')/2; % enforce symmetry
    [V,D] = eig(A);
    radii = sqrt(max(diag(D), 0));

    [x,y,z] = sphere(20);
    xyz = [x(:) y(:) z(:)]';
    ellipsoid_pts = V * diag(radii) * xyz;

    x_e = reshape(ellipsoid_pts(1,:), size(x));
    y_e = reshape(ellipsoid_pts(2,:), size(y));
    z_e = reshape(ellipsoid_pts(3,:), size(z));

    surf(x_e, y_e, z_e, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    camlight; lighting gouraud;
end

function mu = J_ellipsoid_volume(J)
    A = J * J.';
    A = (A + A.')/2;
    mu = sqrt(max(det(A), 0));
end