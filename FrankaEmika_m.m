%% Pt M: Bonus - Robotic System Toolbox Simulation
% rigid body tree object
% link lengths, m
close all;
L = [0.333 0.316 0.384 0.107];

% flange offset, m
A = 0.088;


% home position, https://frankarobotics.github.io/docs/robot_specifications.html#kinematic-configuration

M = [1 0 0 A;
    0 -1 0 0;
    0 0 -1 L(1)+L(2)+L(3)-L(4);
    0 0 0 1];

% Screw axis in space frame

ws = {[0;0;1], [0;-1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qs = {[0;0;0], [0;0;L(1)], [0;0;L(1)], [A;0;L(1)+L(2)], [0;0;L(1)+L(2)+L(3)], [0;0;L(1)+L(2)+L(3)], [A;0;L(1)+L(2)+L(3)-L(4)]};

robot = rigidBodyTree('DataFormat', 'row', 'MaxNumBodies', 8);

% Link lengths
% L = [0.333 0.316 0.384 0.107], A = 0.088
L1 = L(1); L2 = L(2); L3 = L(3); L4 = L(4);

%  7 joints of the FR3
% body is attached to the previous one with a transformation
bodies = cell(7,1);
joints = cell(7,1);

% kinematic chain 
% joint offsets to match L and A variables
joint_configs = {
    [0 0 0], 'revolute';           % J1: Z-axis at base
    [0 0 L1], 'revolute';          % J2: Y-axis
    [0 0 0], 'revolute';           % J3: Z-axis
    [A 0 L2], 'revolute';          % J4: Y-axis
    [0 0 0], 'revolute';           % J5: Z-axis
    [0 0 L3], 'revolute';          % J6: Y-axis
    [A 0 -L4], 'revolute'          % J7: Z-axis
};

% Axis definitions corresponding to 'ws'
joint_axes = [0 0 1; 0 -1 0; 0 0 1; 0 1 0; 0 0 1; 0 1 0; 0 0 1];

for i = 1:7
    bodies{i} = rigidBody(['link' num2str(i)]);
    joints{i} = rigidBodyJoint(['joint' num2str(i)], joint_configs{i,2});
    
    % joint axis
    joints{i}.JointAxis = joint_axes(i,:);
    
    % transform from parent to this joint
    trans = joint_configs{i,1};
    setFixedTransform(joints{i}, trvec2tform(trans));
    
    bodies{i}.Joint = joints{i};
    if i == 1
        addBody(robot, bodies{i}, 'base');
    else
        addBody(robot, bodies{i}, ['link' num2str(i-1)]);
    end
end

% Visualization and Comparison
% Toolbox wants 1x7 row vector
theta = [0 0 0 0 0 0 0];

config = theta; 

% Compute Toolbox FK
T_toolbox = getTransform(robot, config, 'link7');


fprintf('\n Robotic Toolbox Bonus \n');

% Graphical Simulation Window
figure('Name', 'Franka Emika FR3 Simulation');
show(robot, config, 'Frames', 'on', 'PreservePlot', false);
title('FR3 Robot Configuration (Robotic System Toolbox)');
axis equal;
view(3);
grid on;

%% Franka Emika
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

%% 2. Packing to Ready
fprintf('Animating IK: Packing to Ready...\n');
% Added 'robot' as the first argument
animate_transition_IK(robot, T_packing, T_ready, pose_packing, S_space, M, steps, pause_time, eomg, ev);

%% 3. Ready to Singularity
fprintf('Animating IK: Ready to Singularity...\n');
% Added 'robot' as the first argument
animate_transition_IK(robot, T_ready, T_singularity, pose_ready, S_space, M, steps, pause_time, eomg, ev);
%%
function animate_transition_IK(robot, T_start, T_end, theta_init, Slist, M, steps, pause_time, eomg, ev)
    % Precompute trajectory in SE(3)
    T_traj = cell(steps,1);
    for i = 1:steps
        s = (i-1)/(steps-1);
        % Standard SE(3) interpolation
        V = MatrixLog6(inv(T_start) * T_end);   
        se3mat = VecTose3(V);                  
        T_traj{i} = T_start * MatrixExp6(se3mat, s);
    end

    theta_guess = theta_init;
    
    % Setup Figure
    figure('Name', 'IK Animation');
    set(gcf, 'Color', 'w');
    
    for t = 1:steps
        % 1. Solve IK for current desired pose
        Tsd = T_traj{t};
        [theta_sol, success, w_err, v_err] = J_inverse_kinematics(Slist, M, Tsd, theta_guess, eomg, ev);
        
        % 2. Update guess for next iteration (temporal continuity)
        theta_guess = theta_sol;

        % 3. Visualization using the Toolbox
        % 'PreservePlot', false is key to prevent the figure from lagging
        show(robot, theta_guess, 'Frames', 'off', 'PreservePlot', false);
        
        hold on;
        % Optional: Plot the target point as a red star to see if the EE reaches it
        p_target = Tsd(1:3, 4);
        plot3(p_target(1), p_target(2), p_target(3), 'r*', 'MarkerSize', 10);
        
        title(sprintf('Step %d/%d | AngErr: %.4f | LinErr: %.4f', t, steps, w_err, v_err));
        axis([-0.8 0.8 -0.8 0.8 0 1.2]);
        view(3);
        grid on;
        drawnow;
        pause(pause_time);
        
        if t < steps, hold off; end % Clear target for next frame
    end
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