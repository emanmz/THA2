addpath("Functions");

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

% Screw axis in space frame

ws = {[0;0;1], [0;-1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qs = {[0;0;0], [0;0;L(1)], [0;0;L(1)], [A;0;L(1)+L(2)], [0;0;L(1)+L(2)+L(3)], [0;0;L(1)+L(2)+L(3)], [A;0;L(1)+L(2)+L(3)-L(4)]};

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

% Screw axis in body frame

wb = {[0;0;-1], [0;1;0], [0;0;-1], [0;-1;0], [0;0;-1], [0;-1;0], [0;0;-1]};
qb = {[-A;0;-L(4)+L(3)+L(2)], [-A;0;-L(4)+L(3)+L(2)], [-A;0;-L(4)+L(3)+L(2)], [0;0;-L(4)+L(3)], [-A;0;-L(4)], [-A;0;-L(4)], [0;0;0]};

S_body = zeros(6,7);

for i=1:7
    wi = wb{i};
    if norm(wi) == 0
        vi = qb{i};
    else
    vi = cross(-wi, qb{i});
    end
    
    S_body(:, i) = [wi; vi];
end
disp("S_body")
disp(S_body)


%% Pt A-C: Forward Kinematics

% arbitrary configurations
theta = [0 -pi/4 0 -3*pi/4 0 pi/2 pi/4]; % i dont know how to pick a test joint angles and we probably need a couple different ones... 

% Pt A: "manually" compute space FK
Ts_manual = screw_to_exp(S_space(:,1), theta(1))*screw_to_exp(S_space(:,2), theta(2))*...
    screw_to_exp(S_space(:,3), theta(3))*screw_to_exp(S_space(:,4), theta(4))*screw_to_exp(S_space(:,5), theta(5))*...
    screw_to_exp(S_space(:,6), theta(6))*screw_to_exp(S_space(:,7), theta(7))*M;
disp('Ts_manual')
disp(Ts_manual)

% Pt B: Compute space FK with function
Emika_FKS = FK_space(M, S_space, theta);
disp('Emika_FKS')
disp(Emika_FKS)

% Pt C: "manually" compute space FK
Tb_manual = M*screw_to_exp(S_body(:,1), theta(1))*screw_to_exp(S_body(:,2), theta(2))*...
    screw_to_exp(S_body(:,3), theta(3))*screw_to_exp(S_body(:,4), theta(4))*screw_to_exp(S_body(:,5), theta(5))*...
    screw_to_exp(S_body(:,6), theta(6))*screw_to_exp(S_body(:,7), theta(7));
disp('Tb_manual')
disp(Tb_manual)


% Pt C: Compute body FK with function
Emika_FKB = FK_body(M, S_body, theta);
disp('Emika_FKB')
disp(Emika_FKB)

%% Pt D: Jacobian

% Space Form
Emika_JS = J_space(S_space, theta);
disp('Emika_JS')
disp(Emika_JS)

% Body Form
Emika_JB = J_body(S_body, theta);
disp('Emika_JB')
disp(Emika_JB)

%% Pt F: Singularity

thetaSym = sym('theta', [1 7]);
% computationally heavy :P idk if we need to do something aboutthis? 
try
    Emika_sing = singularity(S_space, thetaSym); 
    disp('Singularity analysis initialized symbolically...');
catch
    disp('blah did not work');
end

%% Pt G: Manipulability Ellipsoids
figure('Name', 'FR3 Manipulability');
hold on; grid on; axis equal; view(3);
ellipsoid_plot_linear(Emika_JS, Emika_FKS); 
ellipsoid_plot_angular(Emika_JS, Emika_FKS);
title('FR3 Manipulability Ellipsoids (Linear: Magenta, Angular: Green)');

%% Pt H: Inverse Kinematics
T_target = Emika_FKS;
T_target(1,4) = T_target(1,4) + 0.05; % Move 5cm in X

[theta_ik, success_ik] = J_inverse_kinematics(S_space, M, T_target, theta);
fprintf('Standard IK Success: %d\n', success_ik);

%% Pt i: Jacobian Transpose Algorithm
[theta_trans, success_trans] = J_transpose_kinematics(S_space, M, T_target, theta);
fprintf('Jacobian Transpose Success: %d\n', success_trans);

%% Pt j: Redundancy Resolution
[theta_red, success_red] = redundancy_resolution(S_space, M, T_target, theta);
fprintf('Redundancy Resolution Success: %d\n', success_red);

%% Pt k: Bonus - Damped Least Squares (DLS)
% if the target is near a boundary
[theta_dls, success_dls] = DLS_inverse_kinematics(S_space, M, T_target, theta);
fprintf('DLS IK Success: %d\n', success_dls);


%% Pt M: Bonus - Robotic System Toolbox Simulation
% rigid body tree object
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
config = theta; 

% Compute Toolbox FK
T_toolbox = getTransform(robot, config, 'link7');

% with PoE results
error_toolbox = norm(Emika_FKS - T_toolbox);

fprintf('\n--- Robotic System Toolbox Bonus ---\n');
fprintf('Toolbox FK Equivalence Error: %e\n', error_toolbox);

% Graphical Simulation Window
figure('Name', 'Franka Emika FR3 Simulation');
show(robot, config, 'Frames', 'on', 'PreservePlot', false);
title('FR3 Robot Configuration (Robotic System Toolbox)');
axis equal;
view(3);
grid on;

% Optional: Interaction GUI
gui = interactiveRigidBodyTree(robot);