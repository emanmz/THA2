%% Pt M: Bonus - Robotic System Toolbox Simulation
% rigid body tree object
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

fprintf('\n Robotic Toolbox Bonus \n');
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