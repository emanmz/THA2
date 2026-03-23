%% Test Function for Space Jacobian 
% % PA pt. e & f test functions 
addpath("Functions");

%% SCARA Example (W7-L1 Slide 7-9)

% test values 
theta = [pi/2, pi/2, pi/2, pi/2];
L1 = 2;
L2 = 2;

% Space Jacobian from slide 9

Js_SCARA = [ 0 0 0 0;
            0 0 0 0;
            1 1 1 0;
            0  L1*sin(theta(1))  L1*sin(theta(1))+L2*sin(theta(1)+theta(2))  0;
            0  -L1*cos(theta(1))  -L1*cos(theta(1))-L2*cos(theta(1)+theta(2))  0;
            0 0 0 1];

disp("SCARA Js")
disp(Js_SCARA)

%% J_space Test - SCARA

% Make Screw Axis in Space Frame (W5-L1 Slide 8)
ws = {[0;0;1], [0;0;1], [0;0;1], [0;0;0]};
qs = {[0;0;0], [L1;0;0], [L1 + L2;0;0], [0;0;1]};

S_space = zeros(6,4);

for i=1:4
    wi = ws{i};
    if norm(wi) == 0
        vi = qs{i};
    else
    vi = cross(-wi, qs{i});
    end
    
    S_space(:, i) = [wi; vi];
end

% test function 
Js = J_space(S_space, theta);
disp("J_space")
disp(Js)

%% J_body Test - SCARA

% Make Screw Axis in Space Frame (W5-L1 Slide 8)
ws = {[0;0;1], [0;0;1], [0;0;1], [0;0;0]};
qs = {[-L1-L2;0;0], [-L2;0;0], [0;0;0], [0;0;1]};

S_body = zeros(6,4);

for i=1:4
    wi = ws{i};
    if norm(wi) == 0
        vi = qs{i};
    else
    vi = cross(-wi, qs{i});
    end
    
    S_body(:, i) = [wi; vi];
end

% test function 
Jb = J_body(S_body, theta);
disp("J_body")
disp(Jb)

% convert body to space
M = [1 0 0 L1+L2;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

T_sb = FK_space(M, S_space, theta);
Ad_sb = Adjoint(T_sb);
Js_check = Ad_sb*Jb;

%% Singulairty function test
% 
thetaSym = sym('theta', [1 4]);

Sings = singularity(S_space, thetaSym);
%% Test Function for Space Jacobian 
addpath("Functions");
% example from w7-L1 SL 9
theta = [pi/2, pi/2, pi/2, pi/2];
L1 = 2;
L2 = 2;


ws = {[0;0;1], [0;0;1], [0;0;1], [0;0;0]};
qs = {[0;0;0], [L1*cos(theta(1));L1*sin(theta(1));0], [L1*cos(theta(1))+L2*(cos(theta(1)+theta(2))); L1*sin(theta(1))+L2*(sin(theta(1)+theta(2))); 0]...
    [0;0;1]};

J_RRRP = zeros(6,4);

for i=1:4
    wi = ws{i};

    if norm(wi) == 0
        vi = qs{i};
    else
    vi = cross(-wi, qs{i});
    end
    
    J_RRRP(:, i) = [wi; vi];
end
disp("J_RRRP Example")
disp(J_RRRP)

%% home position
ws = {[0;0;1], [0;0;1], [0;0;1], [0;0;0]};
qs = {[0;0;0], [L1;0;0], [L1 + L2;0;0], [0;0;1]};

S_RRRP = zeros(6,4);

for i=1:4
    wi = ws{i};
    if norm(wi) == 0
        vi = qs{i};
    else
    vi = cross(-wi, qs{i});
    end
    
    S_RRRP(:, i) = [wi; vi];
end
% disp("S Home Positions")
% disp(S_RRRP)


Js = J_space(S_RRRP, theta);
disp("J_space")
disp(Js)
