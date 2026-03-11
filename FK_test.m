%% Barrett WAM arm W6-L2-SL 3-4 

L1 = .550; % m
L2 = .300;
L3 = .060;
W1 = .045;

theta = [0, pi/4, 0, -pi/4, 0, -pi/2, 0];

% Body 
wb = {[0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qb = {[0;0;0], [0;0;-(L1+L2+L3)], [0;0;0], [W1;0;-(L2+L3)], [0;0;0], [0;0;-L3], [0;0;0]};

% Space
ws = {[], [], [], [], [], [], []};
qs = {[], [], [], [], [], [], []};


M = [1 0 0 0; 
     0 1 0 0;
     0 0 1 L1+L2+L3;
     0 0 0 1];


Tb_WAM = FK_body(M, wb, qb, theta)
Ts_WAM = FK_space(M, ws, qs, theta)

% also check if space and body are equivalent using adjoint. W6-L1-SL15 i
% dont know how to extract the relevant stuff. 

%% FK Test Cases

% W6-L1-SL11 3R Example

L1 = 2;
L2 = 3;
M = [0 0 1 L1; 0 1 0 0; -1 0 0 -L2; 0 0 0 1];

w = {[0;0;1], [0;-1;0], [1;0;0]};
q = {[0;0;0], [0;L1;0], [0;0;-L2]};
theta = [pi/2, pi, 3*pi/2];

Ts = FK_space(M, w, q, theta);
Tb = FK_body(M, w, q, theta);

% add the adjoint check eventually...