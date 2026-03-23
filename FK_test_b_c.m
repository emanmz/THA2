% Test function for part B & C
addpath("Functions");
close all; clear all; clc;
%% W6-L1-SL11 3R Example

% 3R dimensions
L1 = 2;
L2 = 3;
M_3r = [0 0 1 L1; 0 1 0 0; -1 0 0 -L2; 0 0 0 1];

% Space Form
ws_3r = {[0;0;1], [0;-1;0], [1;0;0]};
qs_3r = {[0;0;0], [L1;0;0], [0;0;-L2]};

% Convert to Slist
Slist_3r = zeros(6, 3);
for i = 1:3
    wi = ws_3r{i};
    vi = cross(-wi, qs_3r{i}); % Space screw v calculation
    Slist_3r(:, i) = [wi; vi];
end

% set to home position
theta_3rH = [0, 0, 0 ];
Ts_3rH = FK_space(M_3r, Slist_3r, theta_3rH);
disp('Ts_3r Home Position')
disp(Ts_3rH)

% set to simple position
theta_3rS = [0, pi/2, 0 ];
Ts_3rS = FK_space(M_3r, Slist_3r, theta_3rS);
disp('Ts_3r Simple Position')
disp(Ts_3rS)

% Derive Blist for 3R
SBlist_3r = zeros(6, 3);
for i = 1:3
    SBlist_3r(:, i) = Adjoint(inv(M_3r)) * Slist_3r(:, i);
end

% set to home position
Tb_3rH = FK_body(M_3r, SBlist_3r, theta_3rH);
disp('Tb_3r Home Position')
disp(Tb_3rH)

% set to simple position
theta_3rS = [0, pi/2, 0 ];
Tb_3rS = FK_body(M_3r, SBlist_3r, theta_3rS);
disp('Tb_3r Simple Position')
disp(Tb_3rS)

fprintf('3R Space/Body Equivalence Error, Home Position: %e\n', norm(Ts_3rH - Tb_3rH));
fprintf('3R Space/Body Equivalence Error, Simple position: %e\n', norm(Ts_3rS - Tb_3rS));


%% Barrett WAM arm W6-L2-SL 3-4 

L1 = .550; % m
L2 = .300;
L3 = .060;
W1 = .045;

theta_wam = [0, pi/4, 0, -pi/4, 0, -pi/2, 0];

% Body Screws: Convert {wb} and {qb} into a 6x7 Blist
wb = {[0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1], [0;1;0], [0;0;1]};
qb = {[0;0;0], [0;0;-(L1+L2+L3)], [0;0;0], [W1;0;-(L2+L3)], [0;0;0], [0;0;-L3], [0;0;0]};

% Home position
M_wam = [1 0 0 0; 
         0 1 0 0;
         0 0 1 L1+L2+L3;
         0 0 0 1];
Blist_wam = zeros(6, 7);
for i = 1:7
    wi = wb{i};
    vi = cross(-wi, qb{i}); % Body screw v calculation
    Blist_wam(:, i) = [wi; vi];
end


% Forward Kinematics (Body)
Tb_WAM = FK_body(M_wam, Blist_wam, theta_wam);
disp('Tb_WAM')
disp(Tb_WAM)

%  derive Slist from Blist using the Adjoint of M
Slist_wam = zeros(6, 7);
for i = 1:7
    Slist_wam(:, i) = Adjoint(M_wam) * Blist_wam(:, i);
end

% Forward Kinematics (Space)

Ts_WAM = FK_space(M_wam, Slist_wam, theta_wam);
disp('Ts_WAM')
disp(Ts_WAM)

% Check equivalence
error_wam = norm(Ts_WAM - Tb_WAM);
fprintf('WAM Space/Body Equivalence Error: %e\n', error_wam);

%% W6-L1-SL4 RRPRRR Chain

L1 = 2;
L2 = 3;

theta_P = [0, 0, 2, 0, pi/2, 0];

% Space Screws: Convert {ws} and {qs} into a 6x6 Slist
ws = {[0;0;1], [1;0;0], [0;0;0], [0;1;0], [1;0;0], [0;1;0]};
qs = {[0;0;0], [0;0;0], [0;1;0], [0;0;0], [0;L1;0], [0;0;0], [0;0;0]};

% Home position
M_P = [1 0 0 0; 
         0 1 0 L1+L2;
         0 0 1 0;
         0 0 0 1];

Slist_P = zeros(6, length(theta_P));
for i = 1:6
    wi = ws{i};
    if norm(wi) == 0
        Slist_P(:, i) = [wi; qs{i}];
    else
        vi = cross(-wi, qs{i}); 
        Slist_P(:, i) = [wi; vi];
    end
end


% Forward Kinematics (Space)
Ts_P = FK_space(M_P, Slist_P, theta_P);
disp('Ts_P')
disp(Ts_P)

%  derive Blist from Slist using the Adjoint of M
Blist_P = zeros(6, length(theta_P));
for i = 1:length(theta_P)
    Blist_P(:, i) = Adjoint(inv(M_P)) * Slist_P(:, i);
end

% Forward Kinematics (Body)

Tb_P = FK_body(M_P, Blist_P, theta_P);
disp('Tb_P')
disp(Tb_P)

% Check equivalence
error_P = norm(Ts_P - Tb_P);
fprintf('RRPRRR Space/Body Equivalence Error: %e\n', error_P);

