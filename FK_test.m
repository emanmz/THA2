addpath("Functions");

%% Barrett WAM arm W6-L2-SL 3-4 

L1 = .550; % m
L2 = .300;
L3 = .060;
W1 = .045;

theta_wam = [0, pi/4, 0, -pi/4, 0, -pi/2, 0];
% Body Screws: Convert your {wb} and {qb} into a 6x7 Blist
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
    vi = cross(wi, -qb{i}); % Body screw v calculation
    Blist_wam(:, i) = [wi; vi];
end


% Forward Kinematics (Body)
Tb_WAM = FK_body(M_wam, Blist_wam, theta_wam);

%  derive Slist from Blist using the Adjoint of M
Slist_wam = zeros(6, 7);
for i = 1:7
    Slist_wam(:, i) = Adjoint(M_wam) * Blist_wam(:, i);
end
% also check if space and body are equivalent using adjoint. W6-L1-SL15 i

% Forward Kinematics (Space)
Ts_WAM = FK_space(M_wam, Slist_wam, theta_wam);

% Check equivalence
error_wam = norm(Ts_WAM - Tb_WAM);
fprintf('WAM Space/Body Equivalence Error: %e\n', error_wam);
%% FK Test Cases

% W6-L1-SL11 3R Example

L1 = 2;
L2 = 3;
M_3r = [0 0 1 L1; 0 1 0 0; -1 0 0 -L2; 0 0 0 1];
w_3r = {[0;0;1], [0;-1;0], [1;0;0]};
q_3r = {[0;0;0], [0;L1;0], [0;0;-L2]};
theta_3r = [pi/2, pi, 3*pi/2];

% Convert to Slist
Slist_3r = zeros(6, 3);
for i = 1:3
    wi = w_3r{i};
    vi = cross(-wi, q_3r{i}); % Space screw v calculation
    Slist_3r(:, i) = [wi; vi];
end

Ts_3r = FK_space(M_3r, Slist_3r, theta_3r);

% Derive Blist for 3R to check the other way
Blist_3r = zeros(6, 3);
for i = 1:3
    Blist_3r(:, i) = Adjoint(inv(M_3r)) * Slist_3r(:, i);
end
Tb_3r = FK_body(M_3r, Blist_3r, theta_3r);

fprintf('3R Space/Body Equivalence Error: %e\n', norm(Ts_3r - Tb_3r));

function Ad = Adjoint(T)
    R = T(1:3, 1:3);
    p = T(1:3, 4);
    p_sk = [0 -p(3) p(2); p(3) 0 -p(1); -p(2) p(1) 0];
    Ad = [R, zeros(3); p_sk*R, R];
end