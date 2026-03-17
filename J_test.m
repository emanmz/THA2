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
