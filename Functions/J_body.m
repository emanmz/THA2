
function Jb = J_body(Bn, theta) % W7-L1-SL6

% n = number of joints
% Bn = 6xn matrix of screw axis in body frame
% theta = 1xn matrix of joint angles

Adj = eye(6);
n = length(theta);
Jb = zeros(6,n);
Jb(:,n) = Bn(:,n);

for i = n-1:-1:1

    % calculate transformation for ith joint
    Ti_prev = screw_to_exp(-Bn(:,i+1),theta(i+1));
    
    % change to adjoint
    Adi_prev = Adjoint(Ti_prev);

    % multiply it to previous adjoints
    Adj = Adj*Adi_prev; 

    % Jacobian column for i+1 column 
     Jb(:,i) = Adj*Bn(:,i);
end 
