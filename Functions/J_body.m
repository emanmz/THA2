% EE to base
function Jb = J_body(Bn, theta) % W7-L1-SL6

% n = number of joints
% Bn = 6xn matrix of screw axis in body frame
% theta = 1xn matrix of joint angles

n = length(theta);
Jb = zeros(6,n);
T = eye(4);
% last is last screw
Jb(:,n) = Bn(:,n);

for i = n-1:-1:1

    % calculate transformation for ith joint
    Ti_next = screw_to_exp(-Bn(:, i+1), theta(i+1));

    % cumulative transformation
    T = T * Ti_next;

    % trans the current body screw axis
    Jb(:, i) = Adjoint(T) * Bn(:, i);
end
end
