% base 2 joint
function Js = J_space(Sn, theta) % W7-L1-SL6

% n = number of joints
% Sn = 6xn matrix of screw axis in space frame
% theta = 1xn matrix of joint angles

n = length(theta);
Js = zeros(6,n);
% first column is screw
Js(:,1) = Sn(:,1);

for i = 2:n

    % calculate transformation for ith joint
    Ti_prev = screw_to_exp(Sn(:,i-1),theta(i-1));

    T = T * Ti_prev;

    % Column i is the screw axis i transformed by the Adjoint of T

    Js(:, i) = Adjoint(T) * Sn(:, i);
end
end

% function Js = J_space(Sn, theta)
% 
% % Sn = matrix of screw axis in home position
% % theta = matrix of joint angles
% 
% Adj = eye(6);
% n = length(theta);
% Js = zeros(6,n);
% Js(:,1) = Sn(:,1);
% 
% for i = 2:n
% 
%     % calculate transformation for ith joint
%     Ti_prev = screw_to_exp(Sn(:,i-1),theta(i-1));
% 
%     % change to adjoint
%     Adi_prev = Adjoint(Ti_prev);
% 
%     % multiply it to previous adjoints
%     Adj = Adj*Adi_prev; 
% 
%     % Jacobian column for i+1 column 
%      Js(:, i) = Adj*Sn(:,i);
% end 
