function T = FK_space(M, Sn, theta) % W6-2 slide 6 (sort of 2-6) & W4-L1 slide 12 (fixed frame vs body frame)
% M: 4x4 home position configuration matrix
% Slist: 6xn matrix of spatial screw axes [wi; vi]
% theta: 1xn vector of joint variables
% Output:
% T: 4x4 configuration matrix of the end-effector

    T = eye(4);
    for i = 1:length(theta)
        % multiply exponentials from left to right
        T = T * screw_to_exp(Sn(:,i), theta(i));
    end
    
    % Post-multiply by the home position
    T = T * M;
end