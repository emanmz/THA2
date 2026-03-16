function T = FK_body(M, Bn, theta)  % W6-2 slide 6 (sort of 2-6) & W4-L1 slide 12 (fixed frame vs body frame)
% M: 4x4 home position configuration matrix
% Blist: 6xn matrix of body screw axes [wi; vi]
% theta: 1xn vector of joint variables
% Output:
% T: 4x4 configuration matrix of the end-effector

    T = M;
    for i = 1:length(theta)
        % multiply the home matrix by the exponentials from left to right
        T = T * screw_to_exp(Bn(:,i), theta(i));
    end
end