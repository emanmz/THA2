function T = FK_body(M, Bn, theta)
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