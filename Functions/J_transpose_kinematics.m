% pt. i
function [theta, success] = J_transpose_kinematics(Slist, M, Tsd, theta0)
% terative using the Jacobian Transpose method
    max_iter = 1000; % Transpose ??? iterations ?? i know it needs more idk how more
    alpha = 0.5;     % Step size (gain)
    theta = theta0;
    for i = 1:max_iter
        Tsb = FK_space(M, Slist, theta);
        [Vb, err] = calculate_twist_error(Tsb, Tsd);
        
        if ~err, break; end
        
        Js = J_space(Slist, theta);
        Vs = Adjoint(Tsb) * Vb; 
        
        % Update rule using Transpose instead of Pinv
        theta = theta + (alpha * Js' * Vs)';

    end
    success = ~err;
end