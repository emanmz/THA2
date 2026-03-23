function [theta, success] = J_inverse_kinematics(Slist, M, Tsd, theta0, eomg, ev)
% Slist: 6xn matrix of screw axes
% M: Home configuration matrix
% Tsd: Desired end-effector configuration (4x4)
% theta0: Initial guess for joint angles (1xn)
% eomg: Error tolerance for orientation (rad)
% ev: Error tolerance for translation (m)

    max_iter = 50; % Limit iterations to prevent infinite loops
    theta = theta0;
    i = 0;
    
    % calc current FK and the twist required to move toward target
    Tsb = FK_space(M, Slist, theta);
    [Vb, err] = calculate_twist_error(Tsb, Tsd);
    
    % until the error is below tolerance or max its 
    while err && i < max_iter
        % Calculate Space Jacobian at current theta
        Js = J_space(Slist, theta);
        
        %  Body Twist (Vb) to Space Twist (Vs) 
        Vs = Adjoint(Tsb) * Vb; 
        
        % Newton-Raphson Update: theta_new = theta_old + pinv(J) * Vs
        % use pinv (pseudoinverse) because the FR3 is redundant (6x7)
        % (franka is redunt
        theta = theta + (pinv(Js) * Vs)';
        
        % Update current state
        i = i + 1;
        Tsb = FK_space(M, Slist, theta);
        [Vb, err] = calculate_twist_error(Tsb, Tsd);
    end
    
    success = ~err;
    if ~success
        warning('IK did not converge within % d iterations', max_iter); % maybe change this for no max its but likeeeee i think we neeed
    end
end