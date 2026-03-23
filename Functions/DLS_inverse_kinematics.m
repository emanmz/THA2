function [theta, success] = DLS_inverse_kinematics(Slist, M, Tsd, theta0)
% robot through singularities using damping
    max_iter = 100;
    lambda_max = 0.1; % Max damping
    epsilon = 0.05;   % Region of influence near singularity
    theta = theta0;
    
    for i = 1:max_iter
        Tsb = FK_space_no_plot(M, Slist, theta);
        [Vb, err] = calculate_twist_error(Tsb, Tsd);
        if ~err, break; end
        
        Js = J_space(Slist, theta);
        Vs = Adjoint(Tsb) * Vb;
        
        % smallest singular value to detect singularity
        s = svd(Js);
        s_min = min(s);
        
        % Adjust damping: higher damping when closer to singularity (fast
        % away slow closer 
        if s_min < epsilon
            lambda = (1 - (s_min/epsilon)^2) * lambda_max^2;
            fprintf('Iteration %d: Near Singularity (s_min=%.4f). Damping engaged.\n', i, s_min);
        else
            lambda = 0;
        end
        
        % Damped Inverse: J_dls = J' * inv(J*J' + lambda*I)
        J_dls = Js' / (Js * Js' + lambda * eye(6));
        theta = theta + (J_dls * Vs)';
    end
    success = ~err;
end