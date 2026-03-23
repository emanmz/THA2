function [theta, success] = redundancy_resolution(Slist, M, Tsd, theta0)
% maximizes manipulability in the null space
    max_iter = 100;
    k = 0.1; % Gain for the secondary objective
    theta = theta0;
    
    for i = 1:max_iter
        Tsb = FK_space_no_plot(M, Slist, theta);
        [Vb, err] = calculate_twist_error(Tsb, Tsd);
        if ~err, break; end
        
        Js = J_space(Slist, theta);
        Vs = Adjoint(Tsb) * Vb;
        
        % Primary Task: Pseudoinverse
        J_pinv = pinv(Js);
        
        % Secondary Task: Gradient of Manipulability w = sqrt(det(J*J'))
        % define a secondary objective function of the joint variables that maximizing this manipulability measure and exploits redundancy to move away from singularities
        eps = 1e-4;
        grad_w = zeros(length(theta), 1);
        for j = 1:length(theta)
            th_plus = theta; th_plus(j) = th_plus(j) + eps;
            w_plus = J_ellipsoid_volume(J_space(Slist, th_plus));
            
            th_minus = theta; th_minus(j) = th_minus(j) - eps;
            w_minus = J_ellipsoid_volume(J_space(Slist, th_minus));
            
            grad_w(j) = (w_plus - w_minus) / (2 * eps);
        end
        
        % Combined: Progress to goal + Null Space Projection
        null_projection = (eye(size(Js,2)) - J_pinv * Js) * (k * grad_w);
        theta = theta + (J_pinv * Vs + null_projection)';
    end
    success = ~err;
end