function V = MatrixLog6(T) % this is from tha 1 
    % Extracts the twist V [w; v] from a transformation matrix T
    R = T(1:3, 1:3);
    p = T(1:3, 4);
    
    if abs(trace(R) - 3) < 1e-6 % No rotation
        w = [0;0;0];
        v = p;
    else
        theta = acos((trace(R) - 1) / 2);
        w_sk = (1 / (2 * sin(theta))) * (R - R');
        w = [w_sk(3,2); w_sk(1,3); w_sk(2,1)];
        
        % Calculate translation part of the twist
        G_inv = eye(3)/theta - 0.5*w_sk + (1/theta - 0.5*cot(theta/2))*w_sk^2;
        v = G_inv * p;
        V = [w; v] * theta;
        return;
    end
    V = [w; v];
end