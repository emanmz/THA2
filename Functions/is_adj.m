function adj = is_adj(M, S, B)
% S: 6x1 Space Screw
% B: 6x1 Body Screw
% M: 4x4 Home position

    R = M(1:3, 1:3);
    p = M(1:3, 4);
    p_sk = [0 -p(3) p(2); p(3) 0 -p(1); -p(2) p(1) 0];
    
    % Adjoint Matrix of M
    AdM = [R, zeros(3); p_sk*R, R];
    
    % Check if S = AdM * B
    S_calc = AdM * B;
    adj = norm(S - S_calc) < 1e-6;
end