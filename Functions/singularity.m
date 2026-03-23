function SingAnalysis = singularity(Sn, thetaSym)
% Sn = 6xn matrix of screw axes
% thetaSym = 1xn symbolic joint variables
    n = length(thetaSym);
    Js = sym(zeros(6, n));
    T = eye(4);
    
    % Build symbolic Space Jacobian
    Js(:, 1) = Sn(:, 1);
    for i = 2:n
        T = T * screw_to_exp(Sn(:, i-1), thetaSym(i-1));
        Js(:, i) = Adjoint(T) * Sn(:, i);
    end
    
    % correct manipulability measure based on DOF
    if n < 6
        % J'J for robots with fewer than 6 joints
        manipM = simplify(det(Js' * Js));
    else
        % JJ' for robots with 6 or more joints (like the Emika)
        manipM = simplify(det(Js * Js'));
    end
    
    % for singularities (where manipulability is zero)
    SingAnalysis = solve(manipM == 0, thetaSym);
end