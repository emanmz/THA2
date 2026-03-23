% pt. f
function SingAnalysis = singularity(Sn, thetaSym)

% thetaSym = 1xn matrix of symbolic joint angles
% Sn = 6xn matrix of screw axis in space frame

% calcualte symbolic space jacobian
n = length(thetaSym);
Js = sym(zeros(6,n));
Js(:,1) = Sn(:,1);
T = eye(4);
% Build symbolic Jacobian
    Js(:, 1) = Sn(:, 1);
    for i = 2:n
        T = T * screw_to_exp(Sn(:, i-1), thetaSym(i-1));
        Js(:, i) = Adjoint(T) * Sn(:, i);
    end

%use det to make it computationally lighter
manipM = simplify(det(Js * Js')); % manipuability measure (W8-L1 slide 10)

SingAnalysis = solve(manipM == 0, thetaSym);

end


