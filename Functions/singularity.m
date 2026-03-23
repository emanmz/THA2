function SingAnalysis = singularity(Sn, thetaSym)
n = length(thetaSym);
Js = sym(zeros(6, n));
T = eye(4);

% Build Symbolic Jacobian
Js(:, 1) = Sn(:, 1);
for i = 2:n
    T = T * screw_to_exp(Sn(:, i-1), thetaSym(i-1));
    Js(:, i) = Adjoint(T) * Sn(:, i);
end

if n < 6
    % For under-actuated robots check the 4x4 Gramian
    manip_term = simplify(det(Js' * Js));
else
    n = length(thetaSym);
    Js = sym(zeros(6, n));
    T = eye(4);

    % Build Symbolic Jacobian
    Js(:, 1) = Sn(:, 1);
    for i = 2:n
        T = T * screw_to_exp(Sn(:, i-1), thetaSym(i-1));
        Js(:, i) = Adjoint(T) * Sn(:, i);
    end

    % SHORTCUT
    % assume q1, q2, q3, q5, q6, q7 are 'arbitrary' (non-zero constants)
    % to see what q4 needs to be to break the system.


    % core trigonometric term for the manipulability
    manip_term = simplify(det(Js(:, 1:6)));

    fprintf('Analytical term found. Solving for Elbow Singularity (theta4)...\n');

    %  solve for theta4, which is the geometric center of the arm
    SingAnalysis = solve(manip_term == 0, thetaSym(4));
end

fprintf('Analytical term found. Solving for roots...\n');

% Solve for the primary joint (theta2 for SCARA, theta4 for Franka)
if n == 4
    SingAnalysis = solve(manip_term == 0, thetaSym(2));
else
    SingAnalysis = solve(manip_term == 0, thetaSym(4));
end
end