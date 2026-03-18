
function SingAnalysis = singularity(Sn, thetaSym)

% thetaSym = 1xn matrix of symbolic joint angles
% Sn = 6xn matrix of screw axis in space frame

% calcualte symbolic space jacobian
Adj = eye(6);
n = length(thetaSym);
Js = sym(zeros(6,n));
Js(:,1) = Sn(:,1);

for i = 2:n

    % calculate transformation for ith joint
    Ti_prev = screw_to_exp(Sn(:,i-1),thetaSym(i-1));
    
    % change to adjoint
    Adi_prev = Adjoint(Ti_prev);

    % multiply it to previous adjoints
    Adj = Adj*Adi_prev; 

    % Jacobian column for i+1 column 
     Js(:, i) = Adj*Sn(:,i);
end 

manipM = sqrt(det(Js*Js')); % manipuability measure (W8-L1 slide 10)

SingAnalysis = solve(manipM == 0, thetaSym);

end


