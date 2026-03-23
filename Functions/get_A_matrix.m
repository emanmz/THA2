% pt. g - n A = J*J'
function A = get_A_matrix(J)
    % A is symmetric positive semi-definite matrix JJ^T
    A = J * J';
end