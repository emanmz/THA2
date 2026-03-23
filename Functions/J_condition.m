% pt. g
function cond_num = J_condition(J) % W7-L2 slide 7
% Condition: Ratio of max eigenvalue to min eigenvalue of A
% Slide: L_max / L_min
A = get_A_matrix(J);
ev = eig(A);
cond_num = max(ev) / min(ev);
end