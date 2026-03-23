% pt. g
function iso = J_isotropy(J) % W7-L2 slide 7
% Isotropy: Ratio of min eigenvalue to max eigenvalue of A
% Slide: sqrt(L_min / L_max)
A = get_A_matrix(J);
ev = eig(A);
iso = sqrt(min(ev) / max(ev));
end