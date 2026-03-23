% pt. g
function vol = J_ellipsoid_volume(J) % W7-L2 slide 7
% Volume: Proportional to sqrt(det(A))
% Slide: Volume is proportional to the product of the lengths of the semi-axes
A = get_A_matrix(J);

% slide's direct definition:
vol = sqrt(det(A));
end