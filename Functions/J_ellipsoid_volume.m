function vol = J_ellipsoid_volume(J)
% Calculates a value proportional to the volume
% J: 3xn sub-Jacobian (either Jv or Jw)

    % The volume of an ellipsoid with semi-axes a, b, c is (4/3)*pi*a*b*c
    s = svd(J);
    
    % product of the singular values
    vol = (4/3) * pi * prod(s);
end