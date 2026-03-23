% pt. g
function iso = J_isotropy(J) % W7-L2 slide 7
%ratio of the smallest to largest singular value
% J can be the full Jacobian or a subpart (Jv or Jw)

    s = svd(J);
    % Ratio of minimum singular value to maximum
    iso = sqrt(min(s)) / sqrt(max(s));
end