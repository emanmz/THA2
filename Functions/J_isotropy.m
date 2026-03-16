function iso = J_isotropy(J)
%ratio of the smallest to largest singular value
% J can be the full Jacobian or a subpart (Jv or Jw)

    s = svd(J);
    % Ratio of minimum singular value to maximum
    iso = min(s) / max(s);
end