function cond_num = J_condition(J)
% high values indicate closeness to a singularity.

    s = svd(J);
    % Standard definition: ratio of max to min singular values
    cond_num = max(s) / min(s);
    
    % If min(s) is 0, MATLAB returns Inf, which is correct for
% singularities 
end