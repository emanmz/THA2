function T = FK_body_no_plot(M, Bn, theta)
    % T = M * e^[B1]th1 * e^[B2]th2 * ... * e^[Bn]thn
    T = M;
    for i = 1:length(theta)
        T = T * screw_to_exp(Bn(:,i), theta(i));
    end
end