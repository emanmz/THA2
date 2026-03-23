function plot_2link(Slist, M, theta, Tsd)
    L1 = 1; 
    p0 = [0; 0; 0];
    T1 = screw_to_exp(Slist(:,1), theta(1));
    p1 = T1 * [L1; 0; 0; 1];
    Te = FK_space(M, Slist, theta);
    p2 = Te(1:3, 4);
    
    plot([p0(1), p1(1), p2(1)], [p0(2), p1(2), p2(2)], 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    hold on;
    plot(Tsd(1,4), Tsd(2,4), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Target
    grid on; axis equal; xlim([-0.5 2.5]); ylim([-0.5 2.5]);
    title('IK Iterative Convergence');
    hold off;
    drawnow;
end