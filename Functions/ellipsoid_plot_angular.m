% pt. g
function ellipsoid_plot_angular(J, T) % W7-L2 slide 5-6 and https://slama.dev/notes/robotics-1/slides/11_manipulability.pdf
% J: 6xn Jacobian matrix
% T: 4x4 current end-effector transformation matrix

    % Extract angular part (rows 1-3)
    Jw = J(1:3, :);
    
    % A = Jw * Jw'
    A = Jw * Jw';
    
    [U, S, ~] = svd(A);
    widths = sqrt(diag(S));
    
    [X, Y, Z] = sphere(20);
    for i = 1:numel(X)
        point = U * diag(widths) * [X(i); Y(i); Z(i)];
        X(i) = point(1);
        Y(i) = point(2);
        Z(i) = point(3);
    end
    
    % Plotting at end-effector position
    p = T(1:3, 4);
    surf(X + p(1), Y + p(2), Z + p(3), 'FaceAlpha', 0.3, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
    hold on;
    
    % Plot axes
    for i = 1:3
        axis_vec = U(:,i) * widths(i);
        line([p(1), p(1)+axis_vec(1)], [p(2), p(2)+axis_vec(2)], [p(3), p(3)+axis_vec(3)], ...
            'Color', 'b', 'LineWidth', 2);
    end
    title('Angular Velocity Manipulability Ellipsoid');
    axis equal;
end