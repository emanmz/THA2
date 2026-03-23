% pt. g
function ellipsoid_plot_linear(J, T) % W7-L2 slide 5-6 and https://slama.dev/notes/robotics-1/slides/11_manipulability.pdf
% J: 6xn Jacobian matrix
% T: 4x4 current end-effector transformation matrix (for position)

    % Extract linear part (rows 4-6)
    Jv = J(4:6, :);
    
    % ellipsoid is defined by Jv * Jv'
    A = Jv * Jv';
    
    % SVD to get axes and magnitudes
    % U: eigenvectors (directions), S: eigenvalues (squares of semi-axes)
    [U, S, ~] = svd(A);
    
    % Semi-axis lengths are square roots of eigenvalues
    widths = sqrt(diag(S));
    
    % a unit sphere
    [X, Y, Z] = sphere(20);
    
    % Scale and rotate the sphere into an ellipsoid
    % We use the rotation matrix U and the scaling matrix sqrt(S)
    for i = 1:numel(X)
        point = U * diag(widths) * [X(i); Y(i); Z(i)];
        X(i) = point(1);
        Y(i) = point(2);
        Z(i) = point(3);
    end
    
    % Shift to the end-effector position
    p = T(1:3, 4);
    surf(X + p(1), Y + p(2), Z + p(3), 'FaceAlpha', 0.3, 'FaceColor', [1 0 1] , 'EdgeColor', 'none');
    hold on;
    
    % plot the principal axes
    for i = 1:3
        axis_vec = U(:,i) * widths(i);
        line([p(1), p(1)+axis_vec(1)], [p(2), p(2)+axis_vec(2)], [p(3), p(3)+axis_vec(3)], ...
            'Color', 'r', 'LineWidth', 2);
    end
    title('Linear Velocity Manipulability Ellipsoid');
    axis equal;
end