function drawScrewAxis(S, col)
    % Simple visualization of screw axis direction
    w = S(1:3);
    v = S(4:6);
    if norm(w) > 1e-6
        % Revolute: Draw a line through a point on the axis
        q = cross(w, v) / norm(w)^2;
        quiver3(q(1), q(2), q(3), w(1), w(2), w(3), 2, 'Color', col, 'LineStyle', '--');
    else
        % Prismatic: Draw the direction of translation
        quiver3(0, 0, 0, v(1), v(2), v(3), 2, 'Color', 'b', 'LineStyle', ':');
    end
end