% helper func tion to draw frames
function drawFrame(T, scale, label)
    % Draws a coordinate frame at configuration T
    pos = T(1:3, 4);
    R = T(1:3, 1:3);
    colors = ['r', 'g', 'b']; % X=Red, Y=Green, Z=Blue
    for i = 1:3
        quiver3(pos(1), pos(2), pos(3), R(1,i)*scale, R(2,i)*scale, R(3,i)*scale, ...
                'Color', colors(i), 'LineWidth', 2);
    end
    text(pos(1), pos(2), pos(3), label, 'FontSize', 10, 'FontWeight', 'bold');
end
