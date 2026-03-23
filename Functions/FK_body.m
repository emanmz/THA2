function T = FK_body(M, Bn, theta)  % W6-2 slide 6 (sort of 2-6) & W4-L1 slide 12 (fixed frame vs body frame)
% M: 4x4 home position configuration matrix
% Blist: 6xn matrix of body screw axes [wi; vi]
% theta: 1xn vector of joint variables
% Output:
% T: 4x4 configuration matrix of the end-effector

    T = M;
    for i = 1:length(theta)
        % pre multiply the home matrix by the exponentials from left to right
        T = T * screw_to_exp(Bn(:,i), theta(i));
    end

     % Create Plot
    figure; hold on; grid on; axis equal;
    view(3); xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Body Form Forward Kinematics');

    % Plot Body Frame {b}
    drawFrame(eye(4), 0.5, '{b}');
    
    % Plot Space Frame {s}
    drawFrame(T, 0.5, '{s}');


end

% helper func tion to draw frames
function drawFrame(T, scale, label)
    % Draws a coordinate frame at configuration T
    pos = T(1:3, 4);
    R = T(1:3, 1:3);
    colors = ['r', 'g', 'b']; % X=Red, Y=Green, Z=Blue
    for i = 1:3
        quiver3(pos(1), pos(2), pos(3), R(1,i), R(2,i), R(3,i), scale, ...
                'Color', colors(i), 'LineWidth', 2);
    end
    text(pos(1), pos(2), pos(3), label, 'FontSize', 10, 'FontWeight', 'bold');
end


