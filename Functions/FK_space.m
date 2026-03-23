% PA pt. B
function T = FK_space(M, Sn, theta) % W6-2 slide 6 (sort of 2-6) & W4-L1 slide 12 (fixed frame vs body frame)
% M: 4x4 home position configuration matrix
% Slist: 6xn matrix of spatial screw axes [wi; vi]
% theta: 1xn vector of joint variables
% Output:
% T: 4x4 configuration matrix of the end-effector

    T = eye(4);

    % Create Plot
    figure; hold on; grid on; axis equal;
    view(3); xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Space Form Forward Kinematics');

    % Plot Space Frame {s}
    drawFrame(eye(4), 0.5, '{s}');
    
    for i = 1:length(theta)
       
        % multiply exponentials from left to right
        T = T * screw_to_exp(Sn(:,i), theta(i));
        drawScrewAxis(Sn(:,i), 'r');
    end
    
    % Post-multiply by the home position
    T = T * M;

    % Plot Body Frame {b}
    drawFrame(T, 0.5, '{b}');
    
end