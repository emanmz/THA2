% PA pt. C 
function T = FK_body(M, Bn, theta)  % W6-2 slide 6 (sort of 2-6) & W4-L1 slide 12 (fixed frame vs body frame)
% M: 4x4 home position configuration matrix
% Blist: 6xn matrix of body screw axes [wi; vi]
% theta: 1xn vector of joint variables
% Output:
% T: 4x4 configuration matrix of the end-effector

    T = M;

     % Create Plot
    figure; hold on; grid on; axis equal;
    view(3); xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Body Form Forward Kinematics');

    % Plot Body Frame {b}
    drawFrame(eye(4), 0.5, '{b}');
    for i = 1:length(theta)
        % pre multiply the home matrix by the exponentials from left to right
        T = T * screw_to_exp(Bn(:,i), theta(i));
        drawScrewAxis(Bn(:,i), 'r');
    end
    % Plot Space Frame {s}
    drawFrame(T, 0.5, '{s}');


end