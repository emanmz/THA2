function [is_singular, mu, current_rank] = check_singularity(J, threshold)
% check_singularity: Determines if a configuration is singular
% J: 6xn Jacobian matrix (Space or Body)
% Outputs:
% is_singular: (true/false)
% mu: (sqrt(det(JJ')))
% current_rank: The numerical rank of the Jacobian

    if nargin < 2
        threshold = 1e-4; % Standard threshold for numerical noise
    end

    % A matrix (Task Space Ellipsoid)
    % For a 7-DOF robot, J is 6x7, so A is 6x6.
    A = get_A_matrix(J);

    % As mu -> 0, the ellipsoid flattens, indicating a singularity.
    mu = sqrt(det(A));

    % rank() uses SVD internally to see how many independent directions 
    % the robot can move in. For 6D task space, we want rank = 6.
    current_rank = rank(J);

    % A robot is singular if it loses a degree of freedom (rank < 6)
    % OR if the manipulability volume is effectively zero.
    if current_rank < 6 || mu < threshold
        is_singular = true;
    else
        is_singular = false;
    end
end