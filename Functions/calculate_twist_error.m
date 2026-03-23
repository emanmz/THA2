function [Vb, err] = calculate_twist_error(Tsb, Tsd) % helper func for j inverse kinematics 
    % Calculates the twist matrix log and checks if error is within tolerance
    % Tdiff = inv(Tsb) * Tsd (Transformation from current to desired)
    Tdiff = inv(Tsb) * Tsd;
    
    % Se3ToVec: This is essentially the matrix log of the transformation
    % It returns a 6x1 twist [omega; v]
    Vb = MatrixLog6(Tdiff); 
    
    % Define tolerances
    eomg = 1e-3; 
    ev = 1e-3;
    
    % Check if the norm of angular and linear velocity is small enough
    err = norm(Vb(1:3)) > eomg || norm(Vb(4:6)) > ev;
end