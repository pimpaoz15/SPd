function R = rotation_transform(q)
    % rotation_transform - Convert quaternion to rotation matrix or apply specific rotation
    % Input:
    %   q - quaternion (4x1 vector)
    % Output:
    %   R - rotation matrix (3x3 matrix)
    
    % Ensure quaternion is normalized
    q = q / norm(q);
    
    % Quaternion components
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    % Compute the rotation matrix
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end