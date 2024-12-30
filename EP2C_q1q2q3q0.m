function Cm = EP2C_q1q2q3q0(q)
    % Convert a quaternion to a direction cosine matrix (DCM)
    % Input: q - A 1x4 vector representing the quaternion [q0, q1, q2, q3]
    % Output: Cm - A 3x3 direction cosine matrix (DCM)

    % Extract quaternion components
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q0 = q(4);

    % Initialize the direction cosine matrix
    Cm = eye(3);

    % Compute the elements of the DCM
    Cm(1, 1) = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    Cm(1, 2) = 2 * (q1*q2 + q0*q3);
    Cm(1, 3) = 2 * (q1*q3 - q0*q2);
    Cm(2, 1) = 2 * (q1*q2 - q0*q3);
    Cm(2, 2) = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    Cm(2, 3) = 2 * (q2*q3 + q0*q1);
    Cm(3, 1) = 2 * (q1*q3 + q0*q2);
    Cm(3, 2) = 2 * (q2*q3 - q0*q1);
    Cm(3, 3) = q0*q0 - q1*q1 - q2*q2 + q3*q3;
end