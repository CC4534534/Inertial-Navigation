% Davenport’s q-method
clc; clear;
yaw_pitch_roll = [30, 20, -10]*pi/180;
function Euler3212C = Euler3212C(q)
    % Calculate the rotation matrix using a 3-2-1 Euler angle sequence
    % Input: q - A vector of 3 Euler angles [q1, q2, q3]
    % Output: C - The 3x3 rotation matrix

    % Precompute trigonometric values
    st1 = sin(q(1));
    ct1 = cos(q(1));
    st2 = sin(q(2));
    ct2 = cos(q(2));
    st3 = sin(q(3));
    ct3 = cos(q(3));

    % Construct the rotation matrix
    Euler3212C = [
        ct2 * ct1,               ct2 * st1,               -st2;
        st3 * st2 * ct1 - ct3 * st1, st3 * st2 * st1 + ct3 * ct1, st3 * ct2;
        ct3 * st2 * ct1 + st3 * st1, ct3 * st2 * st1 - st3 * ct1, ct3 * ct2
    ];
end
function Cm = EP2C(q)
    % Convert a quaternion to a direction cosine matrix (DCM)
    % Input: q - A 1x4 vector representing the quaternion [q0, q1, q2, q3]
    % Output: Cm - A 3x3 direction cosine matrix (DCM)

    % Extract quaternion components
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

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


BN = Euler3212C(yaw_pitch_roll) %n frame to b frame;
s_N = [1 0 0]';
m_N = [0 0 -1]';
s_B_true = BN*s_N;
rng('default'); %set random
random = 0 + 0.01 * randn(3, 1); %mean=0, std=0.001
s_B = s_B_true+random;
% s_B = [0.8190 -0.5282 0.2242]';
s_B = s_B/norm(s_B);
m_B_true = BN * m_N;
m_B = m_B_true+random;
% m_B = [0.3138 0.1584 -0.9362]';
m_B = m_B/norm(m_B);

v1_N = s_N;
v2_N = m_N;
v1_B = s_B;
v2_B = m_B;
w1=1;
w2=1;

B = w1*v1_B*v1_N' + w2*v2_B*v2_N'
sigma = trace(B)
S = B+B'
Z = [B(2,3) - B(3,2) ; ...
     B(3,1) - B(1,3) ; ...
     B(1,2) - B(2,1) ]

K = [ sigma Z' ; ...
     Z (S-sigma*eye(3))]

[eigvec,eigval]=eigs(K)
beta_max = eigvec(:,1)
norm(beta_max)
barBN = EP2C(beta_max)
barBB = barBN*BN'
ErrorPhi = acos(0.5*(trace(barBB)-1))*180/pi

