% TRIAD
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
C_nb = Euler3212C(yaw_pitch_roll); %n frame to b frame;
s_N = [1 0 0]';
m_N = [0 0 -1]';
s_B_true = C_nb*s_N
rng('default'); %set random
random = 0 + 0.01 * randn(3, 1); %mean=0, std=0.001
s_B = s_B_true+random;
% s_B = [0.8190 -0.5282 0.2242]';
s_B = s_B/norm(s_B);
m_B_true = C_nb * m_N
m_B = m_B_true+random;
% m_B = [0.3138 0.1548 -0.9362]';
m_B = m_B/norm(m_B);
% set T frame to solve C n frame to b frame
t1_B = s_B;
t2_B = cross(s_B,m_B)/norm(cross(s_B,m_B));
t3_B = cross(t1_B,t2_B);
barBT = [t1_B t2_B t3_B];

t1_N = s_N;
t2_N = cross(s_N,m_N)/norm(cross(s_N,m_N));
t3_N = cross(t1_N,t2_N);
barNT = [t1_N t2_N t3_N];

barBN = barBT*barNT'
disp(C_nb);

barBB = barBN*C_nb';

BB = eye(3);

ErrorPhiDeg = acos(0.5*(trace(barBB)-1))*180/pi;
disp(ErrorPhiDeg);



