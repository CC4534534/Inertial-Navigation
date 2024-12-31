% TRIAD
clc; clear;

tic;  % 開始計時
for i=1:100000
% 這裡是您希望計算執行時間的代碼區塊
yaw_pitch_roll = [30, 20, -10]*pi/180;
C_nb = Euler3212C(yaw_pitch_roll); %n frame to b frame;
s_N = [1 0 0]';
m_N = [0 0 -1]';
s_B_true = C_nb*s_N;
rng('default'); %set random
random = 0 + 0.01 * randn(3, 1); %mean=0, std=0.001
s_B = s_B_true+random;
% s_B = [0.8190 -0.5282 0.2242]';
s_B = s_B/norm(s_B);
m_B_true = C_nb * m_N;
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

barBN = barBT*barNT';
ANS = C2Euler123(barBN)/pi*180;
% barBB = barBN*C_nb';
% 
% BB = eye(3);

% ErrorPhiDeg = acos(0.5*(trace(barBB)-1))*180/pi;
% disp(ErrorPhiDeg);
end
disp(ANS);
elapsedTime = toc;  % 結束計時並獲取經過的時間
fprintf('程式執行時間為：%.6f 秒\n', elapsedTime);
