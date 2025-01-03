% Davenport’s q-method
clc; clear;

tic;  % 開始計時
for i=1:100000
yaw_pitch_roll = [30, 20, -10]*pi/180;
BN = Euler3212C(yaw_pitch_roll); %n frame to b frame;
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

B = w1*v1_B*v1_N' + w2*v2_B*v2_N';
sigma = trace(B);
S = B+B';
Z = [B(2,3) - B(3,2) ; ...
     B(3,1) - B(1,3) ; ...
     B(1,2) - B(2,1) ];

K = [ sigma Z' ; ...
     Z (S-sigma*eye(3))];
% % 使用 Lagrange Multiplier 找最大特徵值 錯的
% q = rand(4, 1); % 初始化四元數（假設隨機初始化）
% q = q / norm(q); % 確保單位化
% tol = 1e-8; % 誤差容忍度
% max_iter = 1000; % 最大迭代次數
% lambda = 1; % 初始化特徵值
% 
% for iter = 1:max_iter
%     q_new = K * q; % 計算 K*q
%     lambda_new = q_new' * q; % 更新拉格朗日乘數（特徵值）
%     q_new = q_new / norm(q_new); % 單位化
% 
%     % 檢查收斂
%     if norm(q_new - q) < tol && abs(lambda_new - lambda) < tol
%         break;
%     end
%     q = q_new;
%     lambda = lambda_new;
% end

%% 利用eig跟mad找最大特徵值
% [eigvec,eigval]=eigs(K)
% % 自動找到最大特徵值的位置
% [~, maxIdx] = max(diag(eigval)); % 找到最大特徵值的位置
% beta_max = eigvec(:, maxIdx);    % 提取對應的特徵向量
% 目標函數 q'*K*q，約制條件 norm(q) = 1
fun = @(q) -(q' * K * q); % 負號用於找最大
nlcons = @(q)deal(norm(q)^2 - 1, []); % 非線性約制：q 的模設為 1

% 給定初始值
q0 = rand(4, 1); % 隨機初始化
q0 = q0 / norm(q0); % 單位化初始值

% 邊界值設置(無邊界)
lb = []; ub = [];

% fmincon 優化選項
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

% 應用 fmincon
[q_opt, fval, exitflag, output, lambda] = fmincon(fun, q0, [], [], [], [], lb, ub, nlcons, options);

% 四元數
beta_max = q_opt; % 最佳解就是最大特徵值對應之特徵向量
norm(beta_max);
barBN = EP2C(beta_max);
ANS = C2Euler123(barBN)/pi*180;
% barBB = barBN*BN'
% ErrorPhi = acos(0.5*(trace(barBB)-1))*180/pi
end
disp(ANS);
elapsedTime = toc;  % 結束計時並獲取經過的時間

fprintf('程式執行時間為：%.6f 秒\n', elapsedTime);


