% Davenport’s q-method
clc; clear;
yaw_pitch_roll = [30, 20, -10]*pi/180;

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

W = [s_B m_B]; % the corresponding observation vectors.
V = [s_N m_N]; % the nonparallel reference vectors
a = [0.5 0.5];
B = zeros(3);
l_opt = 0;
for i=1:2
    B = B + a(i) * W(:,i) * V(:,i)';
    l_opt = l_opt + 1;
end

sigma = trace(B);
S = B+B';
Z = [0;0;0];
Z = [B(2,3) - B(3,2) ; ...
     B(3,1) - B(1,3) ; ...
     B(1,2) - B(2,1) ];
% for i=1:2
%     Z = Z + a(i) * cross(W(:,i),V(:,i));
% end
%%用matlab逆乘求矩陣
% a = (l_opt + sigma)*eye(3);
% p = (a-S)\Z;    
% 
% b = 1/sqrt(1+p'*p);
% q = b*[1; p];

K = [ S-sigma*eye(3) Z ; ...
     Z' sigma];

kappa = trace(adjoint(S));
delta = det(S);
a = sigma^2 - kappa;
b = sigma^2 + Z' * Z;
c = delta + Z' * S * Z;
d = Z' * S^2 * Z;
% 牛頓-拉夫森方法初始化
lambda_t = 1;  % 初始值（可以根據實際情況調整）
tolerance = 1e-6;  % 收斂容忍度
max_iter = 100;  % 最大迭代次數
iter = 0;  % 迭代計數器
old_lambda = 1.0;
% 牛頓-拉夫森方法迭代
while iter < max_iter
    % 計算 f(lambda) 和 f'(lambda)
    f_lambda = lambda_t^4 - (a + b) * lambda_t^2 - c * lambda_t + (a * b + c * sigma - d);
    f_prime_lambda = 4 * lambda_t^3 - 2 * (a + b) * lambda_t - c;

    % 更新 lambda
    lambda_t1 = lambda_t - f_lambda / f_prime_lambda;

    % 檢查收斂
    if abs(lambda_t1 - lambda_t) < tolerance
        break;
    end

    lambda_t = lambda_t1  % 更新 lambda
    iter = iter + 1;
end

alpha = lambda_t^2 - sigma^2 +kappa;
gamma = (lambda_t + sigma) * alpha - delta;
beta = lambda_t - sigma;
X = (alpha * eye(3) + beta * S + S^2)*Z;

beta_max = [X; gamma]/sqrt(gamma^2 + norm(X)^2)

barBN = EP2C_q1q2q3q0(beta_max);
barBB = barBN*BN';
ErrorPhi = acos(0.5*(trace(barBB)-1))*180/pi

