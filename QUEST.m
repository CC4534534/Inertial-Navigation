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
m_N = [0 0 1]';
s_B_true = BN*s_N;
rng('default'); %set random
random = 0 + 0.001 * randn(3, 1); %mean=0, std=0.001
% s_B = s_B_true+random;
s_B = [0.8190 -0.5282 0.2242]';
s_B = s_B/norm(s_B);
m_B_true = BN * m_N;
% m_B = m_B_true+random;
m_B = [-0.3138 -0.1584 0.9362]';
m_B = m_B/norm(m_B);

v1_N = s_N;
v2_N = m_N;
v1_B = s_B;
v2_B = m_B;
a1=1;
a2=1;
a = [a1, a2];
B = zeros(3);
W = [s_B, m_B];  % 假設這是 3x2 的矩陣
V = [s_N, m_N];  % 假設這是 3x2 的矩陣
for i = 1:2
    B = B + a(i) * W(:, i) * V(:, i)'; 
end
sigma = trace(B);
S = B+B';

% 初始化 Z 向量
Z = [0, 0, 0]';  % Z 是 3x1 向量

% 計算 Z = sum(a_i * (W_i x V_i))
for i = 1:2
    Z = Z + a(i) * cross(W(:, i), V(:, i));  % cross 函數計算向量的叉積
end
K = [(S-sigma*eye(3)) Z;
    Z' sigma];
sigma = 0.5*trace(S);
kappa = trace(adjoint(S));
delta = det(S);
a = sigma^2-kappa;
b = sigma^2+Z'*Z;
c = delta + Z'*S*Z;
d = Z'*S^2*Z;
% 牛頓-拉夫森方法初始化
lambda_t = 1;  % 初始值（可以根據實際情況調整）
tolerance = 1e-6;  % 收斂容忍度
max_iter = 100;  % 最大迭代次數
iter = 0;  % 迭代計數器

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
    
    lambda_t = lambda_t1;  % 更新 lambda
    iter = iter + 1;
end

% 計算最優四元數
alpha = lambda_t^2 - sigma^2 + kappa;
beta = lambda_t - sigma;
X = (alpha * eye(3) + beta * S + S^2) * Z;
gamma = (lambda_t + sigma) * alpha - delta;

% 計算最優四元數 q_opt
q_opt = [X; gamma] / sqrt(gamma^2 + norm(X)^2);

% 顯示結果
disp('Optimal quaternion (q_opt):');
disp(q_opt);

% 計算最終的方向余弦矩陣 (DCM)
Cm_opt = EP2C(q_opt);

% 顯示方向余弦矩陣
disp('Direction Cosine Matrix (DCM) corresponding to the optimal quaternion:');
disp(Cm_opt);
disp('Difference between Cm_opt and BN:');
disp(norm(Cm_opt - BN));  
