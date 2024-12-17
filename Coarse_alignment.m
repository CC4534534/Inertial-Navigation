%% COARSE ALIGNMENT Simulation 粗對準模擬
% Increment of Angular rate [Gyro_x Gyro_y Gyro_z] Increment of Acceleration [acc_x; acc_y; acc_z]
clc;clear
%% 1. Parameter setting
% 指定文件路徑
filePath = 'D:\simulation\IMU_Simulation\data\Simulated_IMU.txt';

% 使用 readtable 讀取文件
data = readtable(filePath, 'Delimiter', '\t', 'ReadVariableNames', false);

% 將數據轉換為數值矩陣
IMU_data = table2array(data);
% 提取時間 (第一列)
time = IMU_data(:, 1);

% 提取角速度 (第二到第四列)
gyro_XYZ = IMU_data(:, 2:4); % [Gyro_X, Gyro_Y, Gyro_Z]

% 提取加速度 (第五到第七列)
acc_XYZ = IMU_data(:, 5:7); % [Acc_X, Acc_Y, Acc_Z]

% load IMU_real1;
Wibb0 = IMU_data(:, 2:4)';
Fibb0 = IMU_data(:, 5:7)';
t = 200;

L = 22.9998647916*pi/180;  %Local latitude
%lamda = 116.311*pi/180; %Longtitude
omega = 7.2921151467e-5;  %Earth angular rate
g = 9.7803267714*(1+0.00193185138639*sin(L))/sqrt(1-0.00669437999013*sin(L)*sin(L)); %Local Gravatational Acceleration, ingnore h
h = 0.01; %update frequency

t1 = 0.5*t;
t2 = t;
long = length(Wibb0);  %confirm the correponding points of t1 and t2
long_t1 = ceil(long*t1/t);
long_t2 = ceil(long*t2/t);

% data * sampling interval, change to incremental
Gyro = Wibb0(:,1:long_t2)*h;  %3-by-16000
acc = Fibb0(:,1:long_t2)*h;

Gyro_x = Gyro(1,:); %1-by-16000
Gyro_y = Gyro(2,:);
Gyro_z = Gyro(3,:);
acc_x = acc(1,:);
acc_y = acc(2,:);
acc_z = acc(3,:);

%% get C_en
C_en = [0 1 0; -sin(L) 0 cos(L); cos(L) 0 sin(L)];

%% get C_ie  
theta = omega*t2;
C_ie = [cos(theta) sin(theta) 0; 
    -sin(theta) cos(theta) 0; 
    0 0 1];

%% get C_bib0
Q = [1 0 0 0]';
 V = [0 0 0]';
for ii = 1:long_t2
    temp_x = Gyro_x(ii);
    temp_y = Gyro_y(ii);
    temp_z = Gyro_z(ii);
    mod = sqrt(temp_x^2+temp_y^2+temp_z^2);
    q0 = cos(mod/2);
    q1 = temp_x/mod*sin(mod/2);
    q2 = temp_y/mod*sin(mod/2);
    q3 = temp_z/mod*sin(mod/2);
    Q = [q0 -q1 -q2 -q3; q1 q0 q3 -q2; q2 -q3 q0 q1; q3 q2 -q1 q0]*Q; %(P269 9.3.40)
    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);
    Q3 = Q(4);
    C_bib0 = [Q0^2+Q1^2-Q2^2-Q3^2 2*(Q1*Q2-Q0*Q3) 2*(Q1*Q3+Q0*Q2);
        2*(Q1*Q2+Q0*Q3) Q0^2-Q1^2+Q2^2-Q3^2 2*(Q2*Q3-Q0*Q1);
        2*(Q1*Q3-Q0*Q2) 2*(Q2*Q3+Q0*Q1) Q0^2-Q1^2-Q2^2+Q3^2]; %(P251 9.2.34)

    V_x = acc_x(ii);
    V_y = acc_y(ii);
    V_z = acc_z(ii);
    V = V+C_bib0*[V_x; V_y; V_z];

    if ii == long_t1
        V_t1 = V;
    end
    if ii == long_t2
        V_t2 = V;
    end
end
C2 = [V_t1'; 
    V_t2'; 
    cross(V_t1,V_t2)'];

Vg_t1 = [g*cos(L)*sin(omega*t1)/omega; 
    g*cos(L)*(1-cos(omega*t1))/omega; 
    g*sin(L)*t1];
Vg_t2 = [g*cos(L)*sin(omega*t2)/omega; 
    g*cos(L)*(1-cos(omega*t2))/omega; 
    g*sin(L)*t2];
C1 = [Vg_t1';
    Vg_t2';
    cross(Vg_t1,Vg_t2)'];  %3-by-3
C_ib0i = C1\C2;

%% get C_bn
C_bn = C_en*C_ie*C_ib0i*C_bib0; %3-by-3
%% �����̬  
Theta = asind(C_bn(3,2));   
Gamma = atand(-C_bn(3,1)/C_bn(3,3));   %P252 9.2.41
Psi = -atand(C_bn(1,2)/C_bn(2,2));

%% get attitude
%% true value decision

if abs(C_bn(2,2))<10e-6  
    if C_bn(1,2)<0
        Psi = 90;
    else
        Psi = -90;
    end
else
    if C_bn(2,2)<0  %cos(Psi)*cos(Theta)
       if C_bn(1,2)<0  %sin(Psi)*cos(Theta)
           Psi = Psi+180;
       else
           Psi = Psi-180;
       end
    end
end
if C_bn(3,3)<0
    if Gamma>0
        Gamma = Gamma-180;
    else
        Gamma = Gamma+180;
    end
end

% 
real_Psi = 1+5*cos(t*2*pi/7+pi/3);
real_Theta =1+ 7*cos(t*2*pi/5+pi/4);
real_Gamma =1+ 10*cos(t*2*pi/6+pi/7);

Psi = Psi-real_Psi 
Theta = Theta-real_Theta %error
Gamma = Gamma-real_Gamma

