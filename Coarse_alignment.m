function [Psi,Theta,Gamma] = Coarse_alignment(Wibb0,Fibb0,t)
%% COARSE ALIGNMENT 
%% 1. Parameter setting
L = 22.9998647916*pi/180;  %Local latitude from Chin-Hsin's /nspo-tc-insgnss/tasa-share/simulation demo data REF_Vehicle 
omega = 7.2921151467e-5;  %Earth angular rate
g = 9.7803253359*(1+0.001931852652*sin(L))/sqrt(1-0.0066943799901*sin(L)*sin(L)); %Local Gravatational Acceleration, ingnore height
h = 0.01; %update frequency

t1 = 0.5*t;
t2 = t;
long = length(Wibb0);  %confirm the correponding points of t1 and t2
long_t1 = ceil(long*t1/t);
long_t2 = ceil(long*t2/t);

% data * sampling interval, change to incremental
Gyro = Wibb0(:,1:long_t2)*h;  %3 * (sampling interval*data)
acc = Fibb0(:,1:long_t2)*h;

Gyro_x = Gyro(1,:); %1 * (sampling interval*data)
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

%%  get C_bb0
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
    Q = [q0 -q1 -q2 -q3; q1 q0 q3 -q2; q2 -q3 q0 q1; q3 q2 -q1 q0]*Q;
    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);
    Q3 = Q(4);
    C_bb0 = [Q0^2+Q1^2-Q2^2-Q3^2 2*(Q1*Q2-Q0*Q3) 2*(Q1*Q3+Q0*Q2);
        2*(Q1*Q2+Q0*Q3) Q0^2-Q1^2+Q2^2-Q3^2 2*(Q2*Q3-Q0*Q1);
        2*(Q1*Q3-Q0*Q2) 2*(Q2*Q3+Q0*Q1) Q0^2-Q1^2-Q2^2+Q3^2];
        
    V_x = acc_x(ii);
    V_y = acc_y(ii);
    V_z = acc_z(ii);
    V = V+C_bb0*[V_x; V_y; V_z];
    
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
    cross(Vg_t1,Vg_t2)'];  %3*3
C_ib0i = C1\C2;

%% get C_bn
C_bn = C_en*C_ie*C_ib0i*C_bb0; %3*3
%% 轉換後的角度
Theta = asind(C_bn(3,2));   
Gamma = atand(-C_bn(3,1)/C_bn(3,3));
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

% True value
real_Psi = 0;
real_Theta = 0;
real_Gamma = 0;

Psi = Psi-real_Psi; 
Theta = Theta-real_Theta; %error
Gamma = Gamma-real_Gamma;
