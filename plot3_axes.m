function plot3_axes(lamda,L,C_eg,C_nb)
% PLOT3_AXES ������ͼ����������ϵ����������ϵ����������ϵ
% lamda��L ����γ��
% C_eg eϵ��gϵ������任��
% C_nb nϵ��bϵ������任��

% 2017��12��27�� �쳤�� ��
%{
%quiver3(X,Y,Z,U,V,W)������Ϊͬ����С�ľ���,����ϸӦ�òμ�help quiver3��
%���� 2 ����ָ��ʸ������Ϊ2�� 'filled'��䣻 �߿�
% 'k'Ϊ��ɫ==> x�᣻ 'b'Ϊ��ɫ==>z�᣻ 'g'Ϊ��ɫ==>y��
%}

%% ��һ������������earth �����뾶R=2
R = 2;
[e_x,e_y,e_z] = sphere(20);
mesh(R*e_x,R*e_y,R*e_z);
hold on

% 4. �÷��������������Ϊ���������һ�㡣
%% ��������ϵe��% x��Ϊ [1 0 0]'; y��Ϊ[0 1 0]'; z��Ϊ[0 0 1]' 
quiver3(0,0,0,1,0,0,R+2,'k','filled','LineWidth',2);
hold on
quiver3(0,0,0,0,1,0,R+2,'g','filled','LineWidth',2);
hold on
quiver3(0,0,0,0,0,1,R+2,'b','filled','LineWidth',2);
hold on

%% ��������ϵg
C_eg = (C_eg)';   %���������У�inv(C_eg) = (C_eg)';                    
[b_x,b_y,b_z] = sph2cart(lamda,L,R); 

%����ϵ���ܱ�gϵ�ڸǣ�Ϊ�˹۲죬����ע������һ��
%%{
quiver3(b_x,b_y,b_z,C_eg(1,1),C_eg(2,1),C_eg(3,1),2,'k','filled','LineWidth',2);
hold on
quiver3(b_x,b_y,b_z,C_eg(1,2),C_eg(2,2),C_eg(3,2),2,'g','filled','LineWidth',2);
hold on 
quiver3(b_x,b_y,b_z,C_eg(1,3),C_eg(2,3),C_eg(3,3),2,'b','filled','LineWidth',2);
hold on
%%}

%% ��������ϵb
C_eb = C_nb*C_eg;
C_eb = (C_eb)';

quiver3(b_x,b_y,b_z,C_eb(1,1),C_eb(2,1),C_eb(3,1),1,'k','filled','LineWidth',1);
hold on
quiver3(b_x,b_y,b_z,C_eb(1,2),C_eb(2,2),C_eb(3,2),1,'g','filled','LineWidth',1);
hold on 
quiver3(b_x,b_y,b_z,C_eb(1,3),C_eb(2,3),C_eb(3,3),1,'b','filled','LineWidth',1);
hold on

axis equal
title('eϵ��gϵ����任');   
xlabel('X��');      
ylabel('Y��');      
zlabel('Z��');         
view(90+lamda*90/pi,30) %�ӵ�
hold off

%3. ʵ�ָ��� ����lamda��γ��L����������ϵe����������ϵg�ı任--�任����Ϊ C_eg
%ȱ�ݣ�δ���ڼ�ͷ�ϱ�עx,y,z������ɫ����
%{
X = zeros(3,1);
Y = X;
Z = X;
U = [1 0 0]';
V = [0 1 0]';
W = [0 0 1]';
quiver3(X,Y,Z,U,V,W,2,'k','filled','LineWidth',2);
hold on

lamda = pi/4;
L = pi/4;
C_eg = [-sin(lamda) cos(lamda) 0; -sin(L)*cos(lamda) -sin(lamda)*sin(L) cos(L); cos(L)*cos(lamda) cos(L)*sin(lamda) sin(L)];
U_g = C_eg*U;
V_g = C_eg*V;
W_g = C_eg*W;
quiver3(X,Y,Z,U_g,V_g,W_g,2,'b','filled','LineWidth',2);
hold off
%}
%2. ��Ϊ�����ʾ�µĻ���
%{
X = zeros(3,1);
Y = X;
Z = X;
U = [1 0 0]';
V = [0 1 0]';
W = [0 0 1]';
quiver3(X,Y,Z,U,V,W,2,'k','filled','LineWidth',2);
%}
%1. ��Ϊ�ٷ�quiver3����
%{
[x,y] = meshgrid(-2:.2:2,-1:.15:1);
z = x.*exp(-x.^2-y.^2);
[u,v,w] = surfnorm(x,y,z);
quiver3(x,y,z,u,v,w);
hold on
surf(x,y,z)
hold off
%}