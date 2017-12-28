function C_eg = Geography_axes(lamda,L)
% GEOGRAPHY_AXES ������γ�ȣ�����������ϵ����������ϵ������任����
% lamda ����  !!Ҫ���Ϊ���� Radian
% L γ��
% C_eg ����任��

% 2017-12-27 �쳤�� ��(P198)

C_e1 = [cos(lamda) sin(lamda) 0; -sin(lamda) cos(lamda) 0; 0 0 1];
C_12 = [cos(pi/2-L) 0 -sin(pi/2-L); 0 1 0; sin(pi/2-L) 0 cos(pi/2-L)];
C_2g = [0 1 0; -1 0 0; 0 0 1];
C_eg = C_2g*C_12*C_e1;

