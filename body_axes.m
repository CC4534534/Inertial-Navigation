function C_nb = body_axes(Psi,Theta,Gamma)
% BODY_AXES ������̬��Ϣ����õ�������ϵ����������ϵ������ת�ƾ���
% Psi ����� ��Z�� !!Ҫ���Ϊ���� Radian
% Theta ������ ��Y��
% Gamma ����� ��X��

% 2017-12-27 �쳤�� �ࣨP252��

C_n1 = [cos(Psi) -sin(Psi) 0; sin(Psi) cos(Psi) 0; 0 0 1];
C_12 = [1 0 0; 0 cos(Theta) sin(Theta); 0 -sin(Theta) cos(Theta)];
C_2b = [cos(Gamma) 0 -sin(Gamma); 0 1 0; sin(Gamma) 0 cos(Gamma)];
C_nb = C_2b*C_12*C_n1;