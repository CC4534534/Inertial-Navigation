
% main()��������������ָ����ݣ����ö�׼����
% ע�⣺ʮ��IMU_real�����ļ� ������ʮ�β��������ÿһ�β���ʱ����Ϊ300s
% ��˶�׼����ʱ��t ��300s����
%�����0.03�����ڡ�

long = 40000;
t = 300;
n = long*10;
Wibb0 = ones(3,n);
Fibb0 = ones(3,n);
for ii= 1:10
    A = load (strcat('IMU_real',num2str(ii),'.mat'));
    Wibb0(:,long*(ii-1)+1:long*ii) = A.Wibb0(:,2:end);
    Fibb0(:,long*(ii-1)+1:long*ii) = A.Fibb0(:,2:end);
end
Sample = long/300*t;
Sample_Wibb0 = ones(3,Sample);
Sample_Fibb0 = ones(3,Sample);
count = long*10/Sample;
Psi = zeros(1,count);
Theta = zeros(1,count);
Gamma = zeros(1,count);
for jj = 1:count
    Sample_Wibb0 = Wibb0(:,Sample*(jj-1)+1:Sample*jj);
    Sample_Fibb0 = Fibb0(:,Sample*(jj-1)+1:Sample*jj);
    [Psi(jj),Theta(jj),Gamma(jj)] = Coarse_alignment(Sample_Wibb0,Sample_Fibb0,t);
end

subplot(3,1,1);
plot(Psi);
title('Psi(n),�����');
subplot(3,1,2);
plot(Theta);
title('Theta(n)��������');
subplot(3,1,3);
plot(Gamma);
title('Gamma(n)�������');

