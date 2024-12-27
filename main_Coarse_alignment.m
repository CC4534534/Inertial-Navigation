
clc;clear;
% ָ���ļ�·��
filePath = 'D:\simulation\IMU_Simulation\data\Simulated_IMU.txt';

% ʹ�� readtable �xȡ�ļ�
data = readtable(filePath, 'Delimiter', '\t', 'ReadVariableNames', false);

% �������D�Q�锵ֵ���
IMU_data = table2array(data);

% ��ȡ�r�g (��һ��)
time = IMU_data(:, 1);

% ȡǰ 250 ��Ĕ���
sampleRate = 0.01;  % IMU�l��100Hz
align_time = 250;
align_time_imudata = time + align_time;  % ���ʕr�g�� 250��

% �����r�g�����^�V����
time_filter = time <= align_time_imudata;
IMU_data_filtered = IMU_data(time_filter, :);

% ��ȡ�ӱ�����
long = length(IMU_data_filtered);

Gyro = ones(3, long);
Acce = ones(3, long);
Gyro(1,:) = IMU_data_filtered(:, 3);
Gyro(2,:) = IMU_data_filtered(:, 2);
Gyro(3,:) = IMU_data_filtered(:, 4);
Acce(1,:) = IMU_data_filtered(:, 6);
Acce(2,:) = IMU_data_filtered(:, 5);
Acce(3,:) = IMU_data_filtered(:, 7);


i = 1;
[yaw(i),pitch(i),row(i)] = Coarse_alignment(Gyro,Acce,align_time)



