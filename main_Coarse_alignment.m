
clc;clear;
% 指定文件路徑
filePath = 'D:\simulation\IMU_Simulation\data\Simulated_IMU.txt';

% 使用 readtable 讀取文件
data = readtable(filePath, 'Delimiter', '\t', 'ReadVariableNames', false);

% 將數據轉換為數值矩陣
IMU_data = table2array(data);

% 提取時間 (第一列)
time = IMU_data(:, 1);

% 取前 250 秒的數據
sampleRate = 0.01;  % IMU頻率100Hz
align_time = 250;
align_time_imudata = time + align_time;  % 對準時間為 250秒

% 根據時間範圍過濾數據
time_filter = time <= align_time_imudata;
IMU_data_filtered = IMU_data(time_filter, :);

% 提取樣本數量
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



