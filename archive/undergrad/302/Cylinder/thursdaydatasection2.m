%% 9:30 am lab group thursday Matrtin Kaden Liam Jeff



%% DATA SET NO SMOKE NO TAPE

% 130 rpm (5m/s) no smoke no tape
load('Cylinder_Lab_0ms_Set1.mat')
for i = 3:26
    pavg1_130 = mean(P(:,i))
end

% 240 rpm (10m/s) no smoke no tape
load('Cylinder_Lab_10_Set1.mat')
for i = 3:26
    pavg1_130 = mean(P(:,i))
end


%600 rpm (25 m/s) no smoke no tape
load('Cylinder_Lab_25_Set1.mat')
for i = 3:26
    pavg1_130 = mean(P(:,i))
end


%% DATA SET 2 SMOKE NO TAPE

% 130 rpm (5m/s) smoke no tape
load('')
for i = :
    pavg2_130 = mean(P(:,i))
end

% 240 rpm (10m/s) smoke no tape
load('')
for i = :
    pavg2_240 = mean(P(:,i))
end


%600 rpm (25 m/s) smoke no tape
load('')
for i = :
    pavg2_600 = mean(P(:,i))
end



%% DATA SET 3 SMOKE AND TAPE
% 130 rpm (5m/s) smoke  tape
load('')
for i = :
    pavg3_130 = mean(P(:,i))
end

% 240 rpm (10m/s) smoke  tape
load('')
for i = :
    pavg3_240 = mean(P(:,i))
end


%600 rpm (25 m/s) smoke  tape
load('')
for i = :
    pavg3_600 = mean(P(:,i))
end



%% DATA SET 4
load('')
for i = 3:26
    pavg4_130 = mean(P(:,i))
end

% 240 rpm (10m/s) smoke  tape
load('')
for i = 3:26
    pavg4_240 = mean(P(:,i))
end


%600 rpm (25 m/s) smoke  tape
load('')
for i = 3:26
    pavg4_600 = mean(P(:,i))
end
