clear
clc
close all

% Bed profile extraction

file_name_SB = "SWL_SB_With_SLR.nc"

x_grd = ncread(file_name_SB, "globalx"); % 'local x coordinate
data_bedlevel = ncread(file_name_SB, "zb"); % bed level

x_00 = data_bedlevel(:,:,1); % t=0
x_10 = data_bedlevel(:,:,2); % t=10 minutes
x_20 = data_bedlevel(:,:,3); % t=20 minutes
x_30 = data_bedlevel(:,:,4); % t=30 minutes
x_40 = data_bedlevel(:,:,5); % t=40 minutes
x_50 = data_bedlevel(:,:,6); % t=50 minutes
x_60 = data_bedlevel(:,:,7); % t=60 minutes
x_70 = data_bedlevel(:,:,8); % t=70 minutes
x_80 = data_bedlevel(:,:,9); % t=80 minutes
x_90 = data_bedlevel(:,:,10); % t=90 minutes
x_100 = data_bedlevel(:,:,11); % t=100 minutes
x_110 = data_bedlevel(:,:,12); % t=110 minutes
x_120 = data_bedlevel(:,:,13); % t=120 minutes
x_240 = data_bedlevel(:,:,25); % t=240 minutes
x_360 = data_bedlevel(:,:,37); % t=360 minutes
x_480 = data_bedlevel(:,:,49); % t=480 minutes
x_600 = data_bedlevel(:,:,61); % t=600 minutes


Table_bed_profiles(1,:) = x_00;
Table_bed_profiles(2,:) = x_10;
Table_bed_profiles(3,:) = x_20;
Table_bed_profiles(4,:) = x_30;
Table_bed_profiles(5,:) = x_40;
Table_bed_profiles(6,:) = x_50;
Table_bed_profiles(7,:) = x_60;
Table_bed_profiles(8,:) = x_70;
Table_bed_profiles(9,:) = x_80;
Table_bed_profiles(10,:) = x_90;
Table_bed_profiles(11,:) = x_100;
Table_bed_profiles(12,:) = x_110;
Table_bed_profiles(13,:) = x_120;
Table_bed_profiles(14,:) = x_240;
Table_bed_profiles(15,:) = x_360;
Table_bed_profiles(16,:) = x_480;
Table_bed_profiles(17,:) = x_600;


figure(1);
plot(x_grd,x_00)
hold on
plot(x_grd,x_600)
hold off

