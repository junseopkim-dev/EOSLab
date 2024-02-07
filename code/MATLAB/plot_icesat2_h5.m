clear all; fclose all; %close all;
clc;

set_lon1 = -80;
set_lon2 = -75.5;
set_lat1 = 42.6;
set_lat2 = 44.6;

%% Open the HDF5 File.
% FILE_NAME = 'ATL13_20230813045046_08252001_006_01.h5';
% FILE_NAME = 'ATL13_20181019151241_03210101_006_01.h5';
FILE_NAME = 'ATL13_20190804013844_05650401_006_01.h5';
% FILE_NAME = 'ATL13_20181018202112_03090101_006_01.h5';
file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');


%%
% Open the datasets.
gt1l_lat_id=H5D.open(file_id, 'gt1l/segment_lat');
gt1l_lon_id=H5D.open(file_id, 'gt1l/segment_lon');
gt1l_data_id1=H5D.open(file_id, 'gt1l/ht_ortho');

gt1r_lat_id=H5D.open(file_id, 'gt1r/segment_lat');
gt1r_lon_id=H5D.open(file_id, 'gt1r/segment_lon');
gt1r_data_id1=H5D.open(file_id, 'gt1r/ht_ortho');

gt2l_lat_id=H5D.open(file_id, 'gt2l/segment_lat');
gt2l_lon_id=H5D.open(file_id, 'gt2l/segment_lon');
gt2l_data_id1=H5D.open(file_id, 'gt2l/ht_ortho');

gt2r_lat_id=H5D.open(file_id, 'gt2r/segment_lat');
gt2r_lon_id=H5D.open(file_id, 'gt2r/segment_lon');
gt2r_data_id1=H5D.open(file_id, 'gt2r/ht_ortho');

gt3l_lat_id=H5D.open(file_id, 'gt3l/segment_lat');
gt3l_lon_id=H5D.open(file_id, 'gt3l/segment_lon');
gt3l_data_id1=H5D.open(file_id, 'gt3l/ht_ortho');

gt3r_lat_id=H5D.open(file_id, 'gt3r/segment_lat');
gt3r_lon_id=H5D.open(file_id, 'gt3r/segment_lon');
gt3r_data_id1=H5D.open(file_id, 'gt3r/ht_ortho');


% DATAFIELD2_NAME='gt1r/ht_ortho';
% data_id2=H5D.open(file_id, DATAFIELD2_NAME);

% DATAFIELD3_NAME='gt1r/sigma_unit_fit';
% data_id3=H5D.open(file_id, DATAFIELD3_NAME);

%% Read the datasets.
gt1l_lat=H5D.read(gt1l_lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt1l_lon=H5D.read(gt1l_lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt1l_data=H5D.read(gt1l_data_id1,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

gt1r_lat=H5D.read(gt1r_lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt1r_lon=H5D.read(gt1r_lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt1r_data=H5D.read(gt1r_data_id1,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

gt2l_lat=H5D.read(gt2l_lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt2l_lon=H5D.read(gt2l_lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt2l_data=H5D.read(gt2l_data_id1,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

gt2r_lat=H5D.read(gt2r_lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt2r_lon=H5D.read(gt2r_lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt2r_data=H5D.read(gt2r_data_id1,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

gt3l_lat=H5D.read(gt3l_lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt3l_lon=H5D.read(gt3l_lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt3l_data=H5D.read(gt3l_data_id1,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

gt3r_lat=H5D.read(gt3r_lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt3r_lon=H5D.read(gt3r_lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL','H5P_DEFAULT');
gt3r_data=H5D.read(gt3r_data_id1,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
%% Plot world map coast line.
figure, 


box on;
coast = load('coastlines.mat');


% scatter(gt1l_lon, gt1l_lat, 10,gt1l_data);
% hold on
% scatter(gt1r_lon, gt1r_lat, 10,gt1r_data);
% scatter(gt2l_lon, gt2l_lat, 10,gt2l_data);
% scatter(gt2r_lon, gt2r_lat, 10,gt2r_data);
% scatter(gt3l_lon, gt3l_lat, 10,gt3l_data);
% scatter(gt3r_lon, gt3r_lat, 10,gt3r_data);
% colorbar();
% plot(coast.coastlon, coast.coastlat, 'k');

% Coordinates of the polygon
polygon_lon = [-77.2263, -79.1034, -79.8277, -80.2783, -79.2733, -78.4354, -77.4317, -76.2578, -75.8670, -75.9310, -76.6140, -77.2263];
polygon_lat = [44.2292, 43.9984, 43.5443, 43.1037, 43.0268, 43.2274, 43.0619, 43.2659, 43.6979, 44.1323, 44.3777, 44.2292];

% Filter coordinates based on whether they are inside the polygon
gt1l_filter = inpolygon(gt1l_lon, gt1l_lat, polygon_lon, polygon_lat);
gt1r_filter = inpolygon(gt1r_lon, gt1r_lat, polygon_lon, polygon_lat);
gt2l_filter = inpolygon(gt2l_lon, gt2l_lat, polygon_lon, polygon_lat);
gt2r_filter = inpolygon(gt2r_lon, gt2r_lat, polygon_lon, polygon_lat);
gt3l_filter = inpolygon(gt3l_lon, gt3l_lat, polygon_lon, polygon_lat);
gt3r_filter = inpolygon(gt3r_lon, gt3r_lat, polygon_lon, polygon_lat);

% Plot filtered coordinates
scatter(gt1l_lon(gt1l_filter), gt1l_lat(gt1l_filter), 10, gt1l_data(gt1l_filter));
hold on
scatter(gt1r_lon(gt1r_filter), gt1r_lat(gt1r_filter), 10, gt1r_data(gt1r_filter));
scatter(gt2l_lon(gt2l_filter), gt2l_lat(gt2l_filter), 10, gt2l_data(gt2l_filter));
scatter(gt2r_lon(gt2r_filter), gt2r_lat(gt2r_filter), 10, gt2r_data(gt2r_filter));
scatter(gt3l_lon(gt3l_filter), gt3l_lat(gt3l_filter), 10, gt3l_data(gt3l_filter));
scatter(gt3r_lon(gt3r_filter), gt3r_lat(gt3r_filter), 10, gt3r_data(gt3r_filter));
colorbar();
clim([75.00 75.45])
plot(coast.coastlon, coast.coastlat, 'k');


% axis([-180 180 -90 90]);
axis([set_lon1 set_lon2 set_lat1 set_lat2]);

