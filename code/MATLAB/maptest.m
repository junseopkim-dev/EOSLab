clear all; fclose all; %close all;
clc;

%%
% info = shapeinfo('HydroLAKES_polys_v10.shp'); % extract information from shapefile
% crs = info.CoordinateReferenceSystem; % verifies the shapefile uses latitude and longitude coordinates
% S = shaperead('HydroLAKES_polys_v10.shp','UseGeoCoords',true); % read the shapefile using geocoordinates

MAP_NAME = 'HydroLAKES_polys_v10.shp';

S = shaperead(MAP_NAME,'Selector',{@(v1) (v1 == 7),'Hylak_id'});

%%

% Add shapefile to geoscatter plot
hold on
    plot(S.X,S.Y);
hold off



% y = T(:,5);
% 
% figure,
% plot(x,y);