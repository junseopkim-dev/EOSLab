clear all; fclose all; %close all;
clc;

%% NOAA gauge (9/1/2018 to 8/31/2019, 365 days)
% https://tidesandcurrents.noaa.gov/stations.html?type=Water+Levels


% 2018 10 17 - 2020 01 15
% metric
% LST
% IGLD
%https://tidesandcurrents.noaa.gov/waterlevels.html?id=9052000&units=metric&bdate=20181017&edate=20200115&timezone=LST&datum=IGLD&interval=d&action=


% 1. 9052000 Cape Vincent, NY
% 2. 9052030 Oswego, NY
% 3. 9052058 Rochester, NY
% 4. 9052076 Olcott, NY

% Obs_CV = load('CO-OPS_9052000_met.txt');
% Obs_CVx = Obs_CV(:,1);
% Obs_CVy = Obs_CV(:,5);
% 
% figure,
% plot(Obs_CVx,Obs_CVy);
% Obs_Os
% Obs_Ro
% Obs_Ol

%% ICESAT-2 Data
filenames=textread('file_names.txt','%s');

yymmdd=zeros(length(filenames),3);
for j=1:length(filenames),
    yymmdd(j,1)=str2num(filenames{j}(7:10));
    yymmdd(j,2)=str2num(filenames{j}(11:12));
    yymmdd(j,3)=str2num(filenames{j}(13:14));
end

time2=yymmdd(:,1) + (yymmdd(:,2)-1)/12 + (yymmdd(:,3)/365.5) ;

%%
atl13=zeros(length(yymmdd),1);
%atlxx=zeros(length(yymmdd),1);

len = length(atl13);
% len = 10;

f = waitbar(0,'Please wait...');
for j=1:len,
    
    waitbar(j/len,f,'please wait')
    
    % temp=h5readall(strcat('../ATL13_rel002/',filenames{j}));
    temp = h5readall(filenames{j});
    
    temp2=temp.gt1l; 
    clear temp;
    
    lo=find(temp2.inland_water_body_id.Value(:)==7); % lake Ontario (HydroLake)
    
    if length(lo)>1,
          atl13(j,1)=median(temp2.ht_ortho.Value(lo));
    else
        atl13(j,1)=nan;

    end
    
end
close(f)

%% Bias Correction
% https://vdatum.noaa.gov/vdatumweb/vdatumweb?a=163303520191016
% https://www.ngs.noaa.gov/web/about_ngs/history/Zilkoski4.pdf
% datumBias=-0.4790; % NAD88 to EGM2008 (ICESat2), IGLD 1985 ???


% 알아낸 것 : https://tidesandcurrents.noaa.gov/datum_options.html
% IGLD 1985 : International Great Lakes Datum of 1985
% NAVD88 : North American Vertical Datum of 1988
% NAVD 88 and IGLD 85 are identical. (NAVD 88과 IGLD 85는 동일하다.)
% NAVD 88 bench mark values are given in Helmert orthometric height units while IGLD 85 values are in dynamic heights.


%% Comparision (Gauge-ICESAT2)





%%
figure;
x= time2;
y= atl13;
% plot(x,y,'o');
scatter(time2,atl13);

hold on
% plot(time2,temp2.ht_ortho.Value(lo))
% errorbar(x,y,0.1)