%% 초기 세팅
clear all; fclose all; close all;
clc;

FILE_NAME = 'ATL13_20190804013844_05650401_006_01.h5';
MAP_NAME = 'HydroLAKES_polys_v10.shp';
waterbodyid = 7; % 온타리오 호수의 아이디는 7

%%
cd C:\KJS\data\HydroLakeMap\HydroLAKES_polys_v10_shp
S = shaperead(MAP_NAME,'Selector',{@(v1) (v1 == waterbodyid),'Hylak_id'});

%% 날짜 추출
yymmdd=zeros(1,3);
yymmdd(1,1)=str2num(FILE_NAME(7:10));
yymmdd(1,2)=str2num(FILE_NAME(11:12));
yymmdd(1,3)=str2num(FILE_NAME(13:14));

time2=yymmdd(:,1) + (yymmdd(:,2)-1)/12 + (yymmdd(:,3)/365.5) ;

%%
cd C:\KJS\data\20240116
temp = h5readall(FILE_NAME);


gt1l = temp.gt1l;
gt1r = temp.gt1r;
gt2l = temp.gt2l;
gt2r = temp.gt2r;
gt3l = temp.gt3l;
gt3r = temp.gt3r;

gt1l_lo=find(gt1l.inland_water_body_id.Value(:)==waterbodyid); 
gt1r_lo=find(gt1r.inland_water_body_id.Value(:)==waterbodyid); 
gt2l_lo=find(gt2l.inland_water_body_id.Value(:)==waterbodyid); 
gt2r_lo=find(gt2r.inland_water_body_id.Value(:)==waterbodyid); 
gt3l_lo=find(gt3l.inland_water_body_id.Value(:)==waterbodyid); 
gt3r_lo=find(gt3r.inland_water_body_id.Value(:)==waterbodyid); 

%%
figure;

plot(S.X,S.Y,'k');

hold on;
plot_scatter(gt1l,waterbodyid,10);
plot_scatter(gt1r,waterbodyid,10);
plot_scatter(gt2l,waterbodyid,10);
plot_scatter(gt2r,waterbodyid,10);
plot_scatter(gt3l,waterbodyid,10);
plot_scatter(gt3r,waterbodyid,10);


colorbar;
txt = sprintf('%d.%d.%d',yymmdd(1,1),yymmdd(1,2),yymmdd(1,3));
title(txt);

%%
figure;
alpha = 0.05;
edges = [75.0:0.01:75.5];
t = tiledlayout(2,4);

nexttile
histgraph(gt1l,waterbodyid,75.35,175);
nexttile
histgraph(gt1r,waterbodyid,75.4,50);
nexttile
histgraph(gt2l,waterbodyid,75.32,60);
nexttile
histgraph(gt2r,waterbodyid,75.26,30);
nexttile
histgraph(gt3l,waterbodyid,75.16,60);
nexttile
histgraph(gt3r,waterbodyid,75.23,30);
nexttile([1 2])


xname = ["gt1l","gt1r","gt2l","gt2r","gt3l","gt3r","gtAll"];

htortho1 = gt1l.ht_ortho.Value(gt1l_lo);
htortho2 = gt1r.ht_ortho.Value(gt1r_lo);
htortho3 = gt2l.ht_ortho.Value(gt2l_lo);
htortho4 = gt2r.ht_ortho.Value(gt2r_lo);
htortho5 = gt3l.ht_ortho.Value(gt3l_lo);
htortho6 = gt3r.ht_ortho.Value(gt3r_lo);
htorthosum = vertcat(htortho1,htortho2,htortho3,htortho4,htortho5,htortho6);

% m1 = mean(htortho1);
% m2 = mean(htortho2);
% m3 = mean(htortho3);
% m4 = mean(htortho4);
% m5 = mean(htortho5);
% m6 = mean(htortho6);

pd1 = fitdist(htortho1, 'Normal');
pd2 = fitdist(htortho2, 'Normal');
pd3 = fitdist(htortho3, 'Normal');
pd4 = fitdist(htortho4, 'Normal');
pd5 = fitdist(htortho5, 'Normal');
pd6 = fitdist(htortho6, 'Normal');
pdsum = fitdist(htorthosum, 'Normal');

confidence_interval1 = calculate_confidence_interval(pd1.mu, pd1.sigma, alpha);
confidence_interval2 = calculate_confidence_interval(pd2.mu, pd2.sigma, alpha);
confidence_interval3 = calculate_confidence_interval(pd3.mu, pd3.sigma, alpha);
confidence_interval4 = calculate_confidence_interval(pd4.mu, pd4.sigma, alpha);
confidence_interval5 = calculate_confidence_interval(pd5.mu, pd5.sigma, alpha);
confidence_interval6 = calculate_confidence_interval(pd6.mu, pd6.sigma, alpha);
confidence_intervalsum = calculate_confidence_interval(pdsum.mu, pdsum.sigma, alpha);


% ci1 = paramci(pd1);
% ci2 = paramci(pd2);
% ci3 = paramci(pd3);
% ci4 = paramci(pd4);
% ci5 = paramci(pd5);
% ci6 = paramci(pd6);
% 
er1 = pd1.mu-confidence_interval1(1);
er2 = pd2.mu-confidence_interval2(1);
er3 = pd3.mu-confidence_interval3(1);
er4 = pd4.mu-confidence_interval4(1);
er5 = pd5.mu-confidence_interval5(1);
er6 = pd6.mu-confidence_interval6(1);
ersum = pdsum.mu-confidence_intervalsum(1);



x = [1,2,3,4,5,6,7];
y = [pd1.mu,pd2.mu,pd3.mu,pd4.mu,pd5.mu,pd6.mu,pdsum.mu];
err = [er1,er2,er3,er4,er5,er6,ersum];
% scatter(time2,gt1l.ht_ortho.Value(gt1l_lo));
plot(x,y);
errorbar(x,y,err);
set(gca,'xTick',x);
set(gca, 'xTickLabels',xname); 

%%
cd C:\KJS\data\20240116\ObservationStations
CapeVincent = readtimetable('CO-OPS_9052000_met.csv'); % 9052000 CapeVincent
Oswego = readtimetable('CO-OPS_9052030_met'); % 9052030 Oswego
Rochester = readtimetable('CO-OPS_9052058_met'); % 9052058 Rochester
Olcott = readtimetable('CO-OPS_9052076_met'); % 9052076 Olcott
% Latitude	43° 20.3 N
% Longitude	78° 43.6 W

figure,

plot(CapeVincent,"Verified_m_");
hold on,
plot(Oswego,"Verified_m_");
plot(Rochester,"Verified_m_");
plot(Olcott,"Verified_m_");

%%
function plot_scatter(gt,waterbodyid,marker_size)
    lo = find(gt.inland_water_body_id.Value(:)==waterbodyid); 
    scatter(gt.segment_lon.Value(lo),gt.segment_lat.Value(lo), marker_size, gt.ht_ortho.Value(lo));
end

%%

function histgraph(gt,waterbodyid,x_coord, y_coord)
    alpha = 0.05;
    edges = [75.0:0.01:75.5];
    lo = find(gt.inland_water_body_id.Value(:)==waterbodyid);
    htortho = gt.ht_ortho.Value(lo);

    a1 = histogram(htortho, edges);
    a2 = histfit(htortho);

    pd = fitdist(htortho, 'Normal')
    confidence_interval = calculate_confidence_interval(pd.mu, pd.sigma, alpha);

    title(gt.Attributes.groundtrack_id);

    txt = sprintf('#%d \n%.2f ± %.4f', length(htortho), mean(htortho), mean(htortho) - confidence_interval(1));
    text(x_coord, y_coord, txt, 'FontSize', 10);
    hold on;

    % 정규분포의 95% 신뢰구간을 area로 색칠
    x_fill = linspace(confidence_interval(1), confidence_interval(2), 100);
    y_fill = normpdf(x_fill, pd.mu, pd.sigma);
    fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % stem(confidence_interval(1));
    % stem(confidence_interval(2));
    % disp(['95% Confidence Interval: [' num2str(confidence_interval(1)) ', ' num2str(confidence_interval(2)) ']']);
    hold off;
end



% 정규분포의 95% 신뢰구간 계산
function confidence_interval = calculate_confidence_interval(mu, sigma, alpha)
    % z-score 계산
    z = norminv(1 - alpha/2, 0, 1);
    
    % 신뢰구간 계산
    margin_of_error = z * (sigma / sqrt(length(mu))); % 표본 크기가 크면 sqrt(n)을 사용
    confidence_interval = [mu - margin_of_error, mu + margin_of_error];
end

