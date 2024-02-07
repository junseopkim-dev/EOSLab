%% 초기 세팅
clear all; fclose all; close all;
clc;

FILE_NAME = 'ATL13_20190804013844_05650401_006_01.h5';
MAP_NAME = 'HydroLAKES_polys_v10.shp';
waterbodyid = 7; % 온타리오 호수의 아이디는 7

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
alpha = 0.05;
edges = [75.0:0.01:75.5];

histgraph(gt1l,waterbodyid,75.35,175);


xname = ["gt1l","gt1r","gt2l","gt2r","gt3l","gt3r"];

htortho1 = gt1l.ht_ortho.Value(gt1l_lo);

m1 = mean(htortho1);

% pd1 = makedist(htortho1)
pd1 = fitdist(htortho1, 'Normal');

% y = max(pd1)

%%

% function histgraph(gt,waterbodyid,x_coord, y_coord)
%     alpha = 0.05;
%     edges = [75.0:0.01:75.5];
%     lo = find(gt.inland_water_body_id.Value(:)==waterbodyid);
%     htortho = gt.ht_ortho.Value(lo);
% 
%     a1 = histogram(htortho, edges);
%     a2 = histfit(htortho);
% 
%     pd = fitdist(htortho, 'Normal')
%     confidence_interval = calculate_confidence_interval(pd.mu, pd.sigma, alpha);
% 
%     title(gt.Attributes.groundtrack_id);
% 
%     txt = sprintf('#%d \n%.2f ± %.4f', length(htortho), mean(htortho), mean(htortho) - confidence_interval(1));
%     text(x_coord, y_coord, txt, 'FontSize', 10);
%     disp(['95% Confidence Interval: [' num2str(confidence_interval(1)) ', ' num2str(confidence_interval(2)) ']']);
%     hold off;
% end

function histgraph(gt, waterbodyid, x_coord, y_coord)
    alpha = 0.05;
    edges = [75.0:0.01:75.5];
    lo = find(gt.inland_water_body_id.Value(:) == waterbodyid);
    htortho = gt.ht_ortho.Value(lo);
    
    % 히스토그램 및 정규분포 피팅
    a1 = histogram(htortho, edges);
    a2 = histfit(htortho);
    
    % 정규분포의 모수 추정
    pd = fitdist(htortho, 'Normal');
    confidence_interval = calculate_confidence_interval(pd.mu, pd.sigma, alpha);

    % 히스토그램 제목 및 텍스트 출력
    title(gt.Attributes.groundtrack_id);
    txt = sprintf('#%d \n%.2f ± %.4f', length(htortho), mean(htortho), mean(htortho) - confidence_interval(1));
    text(x_coord, y_coord, txt, 'FontSize', 10);
    hold on;

    % 정규분포의 95% 신뢰구간을 영역으로 색칠
    x_fill = linspace(confidence_interval(1), confidence_interval(2), 100);
    y_fill = normpdf(x_fill, pd.mu, pd.sigma);
    factor = length(htortho)/70;
    normalized_y_fill = y_fill * factor; % 정규화
    fill([x_fill, fliplr(x_fill)], [normalized_y_fill, zeros(size(normalized_y_fill))], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

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
