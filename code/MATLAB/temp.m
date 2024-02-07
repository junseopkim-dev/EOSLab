clear all; fclose all; close all;
clc;

% 사용 예제
mu = 30; % 평균
sigma = 1; % 표준편차
alpha = 0.05; % 95% 신뢰수준

confidence_interval1 = calculate_confidence_interval(mu, sigma, alpha);
disp(['95% Confidence Interval: [' num2str(confidence_interval1(1)) ', ' num2str(confidence_interval1(2)) ']']);

% 정규분포의 95% 신뢰구간 계산
function confidence_interval = calculate_confidence_interval(mu, sigma, alpha)
    % z-score 계산
    z = norminv(1 - alpha/2, 0, 1);
    
    % 신뢰구간 계산
    margin_of_error = z * (sigma / sqrt(length(mu))); % 표본 크기가 크면 sqrt(n)을 사용
    confidence_interval = [mu - margin_of_error, mu + margin_of_error];
end