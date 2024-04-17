function [signal_auto_cor,lags_vec, confidence_95] = autocorrelation(signal, NumLags)
%   Compute the sample autocorrelation function
%   Hubald Verzijl - 2021

%   Input: signal (input value), NumLags (number of lags to determine)

lags_vec = 0:NumLags;

signal_n = length(signal);
signal_mean = mean(signal);

signal_auto_cov = zeros(1,NumLags+1);
    for i = 0:NumLags
        sum = 0;
        for j = i+1:signal_n
            sum = sum + (signal(j)-signal_mean)*(signal(j-i)-signal_mean);
        end
        signal_auto_cov(i+1) = sum / signal_n;
    end

signal_auto_cor = signal_auto_cov/signal_auto_cov(1);

% 95% confidence interval
confidence_95 = 1.96/sqrt(length(signal));

end

