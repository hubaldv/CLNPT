function out = hinff(V)
%     Vh = -56.39;
%     slope = -7.22;
%     out = 1/(1+exp((Vh-V)/slope));
    alpha = .07 * exp(-(V+59)/20);
    beta = 1 / (exp(-(V+29)/10) + 1);
    sum = alpha + beta;
    out = alpha/sum;
end
