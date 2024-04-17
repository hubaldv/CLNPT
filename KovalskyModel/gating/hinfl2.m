function out = hinfl2(V)
%     Vh = -55.67;
%     slope = -6.552;
%     out = 0.9827/(1+exp((Vh-V)/slope));
    out = 0.9827/(1 + exp(-((V + 55.67)/-6.552)));
end