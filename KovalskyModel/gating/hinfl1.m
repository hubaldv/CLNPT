function out = hinfl1(V)
%     Vh = -72.50;
%     slope = -8.0;
%     out = 1/(1+exp((Vh-V)/slope));
    out = (1+exp((V+72.5)/8))^-1;
end