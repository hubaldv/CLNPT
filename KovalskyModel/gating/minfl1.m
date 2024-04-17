% Na late 1 functions CHECK DONE!
function out = minfl1(V)
%     Vh = -25.29;
%     slope = 9.052;
%     out = 1/(1+exp((Vh-V)/slope));
    alpha = .1 * vtrap(-(V+25),10);
    beta =  4 * exp(-(V+50)/18);
    sum = alpha + beta;
    out = alpha/sum;
end