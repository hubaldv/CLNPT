% Na fast functions
function out = minff(V)
%     Vh = -34.12;
%     slope = 9.146;
%     out = 1/(1+exp((Vh-V)/slope));
    alpha = .1 * vtrap(-(V+34),10);
    beta =  4 * exp(-(V+59)/18);
    sum = alpha + beta;
    out = alpha/sum;
end