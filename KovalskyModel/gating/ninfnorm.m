% K norm functions
function out = ninfnorm(V)
%     Vh = -34.12;
%     slope = 9.146;
%     out = 1/(1+exp((Vh-V)/slope));
    alpha = .01*vtrap(-(V+55),10);
    beta = .125*exp(-(V+65)/80);
    sum = alpha + beta;
    out = alpha/sum;
end