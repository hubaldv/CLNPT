function out = tauhf(V)
%     out = 0.246 + 1.63*exp(-0.5*((V+61.87)/15.25)^2);
    celsius = 20;
    q10 = 3^((celsius - 6.3)/10);   
    
    alpha = .07 * exp(-(V+59)/20);
    beta = 1 / (exp(-(V+29)/10) + 1);
    sum = alpha + beta;
    out = 1/(q10*sum);
end