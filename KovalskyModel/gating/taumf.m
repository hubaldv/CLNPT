function out = taumf(V)
%     out = 0.1092*exp(-0.5*((V+28.71)/25.5)^1.8);
%     out = 0.11;
    celsius = 20;
    q10 = 3^((celsius - 6.3)/10);

    alpha = .1 * vtrap(-(V+34),10);
    beta =  4 * exp(-(V+59)/18);
    sum = alpha + beta;
    out = 1/(q10*sum);
end