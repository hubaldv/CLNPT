function [X_bar] = map_model(X, sigma, zeta)
%MAP_MODEL Summary of this function goes here
%   Detailed explanation goes here

    global beta mu

    x = X(1);
    y = X(2);
    xBar = fAlpha(x,y+beta) + zeta;
    yBar = y-mu*(x+1-sigma);

    X_bar = [xBar; yBar];

end

function [f, p] = fAlpha(x,y) 
    global alpha
    if (x < (-1-alpha/2))
        f = -alpha^2/4-alpha + y;
        p = 1;
    elseif (((-1-alpha/2) <= x) && (x <= 0))
        f = alpha*x+(x+1)^2+y;
        p = 2;
    elseif ((0 < x) && (x < (y+1)))
        f = y+1;
        p = 3;
    elseif x >= (y+1)
        f = -1;
        p = 4;
    else
        disp('error')
        f = 0;
        p = 5;
    end
end
