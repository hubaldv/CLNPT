function out = vtrap(x,y)
    if (abs(x/y) < 1e-6)
        out = y*(1 - x/y/2);
    else
        out = x/(exp(x/y) - 1);
    end
end