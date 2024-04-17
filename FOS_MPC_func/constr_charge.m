function [c,ceq] = constr_charge(x,M_half,M)

% This function defines an equality constraint for the pulse
% A positive or negative first pulse is possible. The pulse width is 50%. 

[n, ~] = size(x);
c = [];

% Use this part for const charge biphasic pulse
% ceq(1:M_half) = x(1) * ones(M_half,1) - eye(M_half)*x(1:M_half);
% second = -sum(x(1:M_half))/(M-M_half);
% ceq(M_half+1:M) = second*ones(M-M_half,1)-eye(M-M_half)*x(M_half+1:M);
% ceq(M+1:n) = x(M+1:n);

% Use this part for const charge arbitrary pulse
ceq(1) = sum(x(1:M));
if n ~= M
    diff = n-M-1;
    ceq(2:2+diff) = x(M+1:n);
end


end