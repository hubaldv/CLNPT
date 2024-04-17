function [c,ceq] = constr_pulse(x,M_half)

% This function defines an equality constraint for the pulse
% A positive or negative first pulse is possible. The pulse width is 50%. 

[n, ~] = size(x);
c = [];
 
ceq(1:M_half) = x(1) * ones(M_half,1) - eye(M_half)*x(1:M_half);
ceq(M_half+1:M_half*2) = -x(1) * ones(M_half,1) - eye(M_half)*x(M_half+1:M_half*2);
ceq(M_half*2+1:n) = x(M_half*2+1:n);

end