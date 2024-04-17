%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bionomial expansion for Fractional order differencing operator to 
% generate the new time series Z
% Z(t)=(1-L)^d*X(t)
%
% Input arguments:
%
%   X          : the time series input
%
%   d          : the fractional order 
%
%   N          : the length of the expansion

% Output arguments: 
%
%   Z          : the new time series
% 
%
% (c) Yuankun Xue 2015
% University of Southern California
% All copyrights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z]= Frac_bin_expansion(X,d,N)

infit=N;
for i=1:length(X)
   Z(i)=X(i);
   for j=1:infit
      if(i>j)
      Z(i)=Z(i)+(gamma(-d+j)./(gamma(-d).*gamma(j+1))).*X(i-j);
      end
   end
end

end