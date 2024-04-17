%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear regression estimation for vector ARFIMA model
% Data preprocess--Matrix X generation
% generate the new time series Z
% Z(t)=X*W+E
%
% Input arguments:
%
%   Y          : the K-element time series vector
%
%   p          : the autoregression order
%
%   N          : the length of observations

% Output arguments: 
%
%   X          : KP*N matrix
% 
%
% (c) Yuankun Xue 2015
% University of Southern California
% All copyrights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X]= X_gen(Y,p,N)

K=length(Y(:,1)');
X=zeros(N,K*p);
for i=1:N
	for j=1:K
		for t=1:p
			if(i>t)
				X(i,(j-1)*p+t)=Y(j,i-t);
			end
		end
	end
end





end