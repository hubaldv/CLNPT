%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear regression kernel for STF model: (D(L)^a)X=AX+E
% 
%
% Input arguments:
%
%
%   p          : the autoregression order
%
%
%   d          : fractional differencing order vector
%
%   X          : M by N matrix of time series
% Output arguments: 
%
%   Z          : dependency matrix cell
% 
%
% (c) Yuankun Xue 2015
% University of Southern California
% All copyrights reserved
function Z=lin_kernel(X,p,d)
[M N]=size(X);
X_msr=X_gen(X,p,N); % correlation matrix
for i=1:p
	Z{i}=zeros(M,M);
end



for i=1:M
temp=Frac_bin_expansion(X(i,:),d(i),10);  % Get the Z series for one of the time series 
temp_coeff=regress(temp',X_msr);
for j=1:p
	index=zeros(1,M);
	for k=1:M
		index(k)=(k-1)*p+j;
	end
	Z{j}(i,:)=temp_coeff(index)';
end

end


end


