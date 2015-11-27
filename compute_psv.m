function [psa] = compute_psv(x,tau,m)
%UNTITLED3 Summary of this function goes here
%   x is a N X 1 vector (Samples of scalar measurement), tau is OptTau and
%   m is the optimal embedding dimension

% Author : Narayan P Subramaniyam
%           Department of Electronics and Communications Engineering
%           Tampere University of Technology, Finland

N = length(x);
M = N-(m-1)*tau; 


psa=zeros(M,m); 

for i=1:m
   psa(:,i)=x((1:M)+(i-1)*tau)';
end

end
