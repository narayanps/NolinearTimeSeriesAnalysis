function [opd] = get_ordinal_pattern_dist(ts,emb_dim,tau)
%get_ordinal_pattern_dist - compute the mpr-complexity measure outlined in
%Keller and Sin 2005 and Rosso et al. 2012.

%INPUT - ts, time series 
%        emb_dim, embedding dimension
%       tau, delay

%OUTPUT - opd, ordinal pattern distribution

% AUTHOR - Narayan Puthanmadam Subramaniyam, 2015, Tampere Univ. of
% Technology

%begin of change log%

%end of change log%
%initialize some variables
opd = zeros(1,factorial(emb_dim));

%embed the time series
N = length(ts);
M = N-(emb_dim-1)*tau; 
psv_ts=zeros(M,emb_dim); 

for j=1:emb_dim
   psv_ts(:,j)=ts((1:M)+(j-1)*tau)';
end

% get the ordinal pattern
for k=1:1:M
    x = psv_ts(k,:);
    i = zeros(emb_dim,emb_dim);
    
for l=2:emb_dim
    for t=l:emb_dim
        i(t,l) = i(t-1,l-1);
        if (x(t) <= x(t-(l-1))) %|| (abs(x(t)- x(t-(l-1))) < 1e-10)
           i(t,l) = i(t,l)+1;
        end
    end
end

n_d = 0;
for l=2:emb_dim
   n_d = n_d + i(emb_dim,l)*(factorial(emb_dim)/factorial(l));
end

% n_d =i(emb_dim,2);
% for l=3:emb_dim
%     n_d = (l * n_d) + i(emb_dim, l);
% end

opd(n_d+1) = opd(n_d+1) + 1;
end
end
