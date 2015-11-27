function [C_js,H_s] = get_mpr_complexity(opd)
%mpr_complexity - compute the mpr-complexity measure outlined in Rosso et
%al. 2007 and Rosso et al. 2012.

%INPUT - opd, ordinal pattern distribution obtained using
%get_ordinal_pattern_dist function

%OUTPUT - C_js. MPR complexity and H_s, normalized shannon entropy

% AUTHOR - Narayan Puthanmadam Subramaniyam, 2015, Tampere Univ. of
% Technology

%begin of change log%

%end of change log%


% initialize some variables
M = length(opd);
P_e = repmat(1/M,1,M);
P = opd/sum(opd);

%compute the normalized shannon entropic measure H_s as per Rosso et al. 2007, PRL 99, 154102 (2007) 

H_s = -sum(P(P>0) .* log2(P(P>0)))/ log2(length(P));

    

%compute disequillibrium term Q_j = Q_0*J, where J is the Jensen-Shannon divergence. Refer Rosso et al. 2012, Physica A 391(2012)42-55
A = 0.5*(P_e + P);
S_A= -sum(A(A>0) .* log2(A(A>0)));
S_P = -sum(P(P>0) .* log2(P(P>0)));
S_P_e= -sum(P_e(P_e>0) .* log2(P_e(P_e>0)));

J = S_A - 0.5*S_P - 0.5*S_P_e ;
Q_0 = -2*(1/(((M+1)/M)*log2(M+1) - 2*log2(2*M) + log2(M)));
Q_j = Q_0*J;

%MPR statistical complexity C_js = Q_j*H_s
C_js = Q_j*H_s;

end
