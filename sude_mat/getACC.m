function [ACC] = getACC(ref, clus)
%   This function returns the clustering and classification accuracy.
n = length(ref);
p = unique(ref');
c = unique(clus');
P_size = length(p);
C_size = length(c);
Pid = double(ones(P_size,1)*ref' == p'*ones(1,n) );
Cid = double(ones(C_size,1)*clus' == c'*ones(1,n) );
CP = Cid*Pid';
[~,cost] = munkres(-CP);
ACC = -cost/n;