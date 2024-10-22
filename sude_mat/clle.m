function Y_i = clle(X_samp, Y_samp, X_i, N_dis)
%   Constrained Locally Linear Embedding (CLLE)
%   This function returns representation of point X_i
%
%   Parameters are:
%
%   'X_samp'   - High-dimensional features of KNN of point X_i. Each row denotes an observation.
%   'Y_samp'   - Low-dimensional embeddings of KNN of point X_i. 
%   'X_i'      - Current non-landmark point.
%   'N_dis'    - Distance between point X_i and its nearest neighbor in lower-dimensional space.

[n,~] = size(X_samp);
S = (X_samp-X_i)*(X_samp-X_i)';
if(abs(det(S)) <= eps)
    S = S + (0.1^2/n)*trace(S)*eye(n);
end
W = (inv(S)*ones(n,1))/(ones(1,n)*inv(S)*ones(n,1));
Y_0 = W'*Y_samp;
dd = pdist2(Y_samp(1,:),Y_0);
if(dd~=0)
    Y_i = Y_samp(1,:)+N_dis*(Y_0-Y_samp(1,:))/dd;
else
    Y_i = Y_samp(1,:);
end