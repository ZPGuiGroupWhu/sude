function scale = opt_scale(X, Y, k_num)
% This function computes the optimal scales of landmarks.
%
%   Parameters are:
%
%   'X'       - N by D matrix. Each row in X represents high-dimensional features of a landmark.
%   'Y'       - N by d matrix. Each row in Y represents low-dimensional embedding of a landmark. 
%   'k_num'   - A non-negative integer specifying the number of KNN.

[n,~] = size(X);
[samp_knn,~] = knnsearch(Y,Y,'k',k_num+1);
scale = zeros(n,1);
for i=1:n
    XDis = pdist2(X(samp_knn(i,:),:),X(samp_knn(i,:),:));
    YDis = pdist2(Y(samp_knn(i,:),:),Y(samp_knn(i,:),:));
    scale(i) = sum(sum(XDis.*YDis))/max(sum(sum(XDis.^2)),realmin);
end