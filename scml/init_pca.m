function xx = init_pca(X, no_dims, contri)
% This function preprocesses data with excessive size and dimensions.
%   'X'       - N by D matrix. Each row in X represents an observation.
%   'no_dims' - A positive integer specifying the number of dimension of the representation Y. 
%   'contri'  - Threshold of PCA variance contribution.

[~, m] = size(X);
if m < 2001
    [~, map] = pca(X, m);
    bestDim = max(no_dims+1, find(cumsum(map.lambda)/sum(map.lambda) < contri, 1, 'last' ));
    xx = X*map.M(:,1:bestDim);
else
    X = pca(X, 2000);
    xx = init_pca(X, no_dims, contri);
end