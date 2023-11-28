function xx = init_pca(X, no_dims, contri)
% This function preprocesses data with excessive size and dimensions.
%
%   Parameters are:
%
%   'X'       - N by D matrix. Each row in X represents an observation.
%   'no_dims' - A positive integer specifying the number of dimension of the representation Y. 
%   'contri'  - Threshold of PCA variance contribution.

[~, m] = size(X);

% Make sure data is zero mean
mapping.mean = mean(X, 1);
X = bsxfun(@minus, X, mapping.mean);

% Compute covariance matrix C
if size(X, 2) < size(X, 1)
    C = cov(X);
else
    C = (1 / size(X, 1)) * (X * X');        % if N>D, we better use this matrix for the eigendecomposition
end

% Perform eigendecomposition of C
C(isnan(C)) = 0;
C(isinf(C)) = 0;
[M, lambda] = eig(C);

% Sort eigenvectors in descending order
[lambda, ind] = sort(diag(lambda), 'descend');

% Obtain the best PCA dimension
if m < 2001
    bestDim = max(no_dims+1, find(cumsum(lambda)/sum(lambda) < contri, 1, 'last' ));
else
    bestDim = max(no_dims+1, find(cumsum(lambda)/sum(lambda(1:2000)) < contri, 1, 'last' ));
end

% Apply mapping on the data
xx = X*M(:,ind(1:bestDim));  