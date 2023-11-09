function embedX = mds (X, no_dims)
% Multidimensional scaling (MDS).
[n,~] = size(X);
D = pdist2(X,X).^2;
sumd = mean(D);
sumD = mean(D(:));
B = zeros(n,n);
for i=1:n
    for j=i+1:n
        B(i,j) = -0.5*(D(i,j)-sumd(i)-sumd(j)+sumD);
        B(j,i) = B(i,j);
    end
end
[U, value] = eigs(B,no_dims);
embedX = U*sqrt(abs(value));
end