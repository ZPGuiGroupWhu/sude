% Input data
clear;
data = csvread('benchmarks/rice.csv');

% Obtain data size and true annotations
[~, m] = size(data);
ref = data(:, m);
X = data(:, 1:m-1);
clear data

% Perform SUDE embedding
t1 = clock;
[Y, idx, para] = sude(X,'NumNeighbors',10);
t2 = clock;
disp(['Elapsed time:', num2str(etime(t2,t1)),'s']);
plotcluster2(Y, ref);
axis off;