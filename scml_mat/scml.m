function [Y, id_samp, para] = scml(X, varargin)
%   This function returns representation of the N by D matrix X in the lower-dimensional space and 
%   the ID of landmarks sampled by PPS. Each row in X represents an observation.
% 
%   Parameters are: 
% 
%   'NumDimensions'- A positive integer specifying the number of dimension of the representation Y. 
%                    Default: 2
%   'NumNeighbors' - A non-negative integer specifying the number of nearest neighbors for PPS to 
%                    sample landmarks. It must be smaller than N.
%                    Default: adaptive
%   'Normalize'    - Logical scalar. If true, normalize X using min-max normalization. If features in 
%                    X are on different scales, 'Normalize' should be set to true because the learning 
%                    process is based on nearest neighbors and features with large scales can override 
%                    the contribution of features with small scales. 
%                    Default: True
%   'LargeData'    - Logical scalar. If true, the data can be split into multiple blocks to avoid the problem 
%                    of memory overflow, and the gradient can be computed block by block using 'learning_l' function.                    
%                    Default: False
%   'InitMethod'   - A string specifying the method for initializing Y before manifold learning. 
%       'le'       - Laplacian eigenmaps.
%       'pca'      - Principal component analysis.
%       'mds'      - Multidimensional scaling.
%                    Default: 'le' 
%   'AggCoef'      - A positive scalar specifying the aggregation coefficient. 
%                    Default: 1.2
%   'MaxEpoch'     - Maximum number of epochs to take. 
%                    Default: 50                 

% Remove duplicate observations
[X, ~, orig_id] = unique(X, 'rows');

% Obtain size and dimension of data
[n, dim] = size(X);

% Specify default parameters
paramNames = {'NumDimensions','NumNeighbors','Normalize','LargeData','InitMethod','AggCoef','MaxEpoch'};
defaults   = {2,[],true,false,'le',1.2,50};
if(n>20000)
    defaults{2} = 50;
elseif(n>10000)
    defaults{2} = 20;
elseif(n>2000)
    defaults{2} = 10;
else
    defaults{2} = 0;
end
[no_dims, k1, normalize, large, initialize, agg_coef, T_epoch] = internal.stats.parseArgs(paramNames, defaults, varargin{:});
para = [paramNames;defaults];

% Normalize the data
if normalize
    X = mapminmax(X',0,1)';
end

% Perform PPS to obtain the landmarks
if(k1 > 0)
    if(n >= 5000 && dim >= 50)
        xx = init_pca(X, no_dims, 0.8);
        [get_knn, ~]= knnsearch(xx,xx,'k',k1+1);
    else
        [get_knn, ~]= knnsearch(X,X,'k',k1+1);
    end
    count = tabulate(get_knn(:));
    rnn = count(:,2);
    id_samp =  pps(get_knn, rnn, 1);
else
    get_knn = [];
    rnn = [];
    id_samp = 1:n;
end
X_samp = X(id_samp,:);

% Compute embedding of landmarks
if ~large
   [Y_samp, k2] = learning_s(X_samp, k1, get_knn, rnn, id_samp, no_dims, initialize, agg_coef, T_epoch);
else
   [Y_samp, k2] = learning_l(X_samp, k1, get_knn, rnn, id_samp, no_dims, initialize, agg_coef, T_epoch);
end

% Compute embedding of non-landmarks
if(k1 > 0)
    id_rest = setdiff(1:n,id_samp);
    X_rest = X(id_rest,:);
    Y_rest = zeros(length(id_rest),no_dims);
    % Compute the optimal scale for each landmark
    scale = opt_scale(X_samp, Y_samp, k2);
    top_k = no_dims+1;
    if(n >= 5000 && dim >= 50)
        [near_samp, near_dis] = knnsearch(xx(id_samp,:),xx(id_rest,:),'k',top_k);
    else
        [near_samp, near_dis] = knnsearch(X_samp,X_rest,'k',top_k);
    end
    for i=1:length(id_rest)
        near_top_k = near_samp(i,:);
        top_X = X_samp(near_top_k,:);
        top_Y = Y_samp(near_top_k,:);
        N_dis = near_dis(i,1)*scale(near_top_k(1));
        % Perform CLLE
        Y_rest(i,:) = clle(top_X,top_Y,X_rest(i,:),N_dis);
    end
    YY = zeros(n,no_dims);
    YY(id_rest,:) = Y_rest;
    YY(id_samp,:) = Y_samp;
else
    YY = Y_samp;
end

% Generate final result
Y = YY(orig_id,:);
