function [Y, id_samp] = scml(X, varargin)
%   This function returns representation of the N by D matrix X in the lower-dimensional space and 
%   the ID of landmarks sampled by PPS. Each row in X represents an observation.
% 
%   Parameters are: 
%
%   'NumDimensions'- A positive integer specifying the number of dimension of the representation Y. 
%                    Default: 2
%   'NumNeighbors' - A non-negative integer specifying the number of nearest neighbors for PPS to 
%                    sample landmarks. It must be smaller than N.
%                    Default: 20
%   'Normalize'    - Logical scalar. If true, normalize X using min-max normalization. If features in 
%                    X are on different scales, 'Normalize' should be set to true because the learning 
%                    process is based on nearest neighbors and features with large scales can override 
%                    the contribution of features with small scales. 
%                    Default: True
%   'InitMethod'   - A string specifying the method for initializing Y before manifold learning. 
%       'le'       - Laplacian eigenmaps.
%       'pca'      - Principal component analysis.
%       'mds'      - Multidimensional scaling.
%                    Default: 'le' 
%   'AggCoef'      - A positive scalar specifying the aggregation coefficient. 
%                    Default: 1.2
%   'MaxEpoch'     - Maximum number of epochs to take. 
%                    Default: 50
%   'TolVcc'       - Termination tolerance for variation coefficient of the last three KLD costs. 
%                    Default: 1e-7                   
                   
% Specify default parameters
paramNames = {'NumDimensions','NumNeighbors','Normalize','InitMethod','AggCoef','MaxEpoch','TolVcc'};
defaults   = {2,20,true,'le',1.2,50,1e-7};
[no_dims,k1,normalize,initialize,agg_coef,T_epoch,T_vcc] = internal.stats.parseArgs(paramNames, defaults, varargin{:});

% Obatin size and dimension of data
[n, dim] = size(X);

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
if(length(id_samp) <= 40000)
   [Y_samp, k2] = learning_s(X_samp,k1,get_knn,rnn,id_samp,no_dims,initialize,agg_coef,T_epoch,T_vcc);
else
    % Excessive data size may cause the probelm of memory overflow. Split the data into mutiple blocks
    % and compute the gradient block by block using 'learning_l' function
   [Y_samp, k2] = learning_l(X_samp,k1,get_knn,rnn,id_samp,no_dims,initialize,agg_coef,T_epoch,T_vcc);
end

% Compute embedding of non-landmarks
if(k1 > 0)
    id_rest = setdiff(1:n,id_samp);
    X_rest = X(id_rest,:);
    Y_rest = zeros(length(id_rest),no_dims);
    % Compute the optimal scale for each landmark
    scale = opt_scale(X_samp, Y_samp, k2);
    top_k = no_dims+1;
    near_samp = knnsearch(X_samp,X_rest,'k',top_k);
    for i=1:length(id_rest)
        near_top_k = near_samp(i,:);
        top_X = X_samp(near_top_k,:);
        top_Y = Y_samp(near_top_k,:);
        N_dis = pdist2(X_rest(i,:),top_X(1,:))*scale(near_top_k(1));
        % Perform CLLE
        Y_rest(i,:) = clle(top_X,top_Y,X_rest(i,:),N_dis);
    end
    Y = zeros(n,no_dims);
    Y(id_rest,:) = Y_rest;
    Y(id_samp,:) = Y_samp;
else
    Y = zeros(n,no_dims);
    Y(id_samp,:) = Y_samp;
end
end