function [Y, id_samp, para] = sude(X, varargin)
%   This function returns representation of the N by D matrix X in the lower-dimensional space and 
%   the ID of landmarks sampled by PPS. Each row in X represents an observation.
% 
%   This optimized copy uses sparse-P landmark learning by default.
%   Optional optimized parameter:
%   'LearningMode' - 'sparse-p' follows the baseline LargeData routing;
%                    'auto' also switches to block for large landmark sets;
%                    'dense' and 'block' force either path.
%   'BlockSize'    - Positive integer controlling block gradient memory use.
%   'SNNBlockSize' - Positive integer controlling block SNN construction.
%   'MemoryBudgetMB' - Dense learning memory budget. Sparse-p/auto switch to
%                      block when estimated dense memory exceeds this budget.

% Remove optimized-only arguments before parsing the legacy SUDE parameters.
learning_mode = '';
block_size = [];
snn_block_size = [];
memory_budget_mb = [];
[varargin, learning_mode] = extract_name_value(varargin, 'LearningMode', learning_mode);
[varargin, block_size] = extract_name_value(varargin, 'BlockSize', block_size);
[varargin, snn_block_size] = extract_name_value(varargin, 'SNNBlockSize', snn_block_size);
[varargin, memory_budget_mb] = extract_name_value(varargin, 'MemoryBudgetMB', memory_budget_mb);
if ~isempty(learning_mode)
    learning_mode = validatestring(learning_mode, ...
        {'sparse-p','auto','dense','block'}, mfilename, 'LearningMode');
end
memory_budget_mb = resolve_memory_budget_mb(memory_budget_mb);

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
learning_mode = resolve_learning_mode(learning_mode, large, size(X_samp,1), k1, memory_budget_mb);
[Y_samp, k2] = learning(X_samp, k1, get_knn, rnn, id_samp, no_dims, initialize, agg_coef, T_epoch, ...
    'Mode', learning_mode, 'BlockSize', block_size, 'SNNBlockSize', snn_block_size, ...
    'MemoryBudgetMB', memory_budget_mb);

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

end

function mode = resolve_learning_mode(mode, large, n_landmarks, k1, memory_budget_mb)
if isempty(mode) || strcmpi(mode, 'sparse-p')
    if large || exceeds_dense_memory_budget(n_landmarks, k1, memory_budget_mb)
        mode = 'block';
    else
        mode = 'dense';
    end
elseif strcmpi(mode, 'auto')
    if large || should_use_block_learning(n_landmarks) || exceeds_dense_memory_budget(n_landmarks, k1, memory_budget_mb)
        mode = 'block';
    else
        mode = 'dense';
    end
end
end

function tf = should_use_block_learning(n_landmarks)
dense_limit = 12000;
tf = n_landmarks > dense_limit;
end

function tf = exceeds_dense_memory_budget(n_landmarks, k1, memory_budget_mb)
if isempty(memory_budget_mb) || isinf(memory_budget_mb)
    tf = false;
    return
end
tf = estimate_dense_learning_memory_mb(n_landmarks, k1) > memory_budget_mb;
end

function memory_mb = estimate_dense_learning_memory_mb(n_landmarks, k1)
% Dense mode can hold SNN/Dis/sort buffers and several gradient matrices.
dense_matrices = 12;
if k1 <= 0
    dense_matrices = 8;
end
memory_mb = dense_matrices * n_landmarks * n_landmarks * 8 / 1024^2;
end

function memory_budget_mb = resolve_memory_budget_mb(memory_budget_mb)
default_budget_mb = 4096;
if isempty(memory_budget_mb)
    memory_budget_mb = default_budget_mb;
end
if isinf(memory_budget_mb)
    return
end
try
    user_mem = memory;
    available_mb = user_mem.MemAvailableAllArrays / 1024^2;
    memory_budget_mb = min(memory_budget_mb, max(512, 0.45*available_mb));
catch
end
end

function [args, value] = extract_name_value(args, name, default_value)
value = default_value;
idx = [];
for i = 1:2:length(args)
    if ischar(args{i}) && strcmpi(args{i}, name)
        idx = i;
        value = args{i+1};
        break;
    end
end
if ~isempty(idx)
    args(idx:idx+1) = [];
end
end
