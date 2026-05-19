function [Y, k2] = learning (X_samp, k1, get_knn, rnn, id_samp, no_dims, initialize, agg_coef, T_epoch, varargin)
% Sparse-P landmark learning with two execution modes:
% dense keeps the original learning_s vectorized gradient path for speed;
% block keeps the original learning_l memory-bounded gradient path.

paramNames = {'Mode','BlockSize','SNNBlockSize','MemoryBudgetMB'};
defaults = {'dense',[],[],4096};
[mode, block_size, snn_block_size, memory_budget_mb] = internal.stats.parseArgs(paramNames, defaults, varargin{:});
mode = validatestring(mode, {'dense','block'}, mfilename, 'Mode');

[N, ~] = size(X_samp);
k2 = adaptive_landmark_neighbors(N);

if strcmpi(mode, 'dense')
    P = build_sparse_probability_matrix_dense(X_samp, k1, get_knn, rnn, id_samp, no_dims, agg_coef, k2);
else
    if isempty(snn_block_size)
        snn_block_size = memory_safe_snn_block_size(N, memory_budget_mb);
    end
    snn_block_size = max(1, min(N, snn_block_size));
    P = build_sparse_probability_matrix_block(X_samp, k1, get_knn, rnn, id_samp, no_dims, agg_coef, k2, snn_block_size);
end
P = sparse((P + P'))/2;
clear get_knn rnn id_samp

Y = initialize_embedding(P, X_samp, no_dims, initialize, mode);

P = P/(sum(P(:))-N);

if strcmpi(mode, 'dense')
    Y = optimize_dense(P, Y, T_epoch);
else
    if isempty(block_size)
        block_size = memory_safe_block_size(N);
    end
    block_size = max(1, min(N, block_size));
    Y = optimize_block(P, Y, T_epoch, block_size);
end
disp([num2str(T_epoch),' epochs have been computed!']);

end

function k2 = adaptive_landmark_neighbors(N)
if (N < 9)
    k2 = N;
elseif(N > 1000)
    k2 = ceil(log(N)/log(2)) + 18;
elseif(N > 50)
    k2 = ceil(0.02*N) + 8;
else
    k2 = 9;
end
end

function P = build_sparse_probability_matrix_dense(X_samp, k1, get_knn, rnn, id_samp, no_dims, agg_coef, k2)
[N, dim] = size(X_samp);
if(k1 > 0)
    SNN = zeros(N, N);
    G = get_knn(id_samp,:);
    knn_rnn_mat = rnn(G);
    for i = 1:N
        snn_id = ismember(G,G(i,:));
        nn_id = find(max(snn_id,[],2)==1);
        SNN(i,nn_id) = sum(knn_rnn_mat(nn_id,:).*snn_id(nn_id,:),2);
        SNN(i,:) = SNN(i,:)./max(max(SNN(i,:)),realmin);
    end
    Dis = (1-SNN).^agg_coef.*pdist2(X_samp,X_samp);
    [sort_dis, idx] = sort(Dis, 2);
    row = zeros(N*k2,1);
    col = zeros(N*k2,1);
    p_val = zeros(N*k2,1);
    for i = 1:N
        pos = (i-1)*k2+1:i*k2;
        row(pos) = i;
        col(pos) = idx(i,1:k2);
        p_val(pos) = exp(-0.5*sort_dis(i,1:k2).^2./max(mean(sort_dis(i,1:k2)).^2,realmin));
    end
    P = sparse(row,col,p_val,N,N);
else
    if(N >= 5000 && dim >= 50)
        xx = init_pca(X_samp, no_dims, 0.8);
        [samp_knn, samp_dis]= knnsearch(xx,xx,'k',k2);
    else
        [samp_knn, samp_dis]= knnsearch(X_samp,X_samp,'k',k2);
    end
    P = sparse(repmat((1:N)',k2,1),samp_knn(:),exp(-0.5*samp_dis.^2./max(mean(samp_dis,2).^2,realmin)),N,N);
end
end

function P = build_sparse_probability_matrix_block(X_samp, k1, get_knn, rnn, id_samp, no_dims, agg_coef, k2, snn_block_size)
[N, dim] = size(X_samp);
if(k1 > 0)
    row = zeros(N*k2,1);
    col = zeros(N*k2,1);
    p_val = zeros(N*k2,1);
    G = get_knn(id_samp,:);
    knn_rnn_mat = rnn(G);
    for start_row = 1:snn_block_size:N
        stop_row = min(start_row+snn_block_size-1,N);
        [row_block, col_block, p_block] = build_sparse_probability_matrix_block_rows( ...
            X_samp, G, knn_rnn_mat, agg_coef, k2, start_row, stop_row);
        pos = (start_row-1)*k2+1:stop_row*k2;
        row(pos) = row_block;
        col(pos) = col_block;
        p_val(pos) = p_block;
    end
    P = sparse(row,col,p_val,N,N);
else
    if(N >= 5000 && dim >= 50)
        xx = init_pca(X_samp, no_dims, 0.8);
        [samp_knn, samp_dis]= knnsearch(xx,xx,'k',k2);
    else
        [samp_knn, samp_dis]= knnsearch(X_samp,X_samp,'k',k2);
    end
    P = sparse(repmat((1:N)',k2,1),samp_knn(:),exp(-0.5*samp_dis.^2./max(mean(samp_dis,2).^2,realmin)),N,N);
end
end

function [row_block, col_block, p_block] = build_sparse_probability_matrix_block_rows(X_samp, G, knn_rnn_mat, agg_coef, k2, start_row, stop_row)
N = size(X_samp,1);
block_len = stop_row - start_row + 1;
row_block = zeros(block_len*k2,1);
col_block = zeros(block_len*k2,1);
p_block = zeros(block_len*k2,1);
D_block = pdist2(X_samp(start_row:stop_row,:),X_samp);
for local_i = 1:block_len
    i = start_row + local_i - 1;
    snn_id = ismember(G,G(i,:));
    nn_id = find(max(snn_id,[],2)==1);
    snn = zeros(1,N);
    snn(nn_id) = sum(knn_rnn_mat(nn_id,:).*snn_id(nn_id,:),2);
    [sort_dis, idx] = sort((1-snn./max(max(snn),realmin)).^agg_coef.*D_block(local_i,:),2);
    local_pos = (local_i-1)*k2+1:local_i*k2;
    row_block(local_pos) = i;
    col_block(local_pos) = idx(1:k2);
    p_block(local_pos) = exp(-0.5*sort_dis(1:k2).^2./max(mean(sort_dis(1:k2)).^2,realmin));
end
end

function Y = initialize_embedding(P, X_samp, no_dims, initialize, mode)
if strcmp(initialize,'le')
    if strcmpi(mode, 'block')
        degree = full(sum(P,1))';
        Dg = spdiags(degree, 0, size(P,1), size(P,2));
        Dg_sqrt = spdiags(sqrt(degree), 0, size(P,1), size(P,2));
        L = Dg_sqrt*(Dg - P)*Dg_sqrt;
    else
        Dg = diag(sum(P));
        L = Dg - P;
        L = Dg.^(0.5)*L*Dg.^(0.5);
    end
    [Y, ~] = eigs(L,no_dims+1,'smallestabs');
    Y(:,1)=[];
elseif strcmp(initialize,'pca')
    Y = pca(X_samp,no_dims);
elseif strcmp(initialize,'mds')
    Y = mds(X_samp,no_dims);
end
end

function Y = optimize_dense(P, Y, T_epoch)
[N, no_dims] = size(Y);
max_alpha = 2.5*N;
min_alpha = 2*N;
warm_step = 10;
preGrad = zeros(N,no_dims);
epoch = 1;
while epoch <= T_epoch
    if(epoch <= warm_step)
        alpha = max_alpha;
    else
        alpha = min_alpha + 0.5*(max_alpha-min_alpha)*(1+cos(pi*((epoch-warm_step)/(T_epoch-warm_step))));
    end
    D = pdist2(Y,Y).^2;
    Q1 = 1./(1+log(1+D));
    QQ1 = 1./(1+D);
    Q = Q1/(sum(Q1(:))-N);
    ProMatY = 4*(P-Q).*Q1.*QQ1;
    grad = (diag(sum(ProMatY))-ProMatY)*Y;
    Y = Y - alpha*(grad+(epoch-1)./(epoch+2)*preGrad);
    preGrad = grad;
    epoch = epoch + 1;
end
end

function Y = optimize_block(P, Y, T_epoch, block_size)
[N, no_dims] = size(Y);
no_blocks = ceil(N/block_size);
mark = zeros(no_blocks,2);
for i=1:no_blocks
   mark(i,:) = [(i-1)*block_size+1,min(i*block_size,N)];
end

max_alpha = 2.5*N;
min_alpha = 2*N;
warm_step = 10;
preGrad = zeros(N,no_dims);
epoch = 1;
while epoch <= T_epoch
    if(epoch <= warm_step)
        alpha = max_alpha;
    else
        alpha = min_alpha + 0.5*(max_alpha-min_alpha)*(1+cos(pi*((epoch-warm_step)/(T_epoch-warm_step))));
    end
    Pgrad = zeros(N,no_dims);
    Qgrad = zeros(N,no_dims);
    sumQ = 0;
    for i = 1:no_blocks
        idx = mark(i,1):mark(i,2);
        D_block = pdist2(Y(idx,:),Y).^2;
        Q1_block = 1./(1+log(1+D_block));
        QQ1_block = 1./(1+D_block);
        Pmat = -4*P(idx,:).*Q1_block.*QQ1_block;
        Qmat = -4*Q1_block.^2.*QQ1_block;
        len_blk = mark(i,2) - mark(i,1) + 1;
        idPQ = (len_blk*(mark(i,1)-1)+1):(len_blk+1):(len_blk*mark(i,2));
        Pmat(idPQ) = Pmat(idPQ) - sum(Pmat,2)';
        Qmat(idPQ) = Qmat(idPQ) - sum(Qmat,2)';
        Pgrad(idx,:) = Pmat*Y;
        Qgrad(idx,:) = Qmat*Y;
        sumQ = sumQ + sum(Q1_block(:));
    end
    grad = Pgrad - Qgrad/(sumQ-N);
    Y = Y - alpha*(grad + (epoch-1)/(epoch+2)*preGrad);
    preGrad = grad;
    epoch = epoch + 1;
end
end

function block_size = memory_safe_block_size(N)
max_block_entries = 8e6;
block_size = max(1, floor(max_block_entries/N));
block_size = min(N, block_size);
end

function block_size = memory_safe_snn_block_size(N, memory_budget_mb)
if isempty(memory_budget_mb) || isinf(memory_budget_mb)
    max_block_entries = 8e6;
else
    max_block_entries = max(1e6, memory_budget_mb * 1024^2 / 8 / 3);
end
block_size = max(1, floor(max_block_entries/N));
block_size = min(N, block_size);
end
