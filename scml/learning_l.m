function [Y, k2] = learning_l (X_samp, k1, get_knn, rnn, id_samp, no_dims, initialize, agg_coef, T_epoch, T_vcc)
% This function returns representation of the landmarks in the lower-dimensional space and the number of 
% nearest neighbors of landamrks. It computes the gradient using probability matrix P and Q of data blocks.

% Obtain size and dimension of landmarks
[N, dim] = size(X_samp);

% Compute the number of nearest neighbors adaptively
if (N < 9)
    k2 = N;
else
    if(N > 1000)
        k2 = ceil(log(N)/log(2))+18;
    elseif(N > 50)
        k2 = ceil(0.02*N) + 8;
    else
        k2 = 9;
    end
end

% Compute high-dimensional probability matrix P
if(k1 > 0)
    row = [];
    col = [];
    p_val = [];
    knn_rnn_mat = rnn(get_knn(id_samp,:)); 
    for i = 1:N
        snn_id = ismember(get_knn(id_samp,:),get_knn(id_samp(i),:));
        nn_id = find(max(snn_id,[],2)==1);
        snn = zeros(1,N);
        snn(nn_id) = sum(knn_rnn_mat(nn_id,:).*snn_id(nn_id,:),2);
        [sort_dis, idx] = sort((1-snn./max(max(snn),realmin)).^agg_coef.*pdist2(X_samp(i,:),X_samp),2);
        p_val = [p_val;exp(-0.5*sort_dis(1:k2).^2./max(mean(sort_dis(1:k2)).^2,realmin))'];
        row = [row;i*ones(k2,1)];
        col = [col;idx(1:k2)'];
    end
    P = sparse(row,col,p_val);
else
    if(N >= 5000 && dim >= 50)
        xx = init_pca(X_samp, no_dims, 0.8);
        [samp_knn, samp_dis]= knnsearch(xx,xx,'k',k2);
    else
        [samp_knn, samp_dis]= knnsearch(X_samp,X_samp,'k',k2);
    end    
    P = sparse(repmat((1:N)',k2,1),samp_knn(:),exp(-0.5*samp_dis.^2./max(mean(samp_dis,2).^2,realmin)));
end
% Symmetrize matrix P
P = (P + P')/2;
clear xx get_knn rnn id_samp samp_knn samp_dis snn_id knn_rnn_mat

% Initialize embedding Y of landmarks
if strcmp(initialize,'le')
    Dg = diag(sum(P));
    L = Dg - P;
    L = Dg.^(0.5)*L*Dg.^(0.5);
    [Y, ~] = eigs(L,no_dims+1,'smallestabs');
    Y(:,1)=[];
elseif strcmp(initialize,'pca')
    Y = pca(X_samp,no_dims);
elseif strcmp(initialize,'mds')
    Y = mds(X_samp,no_dims);
end

% Normalize matrix P
P = P/(sum(P(:))-N);
clear Dg L

% Compute the start and end markers of each data block
no_blocks = ceil(N/10);
mark = zeros(no_blocks,2);
for i=1:no_blocks
   mark(i,:) = [(i-1)*ceil(N/no_blocks)+1,min(i*ceil(N/no_blocks),N)];
end

% Initialization
alpha = 2.5*N;
preGrad = zeros(N,no_dims);
epoch = 1;
vcc = 1;
cost = [];
len = 2;
while epoch <= T_epoch  && vcc > T_vcc
    Pgrad = zeros(N,no_dims);
    Qgrad = zeros(N,no_dims);
    sumQ = 0;
    KL = 0;
    % Compute gradient
    for i = 1:no_blocks
        idx = mark(i,1):mark(i,2);
        D = pdist2(Y(idx,:),Y).^2;
        Q1 = 1./(1+log(1+D));
        QQ1 = 1./(1+D);
        Pmat = -4*P(idx,:).*Q1.*QQ1;
        Qmat = -4*Q1.^2.*QQ1;
        len_blk = mark(i,2) - mark(i,1) + 1;
        idPQ = (len_blk*(mark(i,1)-1)+1):(len_blk+1):(len_blk*mark(i,2));
        Pmat(idPQ) = Pmat(idPQ) - sum(Pmat,2)';
        Qmat(idPQ) = Qmat(idPQ) - sum(Qmat,2)';
        Pgrad(idx,:) = Pmat*Y;
        Qgrad(idx,:) = Qmat*Y;
        sumQ = sumQ + sum(Q1(:));    
        KL = KL - sum(sum(P(idx,:).*log(Q1)));
    end
    % Update embedding Y
    Y = Y - alpha*(Pgrad - Qgrad/(sumQ-N) + (epoch-1)/(epoch+2)*preGrad);
    preGrad = Pgrad - Qgrad/(sumQ-N);
    % Compute KLD cost
    KL = KL + sum(P(:))*log(sumQ-N);
    cost = [cost;KL];
    % Update learning rate
    if(epoch > 1)&&(cost(end)-cost(end-1))>0
        alpha = alpha*0.99;
    end
    % Compute variation coefficient of the last three KLD costs
    if(epoch > len)
        vcc = (len./(len+1))*var(cost(end-len:end))./mean(cost(end-len:end));
    end
    epoch = epoch + 1;
end
disp([num2str(epoch-1),' epochs have been computed!']);
