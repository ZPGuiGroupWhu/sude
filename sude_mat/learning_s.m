function [Y, k2] = learning_s (X_samp, k1, get_knn, rnn, id_samp, no_dims, initialize, agg_coef, T_epoch)
% This function returns representation of the landmarks in the lower-dimensional space and the number of 
% nearest neighbors of landmarks. It computes the gradient using the entire probability matrix P and Q.

% Obtain size and dimension of landmarks
[N, dim] = size(X_samp);

% Compute the number of nearest neighbors of landmarks adaptively
if (N < 9)
    k2 = N;
else
    if(N > 1000)
        k2 = ceil(log(N)/log(2)) + 18;
    elseif(N > 50)
        k2 = ceil(0.02*N) + 8;
    else
        k2 = 9;
    end
end

% Compute high-dimensional probability matrix P
if(k1 > 0)
    % Compute SNN matrix of landmarks
    SNN = zeros(N, N);
    knn_rnn_mat = rnn(get_knn(id_samp,:));   
    for i = 1:N
        snn_id = ismember(get_knn(id_samp,:),get_knn(id_samp(i),:));
        nn_id = find(max(snn_id,[],2)==1);
        SNN(i,nn_id) = sum(knn_rnn_mat(nn_id,:).*snn_id(nn_id,:),2);
        SNN(i,:) = SNN(i,:)./max(max(SNN(i,:)),realmin);
    end
    % Compute the modified distance matrix
    Dis = (1-SNN).^agg_coef.*pdist2(X_samp,X_samp);
    P = zeros(N,N);
    [sort_dis, idx] = sort(Dis, 2);
    for i = 1:N
        P(i,idx(i,1:k2)) = exp(-0.5*sort_dis(i,1:k2).^2./max(mean(sort_dis(i,1:k2)).^2,realmin));
    end
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
P = sparse((P + P'))/2;
clear xx get_knn rnn id_samp samp_knn samp_dis snn_id SNN knn_rnn_mat

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
clear Dis Dg L

% Initialization
max_alpha = 2.5*N;
min_alpha = 2*N;
warm_step = 10;
preGrad = zeros(N,no_dims);
epoch = 1;
while epoch <= T_epoch
    % Update learning rate
    if(epoch <= warm_step)
        alpha = max_alpha;
    else
        alpha = min_alpha + 0.5*(max_alpha-min_alpha)*(1+cos(pi*((epoch-warm_step)/(T_epoch-warm_step))));
    end
    % Compute low-dimensional probability matrix Q
    D = pdist2(Y,Y).^2;
    Q1 = 1./(1+log(1+D));
    QQ1 = 1./(1+D);
    Q = Q1/(sum(Q1(:))-N); 
    % Compute gradient
    ProMatY = 4*(P-Q).*Q1.*QQ1;
    grad = (diag(sum(ProMatY))-ProMatY)*Y; 
    % Update embedding Y
    Y = Y - alpha*(grad+(epoch-1)./(epoch+2)*preGrad);
    preGrad = grad;
    epoch = epoch + 1;
end
disp([num2str(epoch-1),' epochs have been computed!']);