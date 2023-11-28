function id = pps(knn, rnn, order)
%   Plum Pudding Sampling (PPS)
%   This function returns the ID of landmarks. 
%
%   Parameters are:
%
%   'knn'   - A N by k matrix. Each row represents the KNN of each point. 
%   'rnn'   - A vector of length N. Each row represents the RNN of each point. 
%   'order' - A positive integer specifying the order of KNN. Once a landmark is selected, its 
%             KNN will be removed from the queue. PPS supports the removal of multi-order KNN. 

id = [];
[~, id_sort] = sort(rnn,'descend');
while ~isempty(id_sort)
    id = [id;id_sort(1)];
    rm_pts = [id_sort(1)];
    for i=1:order
        temp_pts = knn(rm_pts,:);
        rm_pts = [rm_pts;temp_pts(:)];
    end
    rm_pts = unique(rm_pts);
    [~,rm_id] = ismember(rm_pts,id_sort);
    id_sort(nonzeros(rm_id)) = [];    
end