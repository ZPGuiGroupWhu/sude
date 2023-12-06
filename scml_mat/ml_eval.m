function [CC, knnACC, svmACC, clusACC] = ml_eval(X, Y, ref)
% This function evaluates the embedding quality using congruence 
% coefficient (CC), knnACC, svmACC, clusACC metrics.

[n, ~] = size(X);
%% Global structure
X_dis = pdist(X);
Y_dis = pdist(Y);
CC = sum(X_dis.*Y_dis)/sqrt(sum(X_dis.^2)*sum(Y_dis.^2));
%% kNN accuracy
train_id = 1:4:n;
test_id = setdiff(1:n, train_id);
train_data = Y(train_id,:);
train_labels = ref(train_id,:); 
test_data = Y(test_id,:);
k = 5;
knn_classifier = fitcknn(train_data, train_labels, 'NumNeighbors', k);
predicted_labels = predict(knn_classifier, test_data);
[knnACC] = getACC(ref(test_id,:), predicted_labels);
%% SVM accuracy
svm_classifier = fitcecoc(train_data, train_labels);
predicted_labels = predict(svm_classifier, test_data);
[svmACC] = getACC(ref(test_id,:), predicted_labels);
%% Clustering accuracy
clus = kmeans(Y,length(unique(ref)),"MaxIter",200,"Replicates",5);
[clusACC] = getACC(ref, clus);
end
