# scML: scalable manifold learning by uniform landmark sampling and constrained locally linear embedding
We propose a scalable manifold learning (scML) method that can cope with large-scale and high-dimensional data in an efficient manner. It starts by seeking a set of landmarks to construct the low-dimensional skeleton of the entire data, and then incorporates the non-landmarks into this skeleton based on the constrained locally linear embedding. 

![image](https://github.com/ZPGuiGroupWhu/scml/blob/main/github.png)

# How To Run

```ruby
function [Y, id_samp] = scml(X,varargin)
%   This function returns representation of the N by D matrix X in the lower-dimensional space and 
%   the ID of landmarks sampled by PPS. Each row in X represents an observation.
% 
%   Parameters are: 
%
%   'NumDimensions'- A positive integer specifying the number of dimension of the representation Y. 
%                    Default: 2
%   'NumNeighbors' - A non-negative integer specifying the number of KNN for PPS to sample landmarks.
%                    It must be smaller than N.
%                    Default: 10
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
%                    Default: 500
%   'TolVcc'       - Termination tolerance for variation coefficient of the last three KLD costs. 
%                    Default: 1e-7
```
