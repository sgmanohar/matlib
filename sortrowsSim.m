function [Xsorted, order] = sortSim(X,varargin)
%   [Xsorted, order] = sortSim( X, ...)
% X is a data matrix: X ( samples, dimensions )
% this computes a distance between each pair samples, 
% and performs hierarchical clustering to sort the samples 
% in order of similarity.
% additional arguments are passed to linkage().
% sgm 2020

D = pdist(X, 'mahalanobis'); % pairwise distances
Z = linkage(D,varargin{:})
order = optimalleaforder(Z,D);
Xsorted = X(order,:);
