function y = pyth (x,dim)
% PYTH (X, DIM)
% Applies pythagorus' theorem along final or specified dimension

x = sq(x); if size(x,2)==1, x=x'; end
if nargin<2, dim=ndims(x); end

y = sqrt(sum(x.^2,dim));

