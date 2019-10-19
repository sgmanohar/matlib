function nd=ndims1(X)
% function nd=ndims1(X)
% number of dimensions 
% just like built-in function ndims
%
% BUT: returns 
%   0 for scalars, 
%   1 for column vectors, 
%   2 for row-vectors or matrices, 
%   and higher values for other arrays.
%

    nd=ndims(X); 
    if     isscalar(X), nd=0; 
    elseif iscolumn(X), nd=1; 
    end
