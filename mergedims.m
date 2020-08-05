function Y = mergedims( X, dims )
% Y = mergedims( X, dims )
%   merge a list of dimensions of an n-D array into a single dimension
% eg: X = rand(  2,3,5,7,11  )      % create a 5-dimensional array
%     Y = mergedims( X, [2,4] )     % merge two dimensions
%     size(Y)
%           2,  21,  5, 11          % results in a 4-dimensional array
%     in otherwords in this case it would be equivalent to: 
%        reshape( permute( X, [1,2,4,3,5] ), 2, 3*7, 5, 11 ) 


% place the new collapsed dimension at the position of dims(1)
D = [ 1:dims(1)-1,  dims ];
% find which dimensions will be trailing - the ones not included so far
trailing = 1:ndims1(X);

dims(dims>ndims1(X)) = []; % ignore dimensions higher than number available

% trailing( any(trailing==D') ) = [];  % new version
trailing( any( bsxfun( @eq, trailing, D' ) ) ) = [];
D = [ D trailing ];                   % put them on the end
sz = size(X); % size of original
Y = permute(X, D);
% compute new size
sz( dims(1) ) = prod( sz(dims) );     % new merged dimension = product
sz( dims(2:end) ) = [];               % squeeze out old merged dimensions
Y = reshape(Y, sz);


