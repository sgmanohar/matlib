function K = gaussKern( sz, falloff )
% K = gaussKern( sz , [falloff] )
% create an n-dimensional gaussian kernel.
% useful as an input to convn.
% example:    sz = [ 5 ]    : create a 11 x 1 vector 
%             sz = [ 2 2 2 ]: create a 5 x 5 x 5 array.
% defaults:   falloff = 0.8, higher values mean steeper falloff
%             sz      = [ 3 ]
% sgm 2019
s = [];
if ~exist('falloff','var'), 
  falloff = 0.8; 
end
if ~exist('sz','var'), 
  sz = 3; 
end
for i=1:length(sz) % for each dimension of the kernel
  Ni = sz(i);
  si = ( [-Ni:Ni]' ./ (falloff * Ni) ); % distance on this dimension
  % now move it to dimension i
  newdims = [ [2:i] 1 ]; % new order for dimensions
  % take care of a matlab quirk which doesn't like 1D variables
  if length(newdims)==1, newdims = [newdims 2]; end 
  si = permute( si,  newdims); 
  % take care of matlab quirk where bsxfun doesn't expand an empty array
  if isempty(s) 
    s = si.^2;
  else
    s  = bsxfun( @plus, s , si.^2 );  % add it, expanding singleton dimensions
  end
end
K = exp( - s / length(sz) );
