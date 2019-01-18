function Y = apply(f,X,varargin)
% Y = apply(f,X)
% Y = apply(DIM, f,X)
% Loop over an array, and apply a function to each part.
% Like arrayfun but allows control over dimension and windowing. 
%
% f is a function that takes a vector, matrix, or nd-array.
% DIM = dimension(s) on which to apply the function.
% e.g. if DIM = 1,     then apply the function to each column if X.
%   or if DIM = [1,2], then the first two dimensions of X are treated as a
%                      matrix, and each matrix is passed to the function.
% note: the size of the result of f must be the consistent, and have the 
% same or fewer number of dimensions as the input. 
%
% options:
% 'conv', C: convolve the function with a given window size, over the
%   dimensions of interest. must have same number of elements as DIM. 
%   e.g. apply([1,2],@norm,X,'conv',[2,2]) - get the norm of the 2x2
%   matrices, scanning over dimensions 1 and 2 of X. 
%
%   The output will be smaller than X, on the indicated dimensions, similar
%   to using the 'conv' function with the option 'valid'. 
%   In this case, since norm returns a scalar, 
%   the result will also have a singleton dimension at dimension 2.
%   if X = 5x5x5, then the result will be 4x1x5
%
% 'comb', D: take all pairwise combinations of 'columns' on dimension D, 
%   and passes both of them to f. Note that in this case, f must take two
%   arrays, each of size [DIMS excluding dimension D], and return a scalar.
%   D must be an integer indicating which index of DIM (i.e. from 1 to
%   length(DIM)), that combinations are taken over.
%   e.g. apply( [1 2], @(x,y)corr(x,y), X, 'COMB',2 )
%   calls corr on each pair of columns of X.


% which dimension(s) to pass to the function
DIM=1;

% if the first parameter is an integer or vector of integers
% and the second parameter is a function
if isnumeric(f) && numel(f)<5 && all(floor(f)==f) ...
    && strcmp(class(X),'function_handle')
  DIM=f; f=X; X=varargin{1}; varargin=varargin(2:end);
end

i=find(strcmp(varargin,'conv'));
if numel(i)==1
  CONV = varargin{i+1};
  varargin(i:i+1)=[];
  if length(CONV)~=length(DIM), error('conv must be same size as DIM'); end
else
  CONV = false;
end
i=find(strcmp(varargin,'comb'));
if numel(i)==1
  if CONV, error('conv and comb cannot be used together'); end
  COMB = varargin{i+1};
  varargin(i:i+1)=[];
  if ~isscalar(COMB) || COMB>length(DIM), error('comb must indicate which dimension of DIM to take combinations'); end
else
  COMB = false;
end

% original size
S = size(X); 
% new order of dimensions
neword = 1:length(S); neword(DIM)=[]; neword=[DIM neword];
% how to get back to old order?
oldord = inversePermutation(1:length(S), neword);
% reorder dims
xp = permute(X, neword);
% flatten on the irrelevant dimensions
Sa = num2cell(S);
x = reshape(xp,Sa{DIM},[]);
% a set of colons, for passing the correct dimensions of x to function
colons = repmat({':'}, length(DIM),1);
for i=1:size(x,length(DIM)+1)
  % extract a subarray for only the dimensions DIM
  M = x(colons{:},i);
  if ~CONV
    if ~COMB
      y(colons{:},i) = f(M);
    else % all combinations for a particular dimension
      allsel = repmat({':'},length(DIM),1);
      for j=1:S(DIM(COMB))
        for k=1:S(DIM(COMB))
          s1 = allsel; s1{COMB} = j;
          s2 = allsel; s2{COMB} = k;
          y(j,k,i) = f( M(s1{:}), M(s2{:}) );
        end
      end
    end
  else % convolution along requested dimensions
    if length(DIM)==1 % M is a vector
      for j=1:S(DIM)-CONV
        y(j,i) = f(M(j:j+CONV));
      end
    elseif length(DIM)==2 % M is a matrix
      for j=1:S(DIM(1))-CONV(1)
        for k=1:S(DIM(2))-CONV(2)
          y(j,k,i) = f(M(j:j+CONV(1), k:k+CONV(2)));
        end
      end
    elseif length(DIM)==3
      for j=1:S(DIM(1))-CONV(1)
        for k=1:S(DIM(2))-CONV(2)
          for l=1:S(DIM(3))-CONV(3)
            y(j,k,l,i) = f(M(j:j+CONV(1), k:k+CONV(2), l:l+CONV(3)));
          end
        end
      end
    end
  end
end

% reshape the irrelevant dimensions back
Sy   = size(y);
newS = S(neword); newS(1:length(DIM)) = Sy(1:length(DIM));
y = reshape(y,newS); 
% rearrange dimensions back
Y = permute(y, oldord); 
