function Y=findn(X,   varargin)
% Y=findn(X)
%
% for an array X, go through each column of X and get FIND( ) the set of 
% elements that are true. Return an array of integer indices.
% e.g. if size(X) = [ 10 11 12 ], then 
%         size(Y) = [  n 11 12 ]
% where n is the maximum number of items which are true in any column of X.
% note that the value is zero if no element is found.
% shortcut to get x(1). was very handy before 'find' accepted a count.
% returns [] if x is empty.
% 
% Y=findn(X, 'rows') : produce output as a matrix of subscripts, with one
%    column for each dimension of X, and a row for each found value. 
%    in other words, output in "long form"
% sgm 2011-2019


if strcmp(varargin{1}, 'rows')
  DIMENSIONAL = false;
  varargin(1)=[];
else
  DIMENSIONAL = true;
end

if DIMENSIONAL
  
  sx = size(X);
  X2 = reshape(X,size(X,1),[]); % make into matrix
  Y2=[];
  for i=1:size(X2,2)
    w = find(X2(:,i), varargin{:} );
    if ~isempty(w)
      if length(w)>1 || size(Y2,1)>1 % multiple items found ?
        if isempty(Y2),  Y2=w;
        else
          Y2 = nancat(2, Y2,w);
        end
      else % single items only
        Y2(1,i) = w;
      end
    else
      Y2(:,i) = 0;  % nothing found
    end
  end
  Y = reshape(Y2,[size(Y2,1),sx(2:end)] ); % restore to correct shape


else  % create output as?
  
  nd = ndims(X); 
  out = repmat({[]}, 1,nd); % create space for ind2sub outputs
  % run find to get linear indices, then convert to subscripts
  [out{:}]=ind2sub( size(X), find( real(X(:))>0 ) );
  % yield subscripts as a row.
  Y = [out{:}];
  
end

