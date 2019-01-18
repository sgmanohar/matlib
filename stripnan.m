function Y = stripnan(X, DIM)
% Y = stripnan(X, DIM)
% takes an array removes the nans 'within' moving values 'up' 
% each column. 
%
% X can be a structure where each field has the same size. In this case,
% operate across all fields together, removing positions with nan for any
% field.
%
% sgm 2018

if ~exist('DIM','var'), DIM=1; end

if isstruct(X)
  % deal with structures where each field is the same shape.
  % apply stripnan across all fields, treating positions with nans 
  % in any field as bad
  fn=fieldnames(X);
  % select bad trials, and go across all fieldss
  b=false(size(X.(fn{1})));
  for i=1:length(fn)
    b = b | isnan(X.(fn{i}));
  end % b now indicates if a value was 
  for i=1:length(fn)
    xi = X.(fn{i});  % get one field at a time
    xi(b)=nan;       % remove bad positions
    Y.(fn{i}) = stripnan(xi,DIM);
  end
  return
end

X=shiftdim(X, DIM-1);
sz=size(X);
X=reshape(X,sz(1),[]);
for i=1:size(X,2)
  % select non-nan values
  Y{i} = X( ~isnan(X(:,i)), i );
end
% combine as columns
Y=nancat(2,Y{:});
% restructure later dimensions
sz2=num2cell(sz);
Y=reshape(Y,[],sz2{2:end});
% re-permute if operating on other dimension
Y=shiftdim(Y,ndims(Y)-DIM+1);



