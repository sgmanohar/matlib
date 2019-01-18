function MATRIX = nanassign(MATRIX, INDICES, VECTOR)
% function MATRIX = nanassign(MATRIX, INDICES, VECTOR)
%  
%  subscripted assignment, with nan-padding
%  set the slice of MATRIX corresponding to the indices which are nan. 
%  e.g. if INDICES is [2 5 nan 7], then do
%
%    MATRIX(2,5,:,7) = VECTOR
% 
%  extending MATRIX if vector is too long,
%  and extending VECTOR if it is too short
% 
%  nanassign works by padding the matrix/vector as needed, then calling 
%  'subsasgn' with the appropriate structure
% 
% sgm

if ~iscell(INDICES)
  col=find(isnan(INDICES));
else
  col = find(strcmp(INDICES,':')); 
end
if(length(col)>1) error('only one nan allowed at present!');end
if(sum(size(VECTOR)>1)>1) error('vector must be 1-dimensional');end
%if(length(size(MATRIX))<length(INDICES)) error('please make MATRIX have as many dimensions as INDEX');end;
sm=size(MATRIX,col); sv=length(VECTOR);
if(sv>sm)
  addsize = size(MATRIX);
  addsize(col) = sv-sm;
  addsize(addsize==0)=1; % if the matrix isempty, make sure it has some elements.
  MATRIX=cat(col, MATRIX, nan*ones(addsize));
elseif(sm>sv)
  dim=find(size(VECTOR)>1);
  if(isempty(dim)) dim=1;end % if single elem, doesn't matter which dir!
  % pad with nans along the dimension of the vector
  reorder_dims =  [2:dim 1];                           % move first dimension to end
  if numel(reorder_dims)==1, reorder_dims = [1 2]; end % unless it is the only dimension
  VECTOR=cat(dim, VECTOR, permute(nan*ones(sm-sv,1), reorder_dims) );
end

if ~iscell(INDICES)
  for(i=1:length(INDICES))
    if(isnan(INDICES(i)))
      subs{i}=':';
    else
      subs{i}=INDICES(i);
    end
  end
else 
  subs = INDICES;
end
substruc.type='()';
substruc.subs=subs;
MATRIX = subsasgn(MATRIX, substruc , VECTOR);