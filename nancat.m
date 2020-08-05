function C=nancat(DIM, varargin)
%   X = nancat( DIM, x1, x2, x3... )
%         Concatenate arrays xi along dimension DIM, with nan padding.
% Works just like cat. Doesn't complain if x2 is longer or shorter than x1,
% it simply pads the one that is shorter with nans. Works with arbitrary
% number of dimensions.
%   X = nancat( [SOURCE_DIM, DEST_DIM], cell_array )
%         Concantenate elements from SOURCE_DIM in the cell array, along DEST_DIM. 
%         e.g.   nancat([1 2], { [1;2], [3;4] ; [5;6], 7 })     i.e. X{2,2}(2,1)
%         gives  { [1 5; 2 6], [3 7; 4 nan] }                   i.e. X{1,2}(2,2)
%   X = nancat( cell_array )
%         Concatenate the arrays in the cell-array, along the next
%         dimension (calculated using the size of the first cell)
%   X = nancat( struct_array ) : Concatenate each field using nancat.
%
% 'padvalue', VAL  (defaults to nan) 
%        VAL specifies a single value that is used to fill unused space. 
%        It should be the same class as the items being concatenated, 
%        and it can't be a vector or matrix.
% 'alignend' , true
%        aligns each item at the ends of its dimensions than
%        beginnings. (default = 0, align at beginning of each dimension)
% 2013 - now deals with cell arrays and character arrays.
%        the default padding for cells = { []  [] ... }, and chars = ' '.
%        If any concatenands are cells, the output is a cell.
%        If any of the concatenands are numeric, they are turned into cell
%        arrays with num2cell. 
%        -- if the first operand isn't a cell, but later ones are, we
%        default to padding with { [nan] }, rather than { [] }.
% 2015 - now has 'cell concatenation' function.
% 2016 - if the input is a structure array, then apply the concatenation to the
%        fields of the stucture array
% (c) sanjay manohar 2007

VERBOSE = true;
if nargin==0
    error('syntax: X=nancat(DIMENSION, X1, X2, ...)');
end
if(nargin==1) % just one argument? It has to be a cell-array.
  C=DIM; 
  if iscell(C) % if a cell array, and no dimensions provided, 
    firstdim = ndims1(C{1}) + 1;  % find next 'free' dimension, using first cell
    % get dimensions of cell array
    sz=size(C); if ndims1(C)==1, sz=sz(1); end 
    for j=1:length(sz); % for each dimension of the cell array, 
      % call nancat to move it into the cells
      C=nancat( [ j, firstdim-1+j ], C );
    end
    % now all dimensions should be transferred into the cell, and 
    % we should be left with a single cell. Return its contents.
    C=C{1}; 
    return
  else % if not a cell, then give an error
    error('syntax: X=nancat(DIMENSION, X1, X2, ...)');
  end
end
% at least two arguments provided:
C=varargin{1}; % a dimension is supplied?
if nargin==2 && length(DIM)==2 && iscell(C), C=cellcat(DIM,C); return ; end; % cell concatenation!
if nargin==2 && isstruct(C) && numel(C)>1 % is it a structure array? 
  % if so, devolve the whole thing for each field
  fn=fieldnames(C); 
  if ndims1(C)>1, error('to deal with n-dimensional struct arrays, use catStruct instead'); end
  for i = 1:length(fn) % for each field:
    y.(fn{i}) = nancat(DIM, C.(fn{i})); % call nancat with that field.
  end
  C=y; return;
end
if nargin==2,   return; end; % nothing to concatenate?  return
% default padding values
if      isnumeric(C) || islogical(C),   padvalue = nan;
elseif  iscell(C)                       padvalue = { [] };
elseif  ischar(C)                       padvalue = ' ';
% if it's a struct, this is a hack to get an empty structure with same fields: 
elseif  isstruct(C)                     padvalue = emptyStructLike(C);
else    error('nancat: unknown type, %s', class(C));
end

% check for options
alignend=0; remove=[];
for(i=1:length(varargin))
    if(strcmpi(varargin{i},'alignend'))
        alignend=varargin{i+1};
        remove=[remove i i+1];
    end
    if strcmpi(varargin{i},'padvalue')
      padvalue=varargin{i+1};
      remove=[remove i i+1];
    end
end;  varargin(remove)=[];

% create padding of a given size
makepadding = @(sz) repmat( padvalue, sz );


for(i=2:length(varargin)) % for each operand
    time_i = tic;
    v = varargin{i}; % concatenate v onto C
    if iscell(C) && isnumeric(v), v=num2cell(v); end % convert v to cell
    if iscell(v) && isnumeric(C), C=num2cell(C); padvalue=cell(padvalue); end % switch to cell mode
    sizeC=size(C); sizeV=size(v);
    if all(sizeC==0),
      C=v; warning('nancat:emptyskipped', 'first concatenand empty, so it was skipped.'); continue;
    end; % earlier operands were empty?
    if    (length(sizeC)>length(sizeV)) sizeV=[sizeV ones(1,length(sizeC)-length(sizeV))]; 
    elseif(length(sizeV)>length(sizeC)) sizeC=[sizeC ones(1,length(sizeV)-length(sizeC))]; 
    end % ensure same number of dimensions.
    uneqdim = find(~(sizeV==sizeC)); % for each of the unequal dimensions
    for(j=1:length(uneqdim))
        sizeC=size(C); sizeV=size(v); % Sgm 2019 - sizes may have changed!
        if length(sizeV)<length(sizeC) % if there are extra dimensions in the current accumulated value
          sizeV=[sizeV ones(1,length(sizeC)-length(sizeV))]; % add on the extra dimensions as size 1
        end
        if length(sizeV)>length(sizeC) % if there are fewer dims in the accumulator:
          sizeC=[sizeC ones(1,length(sizeV)-length(sizeC))]; % add them on as ones
        end
        if uneqdim(j)==DIM; continue; end % (except for the one you're catting along)
        d=uneqdim(j);
        if(sizeC(d)>sizeV(d)) % V is too small: add nans to d if it's shorter
            diff= sizeC(d)-sizeV(d); % difference in size along dimension d
            addsize=sizeV;  
            addsize(d)=diff; 
            addsize(addsize==0)=1; % sgm 2019
            if(alignend)   v=cat(d, makepadding( addsize ), v);
            else           v=cat(d, v, makepadding( addsize ));
            end
        elseif(sizeV(d)>sizeC(d)) % C is too small
            diff= sizeV(d)-sizeC(d); % difference in size along dimension d
            addsize=size(C); addsize(d)=diff;
            if(alignend)    C=cat(d, makepadding( addsize ), C);
            else            C=cat(d, C, makepadding( addsize ));
            end
        end
    end
    C=cat(DIM,C,v);
    time_i = toc(time_i);
    if time_i*nargin > 30 && VERBOSE % will the operation take more than 30 seconds?
      fprintf('concatenating item %g\n', i);
    end
end

function out = cellcat( DIM, cell )
% DIM is a pair, specifying 
%  DIM(1): the cell dimension to cat, and 
%  DIM(2): which dimension ot the output arrays to put it on.
% note that if a cell is empty, it won't give rise to nans.
sz = size(cell);
% arrange so that cat dimension is first
permorder = 1:length(sz); permorder = [DIM(1), permorder(permorder~=DIM(1))];
C=permute(cell, permorder); 
% now collapse all other dimensions
C=reshape(C, sz(DIM(1)),[]); 
for i=1:size(C,2) % for all other dimensions
  if any(cellfun(@isempty,C(:,i))), warning('nancat is leaving out empty cells'); end
  out{i} = nancat( DIM(2), C{:,i} ); 
end
% then make shape of new cell the same as original cell, except collapsed
% on the required the DIM(1) dimension.
endsz = sz; endsz(DIM(1))=1; 
out = reshape(out, endsz); 
