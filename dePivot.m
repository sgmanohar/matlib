function Y = dePivot(X, varargin)
% Y = dePivot ( X )
% flatten an N-dimensional array 
% think of the values in X as a single column of length PROD(SIZE(X))
% 
% return a matrix where the first N columns are categorical indices
% 1..size(X)(i) indicating  where the item was
% 
% imagine that you start with a pivot table, this returns you to the original
% data.
%
% parameters: 'extracategs'  : a matrix (or cell array of matrices) the 
%                              same size as X, specifying extra categories 
%                              for the items of X, which are just expanded
%                              into a new column.
% sgm 2015

[ extraCategs  ]  =parsepvpairs( ...
  {'extraCategs'},...
  {[]}, varargin{:});

    sz=size(X); Y=[];
    for i=1:length(sz) % for each dimension of X
      szi=sz; szi(i)=1; 
      ord=1:length(sz); ord(i)=1; ord(1)=i; % how to permute a column to the dimension of interest
      categ{i}=repmat(permute([1:sz(i)]',ord), szi); % create a matrix the same size as X, with the category value on dimension i
      Y=[Y categ{i}(:)]; % create a matrix indicating the category on each dimension-
    end
    if ~isempty(extraCategs)
      if isnumeric(extraCategs) % scalar, matrix or array
        extraCategs = bsxfun(@times,ones(size(X)),extraCategs); % make it the same size as X
        Y=[Y extraCategs(:)]; % add it as a column
      elseif iscell(extraCategs) % a cell array of matrices
        for i=1:length(extraCategs)
          Y=[Y extraCategs{i}(:)];
        end
      else error('unrecognised category parameter');
      end
    end
    Y=[Y X(:)];
    Y(any(isnan(Y),2),:)=[]; % remove rows where X is nan
