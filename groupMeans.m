function [mu, resid] = groupMeans( X, DIM, categ , meanfunc )
%   mu = groupMeans( X, CATEG )         
%   mu = groupMeans( X, DIM, CATEG, [ meanfunc ] )
%   [mu resid] = groupMeans( X, DIM, CATEG, [ meanfunc ] )
%
% MEAN of the matrix X along dimension specified, 
% but grouping the other dimensions according to the categories in 'categ'.
% 
% Returns: a matrix with same size as X, except for the specified dimension -- 
% Along that dimension, there is one slice per unique level of 'categ'.
% So, size( mu ) = size( X ) for all dimensions except DIM.
% 
% RESID returns the residuals after subtracting the relevant means from X.
% 
% CATEG can either be 
%  1) an indicator variable, that is, a vector of numeric
%     values same length as the specified dimension of X. Each unique
%     value of CATEG generates a set of means, on dimension DIM.
%     Therefore, size( mu, DIM ) = length( unique( CATEG ) )
% 
%     example:
%       X = [1 2 3; 4 5 6; 7 8 9; 10 11 12]
%       C = [ 1;     1;     2;      2 ]; % each row of X has a category
%       groupMeans( X, 1 ,C )  % mean along vertical dimension, but group
%                              % rows according to C
%     ans = 
%       2.5   3.5   4.5
%       8.5   9.5  10.5   
%  
%     i.e. the first row is        mean( X( C==1, : ) , 1 )
%          and the second row is   mean( X( C==2, : ) , 1 )
%
%  2) an array of logicals, with the same number of rows as the specified
%     dimension of X. Each column of CATEG generates a set of means on  
%     dimension DIM.
%     Therefore, size( mu, DIM ) = size( CATEG, 2 ).
% 
%     example:
%       X = [1 2 3; 4 5 6; 7 8 9; 10 11 12]
%       C = logical([ 1 0; 1 0; 0 1; 0 1 ])  % this selects the same as the
%                                            % first example, using logicals 
%       groupMeans( X, 1, C )      % each column of C selects a 
%                                  % subset of the rows of X
%       ans = 
%         2.5   3.5   4.5
%         8.5   9.5  10.5
%
%  3) an array of categories, the same size as X. 
%     Elements of X are selected according to the unique values of CATEG.
%     For each category and each column along dimension DIM, those values
%     are passed to meanfunc, and entered in the output matrix:
% 
%     e.g. for DIM==1  :
%      MU(i,j,k)  =  mean( CATEG(:,j,k) == i ,j,k )
%     where i refer to the unique values of CATEG
%
% MEANFUNC is the function used to calculate the mean.
%     it defaults to @NANMEAN if you have it on the path, otherwise it uses
%     @MEAN. The function takes the parameters MEANFUNC ( DATA, DIM ).
%     Can also be one of the following strings:
%       'std' :  standard deviation
%       'cell':  create a cell array of the grouped values, rather than
%                taking the means. In other words, just groups the values
%                into a one-dimensional cell array, each containing an array 
%                similar to X (except smaller on the specified dimension)
%       'dim' :  add an extra dimension to X, at the end. The specified
%                dimension is the categories, and the new trailing
%                dimension is the different values in that category.


% check for nanmean
saveix_extend = false; % extend by one dimension, by expanding the cell afterwards
meanIsCell    = false; % concatenate result as cells rather than matrix


if isscalar(X) && isinteger(X) % did you swap the Dimension and Data parameters?
  tmp=X;X=DIM;DIM=tmp;
end
if ~exist('categ','var') % DIM must be Categories! 
  categ=DIM; DIM=1;      % assume dimension 1.
end
if isscalar(categ) && isinteger(categ) % did you swap the Dimension and Category?
  tmp=categ; categ=DIM; DIM=tmp;
end

% deal with structures: distribute the operation over fields.
if isstruct(X) && numel(X)==1
  fn=fieldnames(X); 
  %DIM=find( size(X.(fn{1})) == length(categ) );
  for i=1:length(fn) % must supply meanfunc for this
    M.(fn{i}) = groupMeans( X.(fn{i}), DIM, categ, meanfunc );
  end
  mu = M;
  return
end

if iscell(categ) && ~iscell(X)
  mu = X; C=categ;
  for i=1:length(categ)
    mu=groupMeans(mu, DIM+i-1, C{i}, meanfunc);
    for j=i+1:length(categ)
      C{j} = groupMeans( C{j}, DIM+i-1, C{i}, 'dim' ); 
    end
  end
  return 
end

if ~exist('meanfunc','var') || isempty(meanfunc)
  if exist('nanmean','file')
    meanfunc = @nanmean;
  else
    meanfunc = @mean;
  end
elseif ischar(meanfunc) && strcmpi(meanfunc,'std') % shortcut for standard deviation on dimension DIM
  meanfunc=@(x,y)std(x,[],y);
elseif ischar(meanfunc) && strcmpi(meanfunc,'cell') % return a grouped cell-array. 
  meanfunc=@(x,y){x};
  meanIsCell=true;
elseif ischar(meanfunc) && strcmpi(meanfunc,'dim') % "dim": create new dimension
  % note that tries to add a new dimension at the end of the nd-array
  meanfunc=@(x,y){x};
  %permute(x, [ndims1(x)+1,1:ndims1(x)]);
  saveix_extend = true;
  meanIsCell=true;
end


% check supplied parameters are valid
if ~isvector(categ) && ~islogical(categ) 
  % error('groupMeans:params', 'multidimensional categories not yet supported!'); 
end  

% create a cell array like {':',':',':',':',':'}
all_indices = repmat({':'},1,ndims(X));

if nargout>1, resid = X; end

if (isnumeric(categ) || islogical(categ)) && isvector(categ)
  if length(categ) ~= size(X,DIM), error('groupMeans:params', 'categories must be same length as the specified dimension of the data'); end
  % find unique levels of CATEG, ignorning NaN
  u = unique(categ);
  u(isnan(u))=[];
  % if the categories start at 1 and are all integers,
  if isnumeric(u) && u(1)==1 && all(u == floor(u)) 
    INTEGER_CATEGS = true; 
    if any(diff(u)>1)
      warning('using categories as index; some results might be blank');
    end
  else
    INTEGER_CATEGS = false;
  end
  
  % create a row in MU for each level of CATEG
  for i=1:length(u)
    searchix      = all_indices;
    searchix{DIM} = categ==u(i);  % take mean of only elements that match category on dimension DIM
    saveix        = all_indices;
    if INTEGER_CATEGS % use category value as the index for storing
      saveix{DIM} = u(i);
    else
      saveix{DIM} = i;            % store mean to all indices on the current slice of DIM
    end
    % now calculate mean and assign it to the correct part of mu.
    % original: 
    %mu(saveix{:}) = meanfunc(X(searchix{:}),DIM); 
    % replaced by a 'nanassign' in case they are different sizes
    if meanIsCell
      mu(saveix{:}) = meanfunc(X(searchix{:}),DIM); 
    else
      if ~exist('mu','var'), mu=[]; end
      mu = nanassign( mu, saveix, meanfunc(X(searchix{:}),DIM) );
    end
    if nargout==2
      resid(searchix{:}) = bsxfun(@minus, resid(searchix{:}), mu(saveix{:}) );
    end
  end
elseif islogical(categ)
  for i=1:size(categ,2) % for each column of categ
    searchix      = all_indices;
    searchix{DIM} = categ(:,i);
    saveix        = all_indices;
    saveix{DIM}   = i;
    mu(saveix{:}) = meanfunc(X(searchix{:}),DIM);
    if nargout>1
      resid(searchix{:}) = bsxfun(@minus, resid(searchix{:}), mu(saveix{:})); 
    end
  end
elseif size(categ)==size(X) % numeric AND multidimensional category!
  % find unique levels of CATEG, ignorning NaN
  u = unique(categ);
  u(isnan(u))=[];
  
  % new order of dimensions: bring DIM to front.
  oldsz  = size(X);
  neword = 1:length(size(X)); neword(DIM)=[]; neword=[DIM, neword];
  invord = inversePermutation([1:length(size(X))], neword);
  newX   = reshape( permute(X, neword),    oldsz(DIM), []) ; % bring DIM to 1 and flatten other dims
  categ  = reshape( permute(categ, neword),oldsz(DIM), []) ;
  % create a row in MU for each level of CATEG
  for i=1:length(u) % for each category
    f = categ==u(i); % a filter for this category
    for j=1:size(newX,2); % for each column
      col=newX( f(:,j) ,j); % get the items matching category
      mu( i,j )   = meanfunc(col); % calculate the mean, store in appropriate location of mu
      if nargout>1 % do we need to calculate residuals too?
        resid( f(:,j),j )  = col-mu(i,j);
      end
    end
  end
  newsz=oldsz; newsz(DIM)=length(u);
  %sgm 2017
  newsz=oldsz; newsz(DIM)=[]; newsz=[length(u) newsz]; % add dimension on front
  mu=permute(reshape(mu,newsz), invord); 
  %mu = reshape( permute(mu, invord), newsz ); % OLD - collapses dimensions
  if nargout>1
    resid = reshape( permute( resid, invord ), oldsz );
  end
else
  error('Categories does not match data');
end
if saveix_extend
  dims = ndims1(mu);
  dim2 = ndims1(mu{1});
  for i=1:dims
    mu=nancat( [i,i+dim2], mu );
  end
  mu=mu{1};
  % now bring the "new dimension" (items within a condition)  to
  % the end
  mu=permute(mu, [2:dims+dim2, 1]);
  if 1 % should we squeeze it too?
    mu=sq(mu);
  end
end
  