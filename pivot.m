function [Y uc count] = pivot(X, varargin)
% X is a set of categories, and the last column of X is a value.
% X =  [ categ1, categ2, ..., value ]
% e.g. 
%   X =  [ 1 1   0.51 
%          1 1   0.23
%          1 2   0.36
%          2 2   0.62 ]     etc
%
% The categories are then used to index an n-dimensional array of the form 
%   Y ( categ1, categ2, ..., index ) = value
% so, if X is a [ N x 5 ] matrix, there are 4 categories, and so
%        Y will be a [ a x b x c x d x num ] array
% where num is the number of elements in the largest category.
% The other categories will be nan-padded.
%
%   [Y uc count] = pivot(X)
% will also give the ordered list of unique categories, and the number of 
% items found for each category combination.
%   UC{i}               is the list of categories for column i of X.
%   count(i,j,k,...)    is the number of items found where X(:,1)==i and
%                       X(:,2)==j ...
% 
% Y = pivot(X,'front',1) : brings the trailing 'index' dimension to the front, 
%   i.e. dimensions 2 upwards of Y correspond to the categories. This form is
%   useful for passing to 'mean' etc
% sgm 2015

[  FRONT    ]=parsepvpairs( ...
  {'front'  }, ...
  {false     }, varargin{:});


%%%%%%
% POOR-MAN'S VERSION - by me!
if 1
% unfortunately this is a very slow algorithm, looping through each value 
% separately.
    nc=size(X,2)-1;
    for i=1:nc % for each category
      uc{i} = unique(X(:,i)); % labels in each category
      uc{i}(isnan(uc{i}))=[]; % ignore NaNs in the categories
      lc(i) = length(uc{i});  % number of labels in each category
    end
    if numel(lc)==1 lc=[lc 1];end % prevent 'zeros' from creating an extra dimension if only 1D needed
    count=zeros(lc); % create matrix of counts
    for i=1:size(X,1) % for each data point
      if isnan(X(i,end)), continue; end % ignore NaNs in the data column
      catgs = X(i,1:end-1); % its category set, as a cell array
      for j=1:length(catgs) % for each category
        catgs(j)=find(catgs(j)==uc{j}); % convert the label to index
      end % now catgs is a list of indices
      where = num2cell(catgs); % indexing for the category
      % increment counter for this category combination
      count(where{:})=count(where{:})+1; 
      where{end+1} = count(where{:});  % where to store the value
      Y(where{:}) = X(i,end); % store the value in the appropriate place
    end
    % now fill remaining space with nans
    sz=size(Y); % output size
    for i=1:numel(count) % for each category combination
      [where{:}] = ind2sub( size(count), i); % find the index
      for j=count(where{:})+1:sz(end) % elements which should not exist
        where_blank = where;
        where_blank{end}=j;
        Y(where_blank{:})=nan;
      end
    end
    
else

  % a MUCH BETTER algorithm from luiz mendo!
  % but assumes categories are integers. Note the beautiful syntax.
  
  N = size(X,1); %// number of values
  [~, ~, label] = unique(X(:,1:end-1),'rows'); %// unique labels for indices
  cumLabel = cumsum(sparse(1:N, label, 1),1); %// used for generating a cumulative count
  %// for each label. The trick here is to separate each label in a different column
  lastInd = full(cumLabel((1:N).'+(label-1)*N)); %'// pick appropriate values from
  %// cumLabel to generate the cumulative count, which will be used as last index
  %// for the result array
  sizeY = [max(X(:,1:end-1),[],1) max(lastInd)]; %// size of result
  Y = NaN(sizeY); %// initiallize result with NaNs
  ind = mat2cell([X(:,1:end-1) lastInd], ones(1,N)); %// needed for comma-separated list
  Y(sub2ind(sizeY, ind{:})) = X(:,end); %// linear indexing of values into Y
end


if FRONT % move the 'index' dimension to the front - good for taking means.
  sz=size(Y);
  Y=permute(Y, [length(sz)  1:length(sz)-1]);
end