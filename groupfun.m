function [varargout] = groupfun( F, X, DIM, categ )
%   [ a b ... ] = groupfun( @ F(X), X, DIM, CATEG )
% 
% send groups of X to function F.
% F should be a function that collapses across dimension DIM.
% e.g. for DIM==2, you could supply   @(x)mean(x,2)
%
% the matrix X is grouped along the dimension specified, 
% according to the categories in 'categ'.
%   e.g. if size( X ) = [ 5  6  7 ],  and   DIM == 1
%        then F is applied to  matrices Y, with size( Y ) = [ N  6  7 ],
%        depending on the groups in categ
% 
% Along dimension DIM, there is F(X) per unique level of 'categ',
% given by a set of 'slices' of X perpendicular to DIM.
% So, size( Y ) = size( X ) for all dimensions except DIM.
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
%       groupfun( @mean, X, 1 ,C )  % mean along vertical dimension, but group
%                                   % rows according to C
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
%       groupfun( @mean, X, 1, C )      % each column of C selects a 
%                                  % subset of the rows of X
%       ans = 
%         2.5   3.5   4.5
%         8.5   9.5  10.5
%
% sgm 2012


% check supplied parameters are valid
if ~isvector(categ) && ~islogical(categ) error('groupMeans:multidimensional categories not yet supported!'); end
if length(categ) ~= size(X,DIM), error('groupMeans:categories must be same length as the specified dimension of the data'); end

% create a cell array like {':',':',':',':',':'}
all_indices = repmat({':'},1,ndims(X));

if isnumeric(categ)
  % find unique levels of CATEG, ignorning NaN
  u = unique(categ);
  u(isnan(u))=[];
  
  % create a row in MU for each level of CATEG
  for i=1:length(u)
    searchix      = all_indices;
    searchix{DIM} = categ==u(i);  % take mean of only elements that match category on dimension DIM
    saveix        = all_indices;
    saveix{DIM}   = i;            % store mean to all indices on the current slice of DIM
    [ o{1:nargout(F)} ] = F( X(searchix{:}) ); % this is tha magic that takes all output args of F for this category
    for j=1:length(o) % for each output arg of F
      mu{j}(saveix{:}) = o{j};  % compile into a given slice of output
    end
  end
elseif islogical(categ)
  for i=1:size(categ,2) % for each column of categ
    searchix      = all_indices;
    searchix{DIM} = categ(:,i);
    saveix        = all_indices;
    saveix{DIM}   = i;
    [ o{1:nargout(F)} ] = F( X(searchix{:}) );
    for j=1:length(o)
      mu{j}(saveix{:}) = o{j};
    end
  end
end
varargout=mu;
