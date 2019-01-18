function [corrcoef, corrp, varnames Y] = correlateStructFields(X, varargin)
% function correlateStructFields(X)
% attempts to correlate the fields in a struct array X.
%
% [ X(1).field1 , X(2).field1, X(i).field1 ... ]  are to be correlated against
% [ X(1).field2 , X(2).field2, X(i).field2 ... ]  
%
% We iterate through the fields of X, and build a matrix of values, where rows
% are elements of X, and columns are fields of X. 
%
% If 'field' is scalar, then a single value from each element of X is placed
% into a single column of the data matrix.
% if 'field' contains a vector or a matrix, then compute the main effects along
% the directions and use them in the correlation. For example, if X.field is a
% set of values [ a,b,c ], then the mean of abc and the linear slope (c-a) are
% calculted, and entered as two parameters into the regression.
% Resulting variables are named 'm' for the mean value, 'e' for the main effects, 
% and 'i' for interaction term.
% 
% If X is a single (scalar) structure, then the rows of each field are
% treated as individual records. Columns will be named by suffixing their number.
%
% passes a matrix of values to 'corr', and additional arguments are passed to
% corr. Results are
%  corrcoef   correlation coefficient from 'corr' - matrix of NxN values, where
%             N is the number of items correlated
%  corrp      p-value for correlation from 'corr'
%  varnames   a cell array of names for each row/column in the c
%
% sgm 2010

INTERACTS = true; % include interaction terms when X.field is a matrix
DRAW = true;      % show the correlation matrix
LIST = true;      % list the significant correlations
ALPHA = 0.05;    % threshold p-value

if ~isvector(X) | ~isstruct(X), error('X must be a struct vector'); end
f=fieldnames(X);
Y=[]; Yn={}; % the correlation values, Subject x Correland
if ~isscalar(X)
  NS = length(X);
  for i=1:length(f)
    % progress indicator for slow ones
    if length(f)>500 & mod(i,30)==0; fprintf('field %g of %g...\n',i,length(f)); end
    s=size(sq(X(1).(f{i}) ));
    if ~isnumeric(X(1).(f{i})), warning('ignoring field %s - non numeric',f{i}); continue; end
    if prod(s)==1 % scalar
      y=cat(1,X.(f{i}));
      yn=f(i);
    elseif sum(s>1)==1 % vector
      t=sq(nancat(3,X.(f{i}))); % matrix of val x subject
      b = ols(t, [ones(size(t,1),1), zscore([1:size(t,1)]')], [1 0; 0 1]);
      y=b'; % mean value and effect size (slope along vector(
      yn={['m' f{i}], ['e' f{i}]};
    elseif sum(s>1)==2 % matrix
      t=sq(nancat(3,X.(f{i}))); % v1 x v2 x subject
      x1=repmat(zscore([1:size(t,1)]') , [1, size(t,2)]); % regressors for dim1
      x2=repmat(zscore([1:size(t,2)]')', [size(t,1) ,1]); % and dim2
      shape1 = size(t,1)*size(t,2); % num of data points per subject
      regr=[ones(shape1,1) reshape(x1, shape1, 1) reshape(x2, shape1, 1) ];
      yn  = {['m' f{i}], ['e1' f{i}], ['e2' f{i}]};
      if INTERACTS % add interaction term?
        regr = [regr regr(:,end).*regr(:,end-1) ];
        yn   = [yn {['i' f{i}]}];
      end
      t1 = reshape(t, shape1, NS); % combine dim1 and dim2 into dim1.
      bad = any(isnan([regr t1]),2);
      b = ols(t1(~bad,:), regr(~bad,:), eye(size(regr,2))); % calculate multiple regrssion
      y=b';
    else
      warning('ignoring field %s - too big',f{i});
    end
    Y  = [Y y];
    Yn = [Yn yn];
  end
else % X is a single structure (scalar): treat each field as containing rows for subject
  for i=1:length(f)
    vals=X.(f{i});
    rows=size(vals,1);
    cols=size(vals,2);
    if rows>1, NS=rows; 
      Y=[Y vals];
      if cols>1 % several columns? name them
        Yn=[Yn  arrayfun( @(x)sprintf('%s:%g',f{i},x),1:cols ) ] ;
      else % one column
        Yn=[Yn f{i}];
      end
    else
      warning('ignoring field %s - no rows',f{i});
    end
  end
end
if size(Y,2)>2000, fprintf('correlating %g variables...\n',size(Y,2)); end
[corrcoef, corrp] = corr(Y,varargin{:});
varnames = Yn;
if DRAW
  imagep( corrp.*sign(corrcoef), Yn );
end
if LIST
  c=corrp; c=c.*(~tril(ones(size(c))));
  [si,sj]=find(c<ALPHA & c>0);
  T={};
  for i=1:length(si)
    T=[T; {Yn{si(i)}, Yn{sj(i)}, corrcoef(si(i),sj(i)), corrp(si(i),sj(i)) }];
  end
  if ~isempty(T)
    displaytable(T, {'v1','v2','rho','p'},[20,20,12,12]);
  else 
    fprintf('no correlations\n');
  end
end
