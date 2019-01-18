function [p t stats terms anovaY anovaX] = anovanTable(T, varargin)
% [p t stats terms anovaY anovaX] = anovanTable(T, varargin)
%
% this is the same as ANOVAN(Y,X,varargin), and performs an N-way ANOVA.
% where there is one observation per condition, simply pass the array of
% data, with one dimension per factor.
% 
% The input 'T' is an N-dimensional array, where each
% dimension is a factor in the ANOVA. Dimension 1 becomes factor 1,
% dimension 2 becomes factor 2 and so on...
% Elements with 'nan' are ignored, so you can compare groups of different 
% sizes by padding smaller group's data with nans.
% 
% If there are serveral observations per condition, use an extra array
% dimension to store repeated measures, and specify to collapse across this
% dimension with the parameter 
%   'collapse', DIM
% where DIM is a list of dimensions that are collapsed across.
%
% all ANOVAN named-parameters are passed on, e.g. 'varnames', 'random',
% 'nested', 'model' etc.  But of course, 'continuous' predictors will only 
% work as linearly spaced values, as the levels are specified as 1,2,3...
% etc.
% 
% the first 4 return values are the same as anovan.
% the anovaY and anovaX results are the actual values passed to anovan,
% i.e. the vector-flattened data and the corresponding regressor set
%
% Extra Groupings: This allows each datapoint to be grouped in other
% different ways, so that you can specify extra factors in the ANOVA. 
% That is, you can regroup a subset of the factors already specified. Choose a
% couple of dimensions to regroup in a different way, and specify a 
% new matrix of levels for the new factor. 
%
%  extraGroupingsDimension:  a column vector of dimensions
%       i.e. all( extraGroupingsDimension < ndims(T) )
%       The dimensions define a "slice" of the main anovan table that is to
%       be regrouped in a different way. For example if dimensions 1 and 2 
%       represent subject and drug-state, and the drug-state was given on
%       different sessions for different subjects, then a new factor can be
%       added:
%           T ( SUBJ, DRUGSTATE, FACTOR1, FACTOR 2 )
%           extraGroupingsDimension = [ 1; 2 ]
%           extraGroupings ( SUBJECT, DRUGSTATE ) = 1 or 2 indicating which
%           session the corresponding data  T( i,j,:,: )   came from.
%
%  extraGroupings: The extra groupings matrix should then be the same 
%       size as the "SLICE" of T defined by the specified dimensions.  
%       i.e.
%         size( extraGroupings ) == [ size(T,DIM(1)), size(T,DIM(2)), ... ]
%       The levels specified in the matrix are used in the final regressor X
%       in an extra column, and therefore you should specify the  name 
%       of the new regressors as extra elements at the end of 'varnames'.
%
%  ignorenan: This means that all nans in the dependent variable ('y') 
%       will be removed before constructing the data for ANOVAN.
%       This usually has no effect except for the return values in anovaY
%       and anovaX, which will have fewer elements.

IGNORENAN=false;
% check for collapse parameter
remove=[];

for(i=1:length(varargin))
  if(strcmpi(varargin{i}, 'collapse'))
    collapse = varargin{i+1};
    remove=[remove i i+1];
  end
  if(strcmpi(varargin{i}, 'extraGroupingsDimension'))
    extraGroupingsDimension = varargin{i+1};
    remove=[remove i i+1];
  end
  if(strcmpi(varargin{i}, 'extraGroupings'))
    extraGroupings = varargin{i+1};
    remove=[remove i i+1];
  end
  if(strcmpi(varargin{i},'continuous'))
    continuous = varargin{i+1}; % don't discard - keep for anovan
  end
  if(strcmpi(varargin{i},'varnames'))
    varnames = varargin{i+1}; % don't discard.
  end
  if(strcmpi(varargin{i},'ignorenan'))
    IGNORENAN=varargin{i+1};
    remove = [remove i i+1];
  end
end
varargin(remove)=[];
if(~exist('collapse','var')) collapse=0;end
extraGrp = exist('extraGroupings','var');
if(xor(extraGrp, exist('extraGroupingsDimension','var')))
  error('extra groupings not specified - please spcify both dimension and levels');
end
% if there are multiple columns in extraGroupingDimension, this counts as multiple extra groupings!
if extraGrp, extraGrp = size(extraGroupingsDimension,2); end 
if(extraGrp>1) error('multiple regroupings not yet supported'); end

s=size(T);
N=length(s);
% collapse across dimesions => there are fewer anova factors
if(collapse) N=N-length(collapse); end; 
if(extraGrp) N=N+extraGrp; end
if exist('varnames','var') && N<length(varnames) % this happens if the final anovan dimension is a singleton
  warning('Table data contains fewer factors than expected.  Final singleton dimension(s) assumed.');
  N=length(varnames);
end
anovaY=[]; anovaX=[];
for(i=1:prod(size((T)))) % for every datapoint
  % ignore nans?
  if IGNORENAN
    if isnan(T(i)), continue; end;
  end
  anovaY=[anovaY; T(i)]; % create Y data
  % get a level for each factor
  [j(1) j(2) j(3) j(4) j(5) j(6) j(7)]=ind2sub(s, i); 
  % remove factors that are collapsed
  if(collapse) j(collapse)=[]; end
  if(extraGrp) % for k=1:extraGrp
    % dimIndex is the factor levels for the extra-grouping-dimension factors for the current datapoint
    dimIndex = num2cell(j(extraGroupingsDimension)); 
    % these levels are used to index into the extra-groupings factor levels
    j(N+1-extraGrp:N) = extraGroupings(dimIndex{:});
  end
  anovaX=[anovaX; j(1:N)]; % create X data
end
if exist('continuous','var') % zscore each column of the continuous variables
  for i=1:length(continuous)
    anovaX(:,continuous(i))=zscore(anovaX(:,continuous(i)));
  end
end
% do the anova
[p t stats terms]=anovan(anovaY, anovaX, varargin{:});

