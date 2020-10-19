function [ A, Lmodel, t ] = rmanova(T, varnames, varargin)
%    ANOVA = rmanova( T )
% [  ANOVA , model , table ] = rmanova( T, varnames )
% repeated measures anova using matlab's fitlme, using data in 
% n-dimensional array, or in table form.
%   T is either a N-dimensional array (factorial form data), 
%      or an N-column matrix (long-form data).
%     column/dimension 1 is subject 
%     column/dimension N is the data (dependent variable)
%     other colulmns are factors / predictors.
%  'categorical': which variables (columns or dimensions of T) are
%                 categorical. Subject (1) is categorical by default.
%                 Others are linear regressors.
% Extra parmeters are passed to fitgmle, e.g. 'covariancepattern'. 
% to do a pure ANOVA, use 'CovariancePattern','CompSymm'.
% 
% returns: model = glme object with parameters and AIC
%          A     = anova table for F-tests of categorical variables.
%          t     = data formatted into a table
%  -- requires fitglme (Matlab 2015 or later)
% sgm 2017

% which random effects to model?
RE = enum({'FACTORIAL','INTERCEPT'});
% intercept random effect is simplest model, with fewer parameters.
% 'factorial' adds random effects for each fixed effect.
re = RE.INTERCEPT;


NOCATEG = false; % prevent categoricalising factor variables? (unless explicitly specified)

if ndims(T)>2 || ... % if the data is multidimensional, or if it is a table 
    ...% and the first column is 
    ~(all(isinteger(T(:,1))) || ...% not all integers
      all(diff(diff(unique(T(:,1))))==0) ) % and not regularly spaced (like subjects)
  % then, convert the array into a table.
  T=dePivot(T);
  if exist('varnames','var') &&  length(varnames) == size(T,2)-1 % have they omitted the dependent var name?
    varnames{end+1}='y';
  end
end

T(:,[2:end-1]) = nanzscore(T(:,[2:end-1])); % zscore the factor columns
constcol = find(nanvar(T(:,1:end-1))==0); % any constant columns?
if any(constcol)
  warning('Intercept present by default, so constant terms removed');
  T(:,constcol)=[];
end


if ~exist('varnames','var') || ~iscell(varnames) || length(varnames)==0
  % default variable names
  varnames = ['subject', arrayfun(@(i)sprintf('factor%g',i),1:size(T,2)-2,'uni',0), 'y'];
end


subs=unique(T(:,1)); % get first column: subjects
if length(subs)<3  ||  ~all(floor(subs)==subs)
  % is it all integers, at least 2?
  error('first column/dimension should be subjects')
end


% convert to table
t = array2table(T,'variablenames',varnames);

vn = t.Properties.VariableNames; % get var names

i=strcmpi(varargin, 'categorical'); % check for 'categorical' parameter
i=find(i); % any match?
if ~isempty(i) % did they specify categorical variables?
  CATEG = varargin{i+1};
  varargin(i:i+1)=[];
  for i=1:length(CATEG) % make them categorical
    t.(vn{CATEG(i)}) = categorical( t.(vn{CATEG(i)}) );
  end
else % guess which are categorical
  for j=1:length(vn)-1 % don't change the last column to categorical! it's the dependent.
    uj = unique(t.(vn{j}));
    % if there's between 2 and 8 unique values,
    if length(uj)>1 && length(uj)<=8 && (j==1 || ~NOCATEG) % is it categorical?
      t.(vn{j}) = categorical(t.(vn{j}));
      if j~=1
        warning('treating %s as categorical',vn{j});
      end
    end
  end
end
factorial = varnames(2:end-1); % build the factorial design with vars 2:N-1
factorial = flat( [factorial; [repmat({'*'},1,length(factorial)-1), ' ']] );
switch re
  case RE.FACTORIAL % factorial random effects
    reparams = [factorial{:}]; % all random effects included
  case RE.INTERCEPT
    reparams  = '1'; % random intercept only
end
% build the model string
model = [ varnames{end}  ' ~ ' factorial{:} ' + ( ' reparams ' | ' varnames{1} ')' ];
if length(unique(T(:,end))) == 2 % dependent variable is binomial?  do logistic regression instead
  warning('using logistic regression');
  varargin = [varargin 'link','logit','distribution','binomial']
end
%%%%%%%%%%%%% Actually fit the model here %%%%%%%%%%%%%%%  
if exist('fitglme','file'),  fitfun = @fitglme;
else                         fitfun = @fitlme;   end
Lmodel = fitfun(t,model, 'dummyvarcoding','effects', varargin{:});
%model = fitglme(t,model, varargin{:});
try % if we have newer stats toolbox: use satterthwaite correction for DF
  A = anova(Lmodel,'dfmethod','satterthwaite');
catch % otherwise use residual DF 
  A = anova(Lmodel);
end

