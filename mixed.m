function result=mixed( Xa, varargin )
%
% function [beta, psi, stats, B, allmodels ]=mixed( Xa, varargin )
% result =  mixed( X, ... )
%           X has the GROUP in column 1, predictors in subsequent columns, 
%           and the dependent variable (y) as the last column. OR: 
% result =  mixed( y, X, ... )
%
% Fit mixed models, trying to find the best combination of random effects.
% Calls nlmefit with a simple linear model.
%
% Unlike anovanTable, this treats all variables as continuous rather than
% categorical. Use x2fx to convert if requried. 
%
% If X is a n-dimensional array, the first dimension is treated as the grouping
% variable, and the remaining dimensions are the main effects.
% X can also be in 'long form', with a row for each observation. In this case,
% the first column is the grouping variable, and the final column is the
% measurement.
%
%  'varnames':     name of each dimension.
%  'group'   :     the grouping dimension (integer from 1 to ndims(Xa), default=1
%  'interactions': a cell array of interactions, e.g. {[2 3]} includes the
%                  interaction between dimensions 2 and 3
%  'modelfun':     the model function, f(phi,X) - see nlmefit
%  'fecombinations': if true, try out all combinations of regressors as fixed
%                  effects. Return the one with best BIC.
%  'verbose':      0, 1 or 2 for output level.
%  'covpattern':   'rm': repeated measures (compound symmetric; 'cs') (default) 
%                  'eye': diagonal only, i.e. no covariance, i.e. independent 
%                  'full': estimate all covariances separately
%  'plot':         show AIC as bars, and which parms in each model.
%  'randinter':    default true. Ensure a random intercept is always included 
% result: stats.re_model is a boolean array representing which random effects
%         were included in a model, with one value for each fixed effect.
% sgm 2015

% previously:
% function [beta, psi, stats, B, allmodels ]=mixed( Xa, varargin )

BEST = 'aic'; % which 'stats' to minimise: bic, aic, mse or logl

%%%% get design matrix, grouping, and data
if nargin>1 && isnumeric(varargin{1})   % both X and Y provided
  Y=varargin{1}; varargin=varargin(2:end); 
  X=Xa;
  if ~(iscolumn(Y) || isrow(Y)) && (iscolumn(X) || isrow(X))
    % interchange X and Y
    tmp=X;X=Y;Y=tmp; 
  end
else % only X provided: treat it as a n-dimensional table (see anovanTable)
  if ndims(Xa)>2, X=dePivot(Xa); else X=Xa; end
  % get data column = last column
  Y=X(:,end);  X(:,end)=[];
end
nd = size(X,1);  % num data points
[varnames   group   interactions   modelfun     ALL_COMBINATIONS, ALL_RE_COMBINATIONS, verbose, ...
  COVPATTERN   PLOT    ENSURE_RANDOM_INTERCEPT, X2FX]=parsepvpairs(...
{'varnames','group', 'interactions', 'modelfun', 'fecombinations','recombinations',   'verbose', ...
 'covpattern','plot',  'randinter',             'x2fx'}, ...
{     []       1          []          @(b,X)X*b'    false           false               2   ...
  'rm'         false         true               true
  } , varargin{:} );
np = size(X,2);      % num parameters
if isempty(varnames) % name the parameters
  varnames = arrayfun(@(x)sprintf('v%g',x),1:np,'uni',0);
end

%%%% expand any categorical factors, if X2FX is true
newX={};
for f=1:size(X,2) % check for categorical factors
  newX{f} = X(:,f); 
  if f==group, continue; end; % ignore group column (that is supposed to be categorical!)
  uf = unique(X(:,f)); uf(isnan(uf))=[];
  if diff(uf)==floor(diff(uf)) && length(uf)>2
    if ~X2FX
      warning( 'categories of factor %g are categorical with %g levels. Did you mean to use x2fx?', f, length(uf) );
    else  % @TODO  - expand categories into separate columns using x2fx
      % note this is done before interactions are calculated.
      for i=1:length(uf)
        newX{f}(:,i) = X(:,f)==uf(i);
      end
    end
  end
end
if X2FX % expand the new factors
  X=[newX{:}]; 
end


% build interaction terms
for i=1:length(interactions) % for each interaction
  term=ones(nd,1); vars=interactions{i}; termname = ''; % start with constant,
  for j=1:length(vars) % for each of the components of the interaction term
    term = term .* X(:,vars(j)); % multiply by that column
    termname = [termname varnames{vars(j)} '*']; % and create a sensible name
  end
  X=[X term]; varnames=[varnames termname(1:end-1)]; % store the new interaction term
  np=np+1;
end
G=X(:,group); X(:,group)=[]; varnames(group)=[];  % get group variable

% column of ones?
const=find(var(X)==0);
if ~any(const), warning('adding intercept to model'); X=[ones(nd,1), X];varnames=['c',varnames]; end

% turn off fitting warnings
oldwarning = warning('query','stats:nlmefit:UnableToDecreaseSSEinPNLS');
warning('off','stats:nlmefit:UnableToDecreaseSSEinPNLS');

% which regressors to include?
if ALL_COMBINATIONS % try each possible combination 
  % each row is a model, and each column refers to which of the predictors is
  % included.
  fe_models  = fliplr( dec2bin(0:2^(size(X,2)-1)-1)=='1' ); 
  fe_models  = [ones(size(fe_models,1),1) fe_models]==1;
else % Just run a single model with all regressors
  fe_models  = true(1,size(X,2));
end
i=1;
% for each column of X, does it have more than one value for a subject G? 
% (test on subject G==1).
fe_levelsPerSubject = cellfun(@(x)length(unique(x)), num2cell( X( G==first(unique(G)) , : )  ,1)) ;
fe_totalLevels      = cellfun(@(x)length(unique(x)), num2cell( X,1 ) );
if sum(fe_totalLevels==1)>1, error('too many constant terms'); end
% allow random effects on variables with more than one level per subject, and
% also on the constant term
re_allowRandom = (fe_levelsPerSubject>1) | (fe_totalLevels == 1); 
% ensure that there is a random intercept
if ENSURE_RANDOM_INTERCEPT
  re_ensureRandom = fe_totalLevels==1;
end
for j=1:size(fe_models,1) % for each combination of fixed included effects,
  fe_model = fe_models(j,:); 
  % if there are many fixed effects, it is not possible to enumerate the
  % random effect combinations.
  if sum(fe_model)>8 
    warning('too many fixed effects for model %g',j);
    continue;
  end
  % go through each possible combination of random effects
  if ALL_RE_COMBINATIONS
    re_models = fliplr( dec2bin([1:2^sum(fe_model)-1]')=='1' );
  else % if not, just allow all the random effects
    re_models = ones(1,sum(fe_model))==1;  
  end
  for k=1:size(re_models,1)
    time_tmp=tic;
    re_model = re_models(k,:); % which random effects to include
    % check each of the random effects is allowed, and all the required ones are
    % included.
    if ~all(re_allowRandom(re_model)), continue; end
    if any(re_ensureRandom(~re_model)), continue; end
    switch COVPATTERN
      case 'full',        covpattern=ones(sum(re_model));
      case {'cs','rm'},   covpattern= eye(sum(re_model))+1;
      case {'eye','diag'},covpattern= eye(sum(re_model));
      otherwise           if isnumeric(COVPATTERN) covpattern=COVPATTERN;
        else error('invalid covariance pattern');
        end
    end
    % create final design matrix
    bad = any(isnan(X(:,fe_model)),2) | isnan(Y);
    Xmodel = X(~bad, fe_model); 
    
    % run model
    try
      bad = any(isnan(X(:,fe_model)),2) | isnan(Y); 
      [beta{i},psi{i},st,B{i}]=nlmefit( X(~bad,fe_model), Y(~bad), G(~bad), ...
        [],@(b,X)X*b', zeros(1,sum(fe_model)),...
        'REParams',re_model, 'CovPattern', covpattern ,...
        'Vectorization','singlegroup');
      st.pvalue   = tcdf(beta{i}'./st.sebeta,st.dfe); % store result
      % convert the re_model flags, which are relative to the fixed effects
      % being used for this model, to refer to ALL the possible fixed effects.
      fe_indices = find(fe_model);  % which FE are we using?
      re_model_params = zeros(1,length(fe_model)); re_model_params(fe_indices(re_model))=true; 
      st.re_model = re_model_params; 
      st.fe_model = fe_model;
      st.varnames = varnames;
    catch me % error occurred
      if exist('stats','var')
        disp(me);
        st=ensureStructsAssignable(struct(),stats); % blank result
      else
        rethrow(me)
      end
    end
    stats(i)    = st;
    i=i+1;
    time_tmp=toc(time_tmp);
    if time_tmp>2 || (verbose>0)
      fprintf('model %g/%g\n', j*length(re_models)+k, length(re_models)*length(fe_models));
      if verbose>1
        disp(st)
      else
        fprintf('F%s,',varnames(fe_model)); fprintf('R%s,',varnames(re_model));
        fprintf('LL=%g,\tAIC=%g,\tBIC=%g\n',stats.logl,stats.aic,stats.bic);
      end
    end
  end
end
if strcmp(oldwarning.state,'on') 
  warning('on','stats:nlmefit:UnableToDecreaseSSEinPNLS'); 
end

[~,bestmodel]=min([stats.(BEST)]);
allmodels.beta=beta; allmodels.stats=stats; allmodels.psi=psi; allmodels.B=B;
beta=beta{bestmodel}; stats=stats(bestmodel); psi=psi{bestmodel};  B=B{bestmodel}; 

result.beta  = beta;
result.psi   = psi;
result.stats = stats;
result.B     = B;
result.allmodels = allmodels;

% chi-squared model comparison for each pair of models
aLL = [allmodels.stats.logl];   aDF = [allmodels.stats.dfe];  % get LL and dfErr for each model
dDF = bsxfun(@minus, aDF,aDF'); dLL = bsxfun(@minus,aLL,aLL'); % table of differences
dLL(dDF>0)=-dLL(dDF>0); dDF = abs(dDF); % if adding params, LL should be higher
p  = chi2cdf(2*dLL,dDF);    % use Wilks theorem: -2 * log( L1 / L2 ) ~ chi^2
p(p>0.5)=p(p>0.5)-1; p=-p;  % flip the p values so that params making the model better give < 0.05
result.model_comparinson_p = p;

if PLOT
  %%
  clf
  bar([allmodels.stats.(BEST)]); barax=gca; title(BEST)
  dots = vertcat(allmodels.stats.fe_model)';
  dots = [dots; nancat(1,allmodels.stats.re_model)'];
  pvals = nancat(1,allmodels.stats.pvalue)'
  dots = [dots; pvals<0.05|pvals>0.95]; % two one-tailed t-tests!
  imgax = axes('position',get(gca,'position').*[1,1,1,0.2]);
  dotlabels = varnames;
  himg = imagesc(dots);
  set(himg,'alphadata',~isnan(dots));
  set(imgax,'xlim',get(barax,'xlim'),'ytick',1:size(dots,1),'yticklabel',dotlabels);
  %axes(barax);
end
