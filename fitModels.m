function [params, ll, aic,  bic, out, extraInfo] = fitModels( ...
  modelFunctions, modelParamRanges, datasets,   varargin )
% [params, aic, ll, out] = fitModels( modelFunctions, modelParamRanges, datasets, ... )
%
% fit a list of models to a list of datasets, given the likelihood functions.
%
% modelFunctions contains the likelihood functions.
%        It takes a parameter vector and a dataset as inputs, 
%        and it should return a likelihood (0 to 1) as its first output. 
% 
% The fitting process calls    
%        p = modelFunctions{i} ( modelParams{i}, datasets{j} )
% 
% Your model function should return a vector of probabilities of the data 
% given the model. These probabilities are then multiplied together, and 
% the most probable model parameter s are found for each model.
%
% modelParams is drawn from the range specified in modelParamRanges{i}
% which for each model, should have two rows, for the minimum and maximum of
% each parameter. Note that these are not constraints, merely the starting
% points for the optimisation. Optimisation calls fminsearchs by default.
% 
%
% e.g. a logistic fit of binary data:
%
%   % create random dataset
%   predictor  = rand(1000,1);                          % column of independent variable: 0 to 1  
%   response   = ( (predictor + randn(1000,1))  ...     % corresponding binary dependent variable: 
%                  > 0.5 ) * 2 -set 1;                  %  either -1 or +1.
%   data       = [predictor, response];                 % each row is one measurement  
%   
%   % logistic model of data: P(response) = 1/(1+e^-bx)
%   % The likelihood function gives P when response == 1, and 1-P when response == -1.
%    likelihood = @(data,params) data(:,2) .* ...       % positive or negative depending on outcome 
%       ( 1 ./ ( 1+exp(-(data(:,1)-params(2))./params(1)) ) - 0.5   ) + 0.5 ;
%   % run the model fitting, with param1 between 0 and 10, and param2 between -10 and +10. 
%   [p, aic] = fitModels( {likelihood}, {[ 0 -10 
%                                         10  10 ]}, {data} );
%
% returns: 
%    PARAMS { MODEL } ( DATASET, PARAMETER )   Best fitting parameters
%    AIC    ( MODEL , DATASET )                Akaike information criterion
%    LL     ( MODEL , DATASET )                Log likelihood at optimum
%    OUT    { MODEL , DATASET } { NARGOUT }    All outputs of the model function at optimum 
%    BIC    ( MODEL , DATASET )                Bayesian information criterion
%    EXTRA.LEAVEOUT_PARAMS { MODEL } ( SUBJECT, PARAM, SUBSAMPLES )
%    EXTRA.LEAVEOUT_LL     ( MODEL, SUBJECT, SUBSAMPLES )
% 
% options:
%    optimset:   pass this value to fminsearch, as the optimisation settings
%    nargout:    number of arguments to collect from the model function, when
%                the optimum parameters are found (returned as 'out')
%    niterations: number of iterations of fminsearchs 
%    constrain:  if true, then ensure the paramaters lie within the
%                specified modelParamRanges.  Constrain can also be a
%                matrix of parameter ranges, just like modelParamRanges.
%    leaveOut:   Number of trials to leave out for subsampling. This allows
%                more robust estimates, by fitting the data with some rows
%                of each dataset missed out. 
%                This enables us to estimate the robustness of the fit, and
%                also the covariance of parameter estimates. The results
%                will be in extra.leaveOut.
%    leaveOutN:  Number of times to re-try fitting data with subsamples of
%                trials. 
% sanjay manohar 2013


%%%%%% Determine options
CALC_COV            = false;         % calculate covariance of parameters, using mlecov?
MINIMISATION_METHOD = 'fminsearchs'; % either 'fmincon', 'fminsearchs', 'ga' or 'mle'.
DEBUG_NANS          = false;         % break to debugger if any likelihood is 'nan'
SUBSAMP_MAXITER     = [];            % number of iterations in fminsearch when subsampling. 0 = same
VERBOSE             = 3;             % allow printing output during fitting.

[ optimset, nargouts ,niterations, constrain, leaveOut, leaveOutN] = parsepvpairs(...
  {'optimset','nargout','niterations', 'constrain', 'leaveOut', 'leaveOutN'},...
  {  [],        [],        10,              [],       [],          []},varargin{:});
if ~isempty(niterations), optimargs = {niterations}; else optimargs={[]}; end
if ~isempty(optimset),    optimargs = [optimargs ,{optimset}]; end
if ~isvector(datasets) && ismatrix(datasets) % datasets is a 2D matrix?
  szdata = size(datasets);                   % flatten datasets
  datasets=datasets(:);                      % and get size of each one
  datasize = permute( reshape( cellfun(@length, datasets), szdata ) , [5 1 2 3 4] );
else szdata=[];                              % datasets is a vector
    datasize = cellfun(@length,datasets); 
end
% for each dataset: work out subsampling
for j=1:length(datasets) % how many 'leave-out' combinations, for each dataset?
  nTrials(j) = size(datasets{j},1);  % number of trials in this dataset
  % if not specified, use the number of combinations, if <500.
  if leaveOut 
    if isempty(leaveOutN), leaveOutNd(j) = min(nchoosek(nTrials, leaveOut), 500); 
    else                   leaveOutNd(j) = leaveOutN;
                           if leaveOutN==0, leaveOut=false; end
    end
    leaveOutN = mean(leaveOutNd);
  end
end
if ( islogical( constrain ) && constrain ) || ( isnumeric(constrain) && isscalar(constrain) && constrain==1 ),
  constrain = modelParamRanges;
end

% for each model: determine number of parameters and outputs
Nm = length(modelFunctions);         % number of models
for i=1:Nm; 
  if isempty(nargouts)               % if the number of outputs is not specified,
    No(i) = nargout(modelFunctions{i}); % get number of function outputs using nargout
    if No(i)<0, No=abs(No);       % negative means variable number of outputs...
      warning('Variable outputs for model %g. Only some will be captured',i); end
  elseif numel(nargouts)==1,    No(i) = nargouts; % else use provided number of args
  else No(i) = nargouts(i);          % you can specify this for each model, if desired
  end
  Np(i,1) = size(modelParamRanges{i},2); % number of parameters
end

%%%%%% Fitting
%%%%%% To use multiple cores, change the i loop to 'parfor'.
parfor i = 1:Nm                         % for each model ( can be parfor )
  % how many outputs should we collect from the model functions?
  % Now go through each data set and fit it
  leaveOut_params = [];   leaveOut_ll = [];   params_i = []; % collate results from parallel threads
  ll_i = []; timesTaken_i = []; out_i = {}; covar_i = {}; 
  leaveOut_best_ll = []; leaveOut_best_params = []; 
  modelFunction_i = modelFunctions{i};  % store separately to prevent 'communication overhead'
  for j=1:length(datasets)           % for each data set ( *** CAN BE PARFOR ***)
    starttime=tic;
    try 
      if ~isempty(datasets{j})         % is there data in this set?
        [params_ij,   ll_i(j), out_i{j}, covar_i{j}] = runFit( ... % Run fitting...
          modelFunction_i, datasets{j}, modelParamRanges{i}, constrain, ...
          optimargs, No(i) , MINIMISATION_METHOD, CALC_COV );
        if ~isempty(params_ij),        % is it a 'no-parameter model'? (i.e. Np(i)==0)
          params_i(j,:) = params_ij;   % if not, store the fitted parameters
        end
        % if there's a nan-value in model ouptuts, switch to debugger?
        if DEBUG_NANS && isnan(ll_i(j)), keyboard; end
        timetaken       = toc(starttime);        % check how long we took,
        timesTaken_i(j) = timetaken;             % and print message if needed
        if (timetaken>1 && VERBOSE>0) || VERBOSE>2,
          fprintf('model %g/%g, dataset %g/%g (%0.1g): LL=%g\n',i,Nm,j,length(datasets), ...
            100 * ((i-1)*length(datasets)+j) / (Nm*length(datasets)), ll_i(j) );
          if leaveOut, fprintf('beginning sampling model %g, dataset %g, estimated time %g s\n', ...
              i, j, timetaken*leaveOutN); end
        end
      else                             % No data in set.
        out_i{j}=cell(1,No(i)); ll_i(j)=nan;
        covar_i{j} = []; out_i{j}=[]; timesTaken_i(j)=nan;
        params_i(j,:)          = nan( 1,Np(i) );
      end
      % is there enough data to run a leave-some-out analysis?
      % if not, skip this bit
      lo = leaveOut; 
      if leaveOutN+5 > nTrials(j), lo = false; end
      if lo                            % do leave-out subsampling?
        subargs = optimargs;           % optimisation options for the subsampling
        if SUBSAMP_MAXITER>0,          % make the subsamples have fewer iterations?
          subargs{2}.MaxIter  = SUBSAMP_MAXITER;
        end
        leaveOut_ll_j     = nan( leaveOutN,1 );
        leaveOut_params_j = nan( Np(i), leaveOutN );
        for k=1:leaveOutN              % for each subsample
          if ~isempty(datasets{j})     % any data?
            subdata = datasets{j};       % create a subsampled dataset,
            subdata( randsample(nTrials(j), leaveOut), : ) = []; % removing some random trials.
            % fit subsampled data. don't calculate covariance or store extra outputs.
            % use simple 'fminsearch' for maximum speed!
            [leaveOut_params_jk, leaveOut_ll_j(k) ] = runFit( ...
              modelFunctions{i}, subdata, modelParamRanges{i}, constrain, ...
              subargs, 0, 'fminsearchs', false );
            if ~isempty(leaveOut_params_jk)  % is it a 'no-parameter model'?
              %leaveOut_params(j,:,k) = leaveOut_params_jk;
              leaveOut_params_j(:,k) = leaveOut_params_jk;
            end
          else                         % no data in subset
            leaveOut_params_j(:,k) = nan( Np(i), 1 );
            leaveOut_ll_j(k)       = nan;
          end
        end % next subsample k
        leaveOut_params(j,:,:) = leaveOut_params_j;
        leaveOut_ll(  j,:)     = leaveOut_ll_j;
        [~,best_k]             = min(leaveOut_ll_j);
        if ~isempty(leaveOut_params_j)
          leaveOut_best_params(j,:)  = leaveOut_params_j(:,best_k);
        end
        leaveOut_best_ll(j)        = leaveOut_ll_j(k);
      else
        leaveOut_best_ll(j) = nan;
        leaveOut_params(j,:,:) = nan;
        leaveOut_ll(j,:) = nan;
        leaveOut_best_params(j,:) = nan;
      end % if leaveOut
    catch me
      me % skip over errors
    end
  end % next data set
  params{i}                    = params_i;          % store the values for this model 
  ll(i,:)                      = ll_i;              % for all datasets
  extraInfo_leaveOut_params{i} = leaveOut_params;   % including all the "leave-out" 
  extraInfo_leaveOut_ll(i,:,:) = leaveOut_ll;       % subsample fits
  timesTaken(i,:)              = timesTaken_i;
  covar(i,:)                   = covar_i;
  out(i,:)                     = out_i;
  extraInfo_leaveOut_best_ll(i,:)   = leaveOut_best_ll;     % best of the leave-out fits
  extraInfo_leaveOut_best_params{i} = leaveOut_best_params;
  % save intermediate results (for models fitted so far), in case of emergency.
  % don't save modelFunctions - it doesn't seem to work
  if 1 % this should be 0 if you use parallel 'parfor' for the i loop (models)
    try
      % save fitModelsTemp    modelParamRanges datasets   params ll out timesTaken covar extraInfo
    catch mexcp
      warning('Could not save temporary fit results file.');
    end
  end
end; j=[]; % next model
save fitModelsTemp    modelParamRanges datasets   params ll out timesTaken covar 
extraInfo.leaveOut_params      = extraInfo_leaveOut_params;
extraInfo.leaveOut_ll          = extraInfo_leaveOut_ll;
extraInfo.leaveOut_best_ll     = extraInfo_leaveOut_best_ll;
extraInfo.leaveOut_best_params = extraInfo_leaveOut_best_params;
save fitModelsTemp    modelParamRanges datasets   params ll out timesTaken covar extraInfo

%%%%%% look for outliers and correct them
if length(datasets)>10
  THRESH = 2;     % standard deviations for parameter to be considered OK
  for i=1:Nm                % each model
    bad=[];
    for j=1:Np(i)           % each parameter
      parm=params{i}(:,j);  % for each dataset, is it bad?
      bad(j,:)=(parm>nanmean(parm)+THRESH*std(parm)) | (parm<nanmean(parm)-THRESH*nanstd(parm));  
      % bad ( parameter, dataset)
    end
    % also check log likelihood.
    bad = any(bad);         % flag datasets with any bad parameter for this model
    bad = bad' | nanzscore(ll(i,:)') < -THRESH; % also see if Log likelihood is low
    bad = find(bad);        % list of bad datasets
    if any(bad) && timesTaken(i,1)>2 && VERBOSE>0, fprintf('improving %g model fits in model %g\n',length(bad), i); end
    for k=1:length(bad)     % for each bad datset
      warning('improving poor fit for dataset %g',bad(k)); 
      % redo the fitting just as for the main fit. Perhaps second time lucky? 
      [p1,ll1,o,co]=runFit( modelFunctions{i}, datasets{bad(k)}, modelParamRanges{i}, ...
        constrain, optimargs, No(i), MINIMISATION_METHOD, CALC_COV ); 
      if ll1 > ll(i,bad(k)) % is it an improvement on the old fit?
        ll(i,j)             = ll1; % if so, store the new parameters
        params{i}(bad(k),:) = p1;
        covar{i,j}          = co;
        out{i,j}            = o;
      end
    end % next bad dataset in the model
  end % next model
end % if

%%%%%%% Prepare results
% reshape results to original data shape (if datasets is a matrix, rather than vector)
if ~isempty(szdata) 
  for i=1:length(modelFunctions)
    params{i} = reshape( params{i}, [szdata  Np(i)] );
  end
  ll=reshape(ll, [Nm, szdata]);
  out=reshape(out,[Nm szdata]); 
  if leaveOut % extra info to resize from subsampling:
    for i=1:Nm   % for each  model
      % last two elements of size are: nParams, leaveOutN
      szleave = size(extraInfo.leaveOut_params{i}); szleave=szleave(end-1:end); 
      try
        extraInfo.leaveOut_params{i} = reshape(extraInfo.leaveOut_params{i}, [szdata szleave]);
      catch exc
        exc
      end
    end
    
    extraInfo.leaveOut_ll = reshape(extraInfo.leaveOut_ll, [Nm  szdata leaveOutN]);
  end
  timesTaken  = reshape(timesTaken,  [Nm  szdata] );
end

% calculate AIC for each model and each dataset
aic = bsxfun(@plus, -2*ll, 2*Np);
% calculate BIC
% how many datapoints are there in each dataset?
% this assumes that the ROWS of datasets are individual data points.
bic = bsxfun(@plus, -2*ll, bsxfun(@times, Np, log(datasize)));
extraInfo.timesTaken = timesTaken;
extraInfo.covar      = covar;


%%
function [params, ll, out, cov] = runFit( modelFunction, dataset, modelParamRanges, ...
    constrain, optimargs, No, MINIMISATION_METHOD, CALC_COV )
% Run the fitting once, on one dataset.
%  modelFunction: function giving likelihood of each trial, given params
%                 vector and a datasets matrix
%  modelParamRanges / constraints: initial range and constraints for params
%  No: number of outputs from model function
%  

if isempty(modelParamRanges)      % model with no parameters?
  p1=[]; nll1=-sum(log(modelFunction([],dataset))); cov=nan;
else                       % Fit parameters.
  % likelihood function: use only the first output of the function.
  % negative log likelihood of parameters
  if isempty(constrain) % no constraints
    nll = @(params) -sum(log(eps+    modelFunction( params, dataset )   ));
  else % constrain parameters to lie in range
    nll = @(params) -sum(log(eps+    modelFunction( params, dataset ) ...
      .*  prod(1*(   (params >= modelParamRanges(1,:)) ... if any condition not met, return zero.
      & (params <= modelParamRanges(2,:)) ) )        ));
  end
  % optimise
  switch MINIMISATION_METHOD
    case 'fminsearchs'     % use sanjay fminsearchs that tries several starting points
      [p1,nll1]=fminsearchs( nll, modelParamRanges, optimargs{:} );
    case 'mle'             % use matlab built int mle fitting.
      [p1, p1ci] = mle( zeros(size(dataset,1),1) , 'nloglf', @(p,q,r,s) nll(p), 'start', mean(modelParamRanges) );
      nll1 = nll(p1);      % get the value of the likelihood function at the optimum
    case 'fmincon',        % use matlab constrained fminsearch with upper/lower bounds
      % choose random initial parameters within range
      p0 = rand(1,size(modelParamRanges,2)).*diff(modelParamRanges) + modelParamRanges(1,:);
      [p1, nll1] = fmincon( nll, p0,  [], [], [], [], modelParamRanges(1,:), modelParamRanges(2,:), optimargs{2} );
    case 'fminsearch',     % basic fminsearch
      p0 = rand(1,size(modelParamRanges,2)).*diff(modelParamRanges) + modelParamRanges(1,:);
      [p1, nll1] = fminsearch( nll, p0, optimargs{2} );
    case 'ga'              % use genetic algorithm (global optimization toolbox)
      [p1, nll1] = ga( nll, size(modelParamRanges,2),  [], [], [], [], modelParamRanges(1,:), modelParamRanges(2,:), [],[], optimargs{2} );
  end
  if CALC_COV && length(p1)>1 % covariance of estimates, if multiple parameters?
    cov = mlecov(p1, zeros(size(dataset,1),1), 'nloglf',@(p,q,r,s) nll(p) );
  else cov=nan;
  end
end
% store results
params = p1;     % best parameters
ll     = -nll1;  % log likelihood of model
if No % If other outputs requested, re-run the model to get all its outputs, at the optimum
  [out{1:No}] = modelFunction( p1, dataset );
end


