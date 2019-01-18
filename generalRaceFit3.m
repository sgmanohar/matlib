function [parameters logLikelihood solution] = generalRaceFit2( ...
  RTs, conditionVector,...
  baseRaceParameters, varyParametersFunction, initialParameterVector, ...
  varargin)
% function [parameters, logLikelihood, solution] = generalRaceFit( 
%      RTs, conditionVector,...
%      baseRaceParameters, varyParametersFunction, 
%      initialParameterVector, [options] )
% 
% Fit a general race model to the reaction time data, using 
% maximization of likelihood.
%
% RTs: a matrix of reaction times. Each row is a trial, and each column
%      corresponds to one particular motor output from the race model.
% 
% conditionVector: One value for each data trial. 
%      The first (or only) column must be an integer, 
%      indicating the condition for that trial; it signifies which element
%      of baseRaceParameters to use.
%      The condition numbers must start at 1.
%      If there is only one condition in the given data, you can use [].  
% 
% baseRaceParameters: a structure array where baceRaceParameters(condition)
%      contains the base parameters for the race model for that condition.
%      This allows different (unfitted) race parameters for each condition.
%      If there is only one condition in the data, then an simple (scalar) 
%      structure will do. 
%      There must be as many structures as there are conditions (as given
%      in conditionVector).
%
% varyParametersFunction: a @function that creates a race model parameters
%      structure based on a base structure, using a parameter vector.
%      This determines what values are varied in the race model.
%      Parameters can be changed in a condition-dependent manner if needed,
%      or can be condition-independent.
%      varyParametersFunction = @(baseParameterStructure, parameterVector)
%
% initialParameterVector: a vector of the initial values of the parameters.
%      This should be of the length accepted by your
%      varyParametersFunction.
% 
% options:
%      'fittingMode': optional parameter to specify fitting mode, can be:
%                     'Genetic','Anneal','Fminsearch', or 'None'. default
%      'fiddle': set this to 1 to allow manual fiddling with the parameters
%                     after the fitting has been done. In 'fiddle' mode,
%                     type 'p1 = [ 1 1 1 1 etc ]; dbcont' to try out values
%                     of parameters in the fit.
%      'N_MODEL_TRIALS' (2000): number of trials in each run of the model. 
%                     Much larger sample size is needed for conditions with joint
%                     probability distributions, e.g. 10000.
%      'N_MODEL_REPEATS' (5): number of trials to repeat runs of the model when
%                     taking the average of log-likelihood. Reduces noise,
%                     useful when using a fitting method that is intolerant
%                     of noise.
%
% hint: if you have the parallel processing toolbox installed, ensure you 
% call 'matlabpool' beforehand to take advantage of multiple CPUs.
%
% RETURN:
%     parameters =    the solution of best fit
%     logLikelihood = the log likelihood of these parameters given the data
%                     i.e. sum(log(probability(data)))
%     solution =      extra information about the method used to solve the 
%                     fitting.
% sanjay manohar 2009

  VARIATION = 2; % the amount by which parameters vary stochastically
  

  %%%%%%%%%%%%%%
  % 
  % additional parameters for fittingMode = 'Genetic' :
  %    POPULATION = 10000; GENERATIONS = 60;
  %
  %%%%%%%%%%%%%%

% check parameters:
if(length(conditionVector)<2) conditionVector = ones(size(RTs,1),1); end % allow empty condition vector
if(size(conditionVector,1)~=size(RTs,1)) error('conditionVector must have same number of rows as RTs');end;
ucond = unique(conditionVector(:,1));
if(ucond(1)~=1) error('conditions must start with 1'); end
if(any(diff(ucond)~=1)) error('conditions must be consecutive integers starting at 1'); end;
if(length(baseRaceParameters)~=length(ucond)) error('baseRaceParameters must have one struct per condition'); end;
%if(~isvector(conditionVector)) error('conditionVector must be a vector');end;
[fittingMode fiddle MAX_ITER GENERATIONS N_MODEL_TRIALS, N_MODEL_AVERAGES, FITTINGARGS, VARIATION] = parsepvpairs( ...
    {'fittingMode', 'fiddle', 'MAX_ITER', 'GENERATIONS', 'N_MODEL_TRIALS', 'N_MODEL_AVERAGES','fittingArgs', 'VARIATION', }, ...
    {'GWMCMC',      0 ,       500,          30,            2000,             5,                 {} ,  0.1  }, varargin{:});
if(strcmpi(fittingMode,'Genetic')) GENETIC=1; % select fitting mode
elseif(strcmpi(fittingMode,'Anneal')) ANNEAL=1;
elseif(strcmpi(fittingMode,'Fminsearch')) FMINSEARCH=1;
end

% longest possible reaction time permitted. times longer than this are 
% counted the same as having no response.
L=10; % inf

solution=struct();
  
  global iter
  
  iter = 1; LL=-inf;
  % mean rate of rise for process 1 and 2
  p0 = initialParameterVector;  % starting parameters
  
  % function to display a set of parameters - use a CLOSURE in a GLOBAL!
  global displayFunction
  displayFunction=@(par) plotParamsFunction(RTs, conditionVector,...
                        baseRaceParameters, varyParametersFunction,  par);

                      
  %%%%%%%%%%%%%%
  % run minimum search
  %
  switch fittingMode
    case 'FMINSEARCH' % FMINSEARCH (original version)
    [p1] = fminunc( @(p) -loglikelihood(p,RTs, ...
      conditionVector, baseRaceParameters, varyParametersFunction, ...
      N_MODEL_TRIALS, N_MODEL_AVERAGES), ...
      p0, optimset('MaxIter',300, 'TolX', 0.05, 'TolFun', 100) );
    case 'GENETIC'  % GENETIC ALGORITHM
    [p1,LL,flags,result,population, scores] = ga( ...
                 @(p) -loglikelihood(p, RTs, conditionVector, ...
                                     baseRaceParameters, varyParametersFunction, ...
                                     N_MODEL_TRIALS, N_MODEL_AVERAGES),...
                 length(p0),...
                 gaoptimset('PopulationSize',MAX_ITER, ...
                     'EliteCount', floor(MAX_ITER/10), ... % 10% of population are elite
                     'InitialPopulation',[p0], ... % parameters initially range -2 to +2
                     'UseParallel','always',...
                     'StallGenLimit',3, ...
                     'Generations',GENERATIONS,...
                     'Display','iter' ...
                     ,'PlotFcns',{@displayFunctionGA}) ...
    );
    solution.population = population;
    solution.scores = scores;
    case 'ANNEAL'     % SIMULATED ANNEALING
    % this stochastically generates a new set of parameters
    %Generator = @(x)(x+(randperm(length(x))==length(x))*randn*VARIATION);
    Generator = @(x,T)( x + randn(size(x)) * VARIATION * 1 );
  
    % run the simulated annealing
    [p1, LL] = anneal_sanjay ...
      (@(p) -loglikelihood(...
          p,... the parameters for this attempt
          RTs,... the data
          conditionVector, ... the condition for each trial
          baseRaceParameters, ... the basic model parameters for the race for each condition
          varyParametersFunction,... a function to generate race parameters
          N_MODEL_TRIALS, N_MODEL_AVERAGES), ... 
      p0, ... initial values of parameters
      struct(   'InitTemp',VARIATION , ... annealing parameters
                'StopTemp',0.01 , ... (1e-8)
                'Generator', Generator, ...
                'Verbosity',2, ...
                'MaxTries', MAX_ITER, ...  maximum tries per temperature (300)
                'MaxSuccess', 20, ... maximum successful tries within one temp (20)
                'CoolSched', @(T)( 0.8*T ) ... cooling rate (0.8)
      ));
    case 'GWMCMC'
      NW = length(p0)*2;
      models0 = repmat(p0', 1, NW) + VARIATION'.*randn(length(p0),NW);
      likelihoodFun = @(p) loglikelihood(...
          p,... the parameters for this attempt
          RTs,... the data
          conditionVector, ... the condition for each trial
          baseRaceParameters, ... the basic model parameters for the race for each condition
          varyParametersFunction,... a function to generate race parameters
          N_MODEL_TRIALS, N_MODEL_AVERAGES);
      i=find(strcmpi(FITTINGARGS,'logprior'));
      if i>0, logPrior = FITTINGARGS{i+1}; FITTINGARGS(i:i+1)=[]; 
      else    logPrior = @(x)0; end
      [solution.models solution.logp, solution.result] = gwmcmc(...
        models0, { logPrior, likelihoodFun }, MAX_ITER, FITTINGARGS{:}...
      );
      p1 = solution.result.optimal;
      LL = solution.result.optimalLP;
      solution.llFun = likelihoodFun; 
    case 'NONE'% NO FITTING - TEST OUT INITIAL PARAMETERS
    p1=p0;  % No Fitting: use the initial parameter values, and just plot 
            % the model compared to the data. Useful for testing the
            % results of a fit.
    LL = loglikelihood1(p1, RTs, conditionVector, baseRaceParameters, varyParametersFunction, N_MODEL_TRIALS);
  end
  
  
  %%%%%%%%%%%%%%%% plot the solution
  fidd = 0;
  while(fidd) % loop to allow fiddling with parameters

    plotParamsFunction(RTs, conditionVector,...
                        baseRaceParameters, varyParametersFunction,  p1);
    if(fiddle) fprintf('[type "fiddle=0;dbcont" to exit.'); % break for fiddling
      keyboard
    else fidd=0; end
  end
  
  parameters    = p1; % assign output
  logLikelihood = LL;

  function LL = loglikelihood(params, RTs, conditions, modelParams, paramFn, N_MODEL_TRIALS, N_MODEL_AVERAGES)
    LL=nan(N_MODEL_AVERAGES,1); % preallocate
    for i=1:N_MODEL_AVERAGES 
      LL(i)=loglikelihood1(params,RTs,conditions,modelParams,paramFn,N_MODEL_TRIALS);
    end
    LL(LL==-inf)=-1e-100;
    LL=nanmean(LL);

    % VERBOSE???
    % fprintf('%g:\t[',N);fprintf('%g ',params);fprintf(']\t(LL = %g) \r', LL);
  
    
% calculate log likelihood of the paramteter vector 'params' given the data,
% and some baseline modelParams. 
% paramFn should generate a race-parameter structure from the vector of
% params.
function LL = loglikelihood1(params, RTs, conditions, modelParams, paramFn, N_MODEL_TRIALS)
  global iter 
  iter = iter + 1;     % keep track of iterations
  L=5;                % Long cutoff - the longest latency allowed
  DISPLAY=0;           % 1 = text, 2 = diagrams
  
  % run the model with the given parameters
  %%%%%% NUMBER OF DATA POINTS PER ESTIMATE: determines how long it takes!
  N = N_MODEL_TRIALS;  %*(1+floor(iter/10))
  %%%%%%
  ucond=unique(conditions(:,1));
  ncond=length(ucond);
  for(c=1:length(ucond))
    % create a parameter set specific to this condition and the current
    % value of params, using the user-provided paramFn.
    if size(conditions,2)>1  % just use the one condition
      p = paramFn(modelParams(c), params, ucond(c));
    else
      p = paramFn(modelParams(c), params);
    end
    % then run the race model 
    [x t] = generalRaceModel4(N,p);
    if all(isnan(t) ) % bad model parameters!
      LL=-inf;
      return;
    end
    % extract the data for this condition, for comparison
    data=RTs( conditions(:,1)==c ,:);
    
    
    ok_model = all(t>0 & t<L,2);
    ok_data  = all(data>0 & data<L,2);
    clear p1_j;
    if sum(ok_model)<5
      Lc(c) = -inf; continue;
    end
    for j=1:size(t,2); % marginal of just response 1 
      p1_j(:,j) = ksdensity( t( ok_model ,j), data( ok_data,j), 'function', 'pdf');
      p1_j(:,j) = p1_j(:,j) * mean(ok_model); % normalise by num valid samples in model
    end
    p1=p1_j(:,1); % use just response 1 marginal
    
    if(DISPLAY>1)
      subplot(ncond,2,2*(c-1)+1)
      reciprobit(nancat(2,t, data),[] );
      xlim([-10,0]); ylim([-5 5]);
      subplot(ncond,2,2*c)
      if(sum(f)>0)
        hist3(t(f,:)/sum(f));
      end
      drawnow
    end
    
    % log likelihood of the params given the dataset, for each condition
    Lc(c) = nansum(log( p1 + eps )) ...
          + nansum(~ok_data)*log(mean(~ok_model) + eps);
    % for the bad trials, use the probability of a bad model sample.
    
    % this may seem odd : I am adding a log probability to a log
    % density. The weighting looks arbitrary, as the latter has units
    % "probability per unit RT" and the former has units probability.
    % However it is absolutely fine, because this is effectively
    % multiplying the density by a probability, just as you would for a
    % censored distribution. 
  end
  % likelihood for all conditions
  LL = sum(Lc);
  if(DISPLAY)
    fprintf('%g:\t [%g %g %g] \t (LL = %g) \r', N , params , L);
  end


  
    
  function plotParamsFunction(  RTs, conditionVector,...
                        baseRaceParameters, varyParametersFunction,  p1)
    AUTO_SUBPLOT=0;
    N_DEMO_TRIALS=500;
    L=3;
    ucond = unique(conditionVector(:,1));
    ncond = length(ucond);
    for(i=1:ncond) % for each condition
      if size(conditionVector,2)>1
        p(i)=varyParametersFunction( baseRaceParameters(i),p1, ucond(i) );
      else
        p(i)=varyParametersFunction(baseRaceParameters(i),p1);
      end
      [x t] = generalRaceModel4( N_DEMO_TRIALS , p(i) );
      if(~all(isnan(x))) 
        d=RTs(conditionVector(:,1)==i,:); d(isnan(d))=inf; % filter appropriate conditions
        % split into 3 columns : Target, Distrac, Correc, for model and data 
        mT = nancat(2,t( t(:,2)>L ,1), t( (t(:,2)<L) & (t(:,2)<t(:,1)), 2 ),  t( (t(:,2)<L) & (t(:,2)<t(:,1)), 1 ));
        dT = nancat(2,d( d(:,2)>L, 1), d( (d(:,2)<L) & (d(:,2)<d(:,1)), 2 ),  d( (d(:,2)<L) & (d(:,2)<d(:,1)), 1 ));
        mT(isnan(mT))=inf; dT(isnan(dT))=inf;
        if(AUTO_SUBPLOT)
          subplot(ncond,2,(i-1)*2+1); % use a new subplot
        end
        if(1) %% PLOT RECIPROBITS
          reciprobit(nancat(2, mT, dT ),...
            [],... criteria: just use the columns
            [mean(t(:,2)>L), mean(t(:,2)<L), mean(t(:,2)<L), ...
            mean(d(:,2)>L), mean(d(:,2)<L), mean(d(:,2)<L)  ] , ...
            'plotfit',0); % don't draw linear fit
          xlim([-15,0]);
          legend off
        end
        if(AUTO_SUBPLOT)
          subplot(ncond,2,(i-1)*2+2);
        end
        if(0) %% PLOT CONDITIONAL CORRECTION RT
          mTok=mT(:,2)<2 & mT(:,3)<2;  dTok=dT(:,2)<2 & dT(:,3)<2;
          %scatter(mT(mTok,2),mT(mTok,3)-mT(mTok,2),'y.'); hold on
          %scatter(dT(dTok,2),dT(dTok,3)-dT(dTok,2),'ro');
          if(1) % DIFFERENCE BETWEEN CORRECTION TIME AND ERROR TIME?
            plotBinsQuantiled(mT(mTok,2),mT(mTok,3)-mT(mTok,2),7,'y');hold on;
            plotBinsQuantiled(dT(dTok,2),dT(dTok,3)-dT(dTok,2),7,'r');
          else
            plotBinsQuantiled(mT(mTok,2),mT(mTok,3),5,'y');hold on;
            plotBinsQuantiled(dT(dTok,2),dT(dTok,3),5,'r');
          end
          hold off;
          title('conditional RTs');
        end
      end;
      title(sprintf('condition %g',i)); xlabel('error RT'); ylabel('correction time');
    end
    if(AUTO_SUBPLOT)
      makeSubplotScalesEqual(ncond,2,[1:2:(ncond*2)]);
      makeSubplotScalesEqual(ncond,2,[2:2:(ncond*2)]);
      for(i=1:2:(ncond*2))
        subplot(ncond,2,i); set(gca,'xticklabel',-1./get(gca,'xtick'));
      end
      title('fit');
    end
    drawnow;
    p1       % display result 
    LL=loglikelihood1(p1, RTs, conditionVector, baseRaceParameters, varyParametersFunction, N_DEMO_TRIALS);
    fprintf('LL=\t%g\tBIC=\t%g\n', LL, -2*LL + length(p1)*log(size(RTs,1)) );
    
    
    function [state]=displayFunctionGA(options,state,flag,interval)
      topix = find(state.Score==min(state.Score),1);
      p=state.Population(topix,:);
      global displayFunction
      displayFunction(p);
      
