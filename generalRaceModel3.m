
function [responseType  RTs distr seqtrig seqt ] = generalRaceModel3( ...
    nTrials, ...
    p, ...
    data, ...
    options )
% [responseType RT] = generalRaceModel(nTrials, parameters, data, options );
% Run Simulation
% Simulates an arbitrary number of rise-to-threshold units.
% They have independent distributions of rate-of-rise and s.d.
% The threshold is assumed to be the same for all units.
% 
%   nTrials:     = number of trials to simulate
%   parameters: 
%     N          = number of competing signals    
%     prior      = the initial value of each signal's accumulator
%                  Default = zeros(1,N)
%     startTime  = the moment at which each signal begins its rise to
%                  threshold. If it is +inf, the signal doesn't occur.
%                  Default = zeros(1,N)
%     meanRate   = the mean initial rate of rise of each signal
%                  Default = ones(1,N)
%     stdRate    = standard deviation of the rate of rise of each signal
%                  Default = ones(1,N). 
%                  To make the rates of rise of one process covary with 
%                  those of another process, specify the square covariance 
%                  matrix as stdRate.
%     threshold  = the threshold. Constant and Identical for all processes.
%     stdPrior   = the variability of the prior, currently implemented as 
%                  a gaussian with mean 'prior'
%                  Default = 0, no variation in prior.
%                  use rectangularPrior = 1 to use a rectangular
%                  distribution instead of gaussian.
%     
%     trigger    = N x N matrix
%        triggers(i,j) is used when process i completes. If process i
%        reaches threshold, then process j's rate is set to this value. 
%        * If the trigger value is 0, then process j is halted.
%        * If the trigger value is nan, then the rate of process j is
%          unaffected
%        * The diagonal values (i,i) should be 0 for processes that cause
%          triggering. This response-type, '0', indicates no response (see
%          censoring, below)
%        * If the diagonal value (i,i) is positive, then a response is
%          emitted when this process i reaches threshold. giving  a 
%          responseType of this value. The simulation continues if the
%          other values are nan.
%        Default = diag(1:N), i.e. no triggering and each process ends the
%        race with a different response.
%        
%     mode      = SET_RATE = 1   triggers set the rate of rise of the specified signal
%                 MUL_RATE = 2   triggers multiply the rate of rise
%                 ADD_RATE = 3   triggers add on to the rate of rise
%                 SET_ACCU = 4   triggers set the accumulator value
%                 ADD_ACCU = 5   triggers add on to the accumulator value
%
%   data:        = data to fit to the model
%                  NOT YET IMPLEMENTED
%
%  for no stop signal, use ssd = inf
% return: 
%   for each trial (i.e. rows in the output)
%     responseType = the chosen action. This is a value from the main 
%        diagonal of 'triggers' (or an integer from 1 to N,  if no triggers
%        specified). This represents action signal that reaches threshold.
%     RT = the reaction time (nan for successful stop)
% 
% Options:
%   options.trialLength = maxiumu duration of a trial. 
%       Responses that are later than this cutoff are truncated as 
%       non-responses. This is needed because response processes could 
%       have a zero or negative rate  rise, and never reach threshold.
%       Default = 2
%   options.censor      = whether to use right-censoring for no-response
%       trials. A response type of 0 is considered a no-response.
%       if censoring is on, return the 'trialLength' as the RT for all 
%       trials where there is no response.
%
% Version
% 2. allow more than one trigger in sequence
% 3. allows more than one response to be made.
% 4. allows triggers to induce a multiplicative or additive change to 
%    other signals

% find number of signals
N  = nTrials;
if isfield(p,'N')
    NS = p.N;
elseif isfield(p,'prior')
    NS = length(p.prior);
elseif isfield(p,'meanRate')
    NS = length(p.meanRate);
elseif isfield(p,'stdRate')
    NS = length(p.stdRate);
elseif isfield(p,'startTime')
    NS = length(p.startTime);
else
    error('couldnt determine number of competing signals');
end
% set defaults
if ~isfield(p,'prior')
    p.prior = zeros(1,NS);
elseif length(p.prior)==1
    p.prior = p.prior * ones(1,NS);
end
if ~isfield(p,'meanRate')
    p.meanRate = ones(1,NS);
elseif length(p.meanRate)==1
    p.meanRate=p.meanRate * ones(1,NS);
end
if ~isfield(p,'stdRate')
    p.stdRate = ones(1,NS);
elseif length(p.stdRate)==1;
    p.stdRate=p.stdRate * ones(1,NS);
elseif ~isvector(p.stdRate) % is covariance matrix?
    if(det(p.stdRate)<=0 || p.stdRate(1)<=0)
      % p.stdRate=diag(abs(diag(p.stdRate))); % one way of resolving it!
      responseType = nan;
      RTs = nan;
      distr = struct(); seqtrig=RTs; seqt=RTs;
      return;
    end
end
if ~isfield(p,'startTime')
    p.startTime = zeros(1,NS);
elseif length(p.startTime)==1;
    p.startTime = p.startTime * ones(1,NS);
end
if ~isfield(p,'threshold')
    p.threshold = 1;
end
if ~isfield(p,'trigger')
    p.trigger = diag(1:NS);
else
    % don't clear out nans in this case! see preamble above...
    % signalsThatTerminate = diag(p.trigger)>0; % zero all trigger effects for signals that terminate
    % clearRegion  = repmat(signalsThatTerminate,1,NS) & eye(NS)==0;
    % if any(any(p.trigger(clearRegion))) warning('trigger matrix has nonzero effects for terminating signals'); end
    % p.trigger(clearRegion) = 0; 
    if any(isnan(diag(p.trigger)))
        warning('response types of "trigger" should not be nan');
        p.trigger(isnan(p.trigger))=0;
    end
end
if ~isfield(p, 'stdPrior')
    p.stdPrior = zeros(1,NS);
elseif length(p.stdPrior)==1
    p.stdPrior = p.stdPrior * ones(1,NS);
end
if(~isfield(p, 'rectangularPrior'))
    p.rectangularPrior = 0;
end
  

SET_RATE=1;  % triggers set the rate of rise of the specified signal
MUL_RATE=2;  % triggers multiply the rate of rise
ADD_RATE=3;  % triggers add on to the rate of rise
SET_ACCU=4;  % triggers set the accumulator value
ADD_ACCU=5;  % triggers add on to the accumulator value

if ~isfield(p,'mode')
    p.mode = SET_RATE;
end
if ~exist('options','var')
    options=struct();
end
if ~isfield(options,'censor')
    options.censor=1;
end
if ~isfield(options,'trialLength')
    options.trialLength=2;
end
if  nargout>2 
    if ~isfield(options,'distributionBins')
        options.distributionBins=100;
    end
    if length(options.distributionBins) == 1
        options.distributionBins = linspace(0,options.trialLength,options.distributionBins);
    end
end


% run the model
% spmd % parallelise!
    
    % initial rate of rise for each process
    if(isvector(p.stdRate))
      R = repmat(p.stdRate,N,1) .* randn(N,NS) + repmat(p.meanRate,N,1);
    else % hopefully is covariance matrix.
      R = randn(N,NS) * chol(p.stdRate) + repmat(p.meanRate, N, 1);
    end
    T = zeros(N,1);             % cumulative time in each trial
    A = repmat(p.prior, N,1);   % cumulative activity of each process
    if(~p.rectangularPrior)
      A = A + randn(N,NS) .* repmat(p.stdPrior, N,1); % variability of prior
    else
      A = A + (rand(N,NS)-0.5)*2 .* repmat(p.stdPrior, N,1); 
    end
    
    Y = zeros(N,1);             % first response emitted
    RT = zeros(N,1);            % time of first response
    nY = zeros(N,1);            % total number of responses made
    nResponseTypes = max(diag(p.trigger)); % max index of responses
    RTs = nan*ones(N,nResponseTypes+1);       % the time that each signal reaches threshold
    % allow a column for zeros (convenient for storing; will be deleted
    % afterwards)
    
    % deal with start-time by compensating the initial value of accumulator
    % (a bit hacky!)
    A = A - repmat(p.startTime,N,1) ./ R;
    % sequential step times and response types
    seqt=[]; seqtrig=[];
    allEnded=0;
    timegapbefore=p.startTime;
    
    % calculate time for each signal to reach threshold
    while ~allEnded
        tsig = (p.threshold - A) ./ R + repmat(timegapbefore, N,1);   
        timegapbefore=zeros(1,NS); % change this if a delay is needed between stages of a process
        tsig(tsig<eps)=inf;        % processes with negative/zero rate never cross threshold
        tsig(R<eps)   =inf;
        % if all( Y>0 | all(tsig,2)==inf ) break;end; % exit if all either complete or nonterminating
        if all( all(tsig,2)==inf ) break;end; % exit if all either complete or nonterminating
        [a,b] = find(tsig==repmat(min(tsig,[],2),1,NS)); % process that completes quickest
        iTrig(a,1) = b;         % index of which signal triggers first?
        % time at which this signal triggers
        t = tsig(sub2ind(size(tsig),[1:N]',iTrig));
        % response type of this signal triggering: diagonal elements of
        % trigger for the winning signal
        r = p.trigger(sub2ind(size(p.trigger), iTrig,iTrig));
        r(t==inf) = 0;          % if no processes reach threshold, response is 0.
    
        % prepare for next iteration:
        deltaA = repmat(t,1,NS) .* R;
        deltaA(isnan(deltaA)) = 0;    % e.g. infinite time and zero rate --> no rise
        A = A + deltaA; % value of all accumulators after the previous trigger
        newR = p.trigger(iTrig,:);   % for each trial, the effect of the trigger on each signal
        iTrig(t==inf)=nan;
        nochange = isnan(newR) | repmat(t==inf | isnan(t),1,NS); % when not to change trigger
        switch p.mode 
            case SET_RATE
                newR(nochange) = R(nochange); % don't change rate if trigger says nan, or if threshold not crossed for that row
                R = newR;           % rate after trigger
            case MUL_RATE
                newR = R .* newR;
                newR(nochange) = R(nochange);
                R = newR;
            case ADD_RATE
                newR = R + newR;
                newR(nochange) = R(nochange);
                R = newR;
            case SET_ACCU
                A(~nochange) = newR(~nochange);
            case ADD_ACCU
                A(~nochange) = A(~nochange) + newR(~nochange);
        end
         
        % store results of each step
        seqtrig = [seqtrig iTrig];
        Y = Y + (Y==0) .* r;     % change first response if it was zero
        T=T+t;                  % increment cumulative time
        seqt = [seqt T];
        RT = RT + (RT==0) .* (r>0) .* T; % set RT if first response occurred just now
        RTs(sub2ind(size(RTs),[1:N]',r+1)) = T; % set thresh-time of this signal
            % use the +1 hack to allow zeros. the zero column can be
            % deleted later
        if all( t==inf | t<=eps | all(~isnan(RTs),2) | any(isnan(A),2) )  % exit when all signals timed out or responded or crashed
            allEnded=1;
        end
    end
    RTs(isnan(RTs))=inf;
    RTs=RTs(:,[2:end]); % remove first column which corresponds to a zero response
    responseType=bi2de(RTs<inf);                % 1 in each column for response
    responseType(RT>options.trialLength)=0;         % or response is too late: = "No response"
    if(options.censor)
        %RTs(RTs>options.trialLength) = options.trialLength;
        RT(responseType==0) = options.trialLength;  % censoring: just set RT to TRIAL_LENGTH
    else
        %RTs(RTs>options.trialLength) = nan;
        RT(responseType==0) = nan;           % use RT=inf and responseIdx=0.
    end

% end                                         % end parallelise
% allRT=vertcat(RT{:});                       % combine parallelised data
% allResponseType=vertcat(responseType{:});

% nonparallelised output
RT;
allResponseType=responseType;

if isfield(options,'distributionBins') || nargout>2
    distr.pResponse=mean(RTs<options.trialLength,1); % proportion of trials where each response occurs
    NR=size(RTs,2);      % number of response types
    k=zeros(NR,length(options.distributionBins));
    for(j=1:NR)          % for each response type
        k(j,:)=ksdensity(RTs(:,j),options.distributionBins).*distr.pResponse(j); % calculate smooothed distribution
    end
    distr.pRT=k;
end
