function reward = adaptiveReward( t, parameters )
% reward = adaptiveReward( t )
%   calculate an adaptive reward given a reaction time
%   reward increases with shorter reaction times
%   but adapts to subjects' mean speed over last 20 trials
% reward = adaptiveReward( t , thread )
%   use specified thread (default 1) for updating rewards.
%
% use adaptiveReward('Reset') to reset the RT staircase
%
% use 'GetParameters' and 'SetParameters' to get the parameters:
%   for one thread,   [historyLength, increment, maxPercentile, minPercentile, tMin, tMax]
%   for many threads, { [historyLengths], [increments], [maxPercentiles],
%                       [minPercentiles], [tMins], [tMaxs]}
%
% history length: how many previous trials to average
% increment:      proportion to change the limits by when subject
%                 improves/worsens
% minPercentile:  the desired proportion of RTs giving maximum reward = 1
% maxPercentile:  the desired proportion of RTs where reward falls 
%                 below = 1/e = 0.368
% tMin:           the RT currently at the minimum percentile
% tMax:           the RT currently at the maximum percentile
%
% e.g. if minPercentile = 0.1, about 10% of rewards will be at ceiling
%  and if maxPercentile = 0.6, about 60% of rewards will be below 37% maximal
%  so for reaction times, the following is an example
%  adaptiveReward('SetParameters', [20, 0.05, 0.80, 0.10, 0.350, 0.420]);

global threads

%% parameters

if strcmpi(t,'SetParameters') % set parameters for all threads
  if(exist('parameters','var'))
      parameters = num2cell(parameters);
    if isnumeric(parameters)
       [historyLength, increment, maxPercentile, minPercentile tMin tMax]=items(parameters);
    elseif iscell(parameters)
      [historyLength, increment, maxPercentile, minPercentile tMin tMax]=deal(parameters{:}); 
    end
    return;
  end
  error('syntax: adaptiveRewardThreads(''SetParameters'',{[],[]...})');
end

if strcmpi(t,'Reset') % reset the thread's memory (but retain parameters)
  if exist('parameters','var') % if the thread was specified
    resetThreads = parameters;
    return;
  else % no thread specified: reset all
    resetThreads = 1:length(threads);
  end
  for i=1:length(resetThreads)
    threads(resetThreads(i)).tRecent = [];
    threads(resetThreads(i)).tMin    = 0.300;
    threads(resetThreads(i)).tMax    = 0.350;
  end
  return;
end;


% which thread are we on?
if isnumeric(parameters) && length(parameters)==1
  thread = parameters;
else
  thread = 1;
end

% ensure the thread is initialised
if length(threads) >= thread && ~isempty(threads(thread).historyLength)  % already initialised: get parameters
  h = threads(thread); 
elseif length(threads)>0 % if one thread is initialised, copy its parameters
  for i=1:length(threads) % find first thread that has been initialised
    if ~isempty(threads(i).historyLength);
      break
    end
  end
  h = struct();
  h.historyLength = threads(i).historyLength;
  h.increment     = threads(i).increment;
  h.minPercentile = threads(i).minPercentile;
  h.maxPercentile = threads(i).maxPercentile;
  h.tRecent       = [];
  h.tMin          = threads(i).tMin;
  h.tMax          = threads(i).tMax;
%   h.thread_number = threads(1).thread_number;
else % no threads initialised: create new default parameters 
  h = struct();
  h.historyLength = 20;
  h.increment     = 0.05;
  h.minPercentile = 0.10;
  h.maxPercentile = 0.70;
  h.tRecent       = [];
  h.tMin          = 0.300;
  h.tMax          = 0.350;
  threads         = h; % initialise threads to be a structure
end

% get parameters for a specifc thread
if strcmpi(t,'GetParameters')
  reward=[h.historyLength h.increment h.maxPercentile h.minPercentile h.tMin h.tMax h.sequenceLength];
  return;
end;

%% calculate reward, using current thread 'h'
n=min(h.historyLength,length(h.tRecent));
if(n>0) 
% calculate new reward parameters
    pmax=sum(h.tRecent>=h.tMax)/n;
    pmin=sum(h.tRecent<h.tMin)/n;
    dma=1; dmi=1;
    if(pmax>h.maxPercentile) dma=1+h.increment;
    elseif (pmax<h.maxPercentile) dma=1/(1+h.increment);
    end
    if(pmin>h.minPercentile) dmi=1/(1+h.increment);
    elseif (pmin<h.minPercentile) dmi=(1+h.increment);
    end;
    newTMin=h.tMin + (dmi-1) * (h.tMax-h.tMin);
    newTMax=h.tMax + (dma-1) * (h.tMax-h.tMin);
    if(newTMin<newTMax)
        h.tMin=newTMin;
        h.tMax=newTMax;
    else
        if(newTMax>h.tMax) h.tMax=newTMax; end; % continue changing if
        if(newTMin<h.tMin) h.tMin=newTMin; end; % they are moving apart
        % or if both moving towards each other, don't replace either;
    end;

    reward=exp(-max(0,(t-h.tMin))/(h.tMax-h.tMin));
else % not enough samples yet to do calculate the reward
    reward = 0.5; 
end;
% keep track of history of RTs
h.tRecent=[t h.tRecent(1:n)]; 
% store new thread state
threads(thread)=h;

