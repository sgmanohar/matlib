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
%      [historyLength, increment, maxPercentile, minPercentile, tMin, tMax]
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

global tRecent tMin tMax
global historyLength increment maxPercentile minPercentile
if(strcmpi(t,'Reset'))
    tRecent=[]; tMin=0.300; tMax=0.350;
    clear historyLength
    return;
end;
if(isempty(tMin)) tMin=0.300; tMax=0.350; end;

%% parameters
if(exist('parameters','var'))
  [historyLength, increment, maxPercentile, minPercentile tMin tMax]=items(parameters);
  return;
elseif(~exist('historyLength','var') || length(historyLength)~=1 || historyLength<1)
  historyLength=20;
  increment=0.05;
  minPercentile=0.10;
  maxPercentile=0.70;
end;
if(strcmpi(t,'GetParameters'))
  reward=[historyLength increment maxPercentile minPercentile tMin tMax];
  return;
end;

%% calculate reward
n=min(historyLength,length(tRecent));
if(n>0) 
% calculate new reward parameters
    pmax=sum(tRecent>=tMax)/n;
    pmin=sum(tRecent<tMin)/n;
    dma=1; dmi=1;
    if(pmax>maxPercentile) dma=1+increment;
    elseif (pmax<maxPercentile) dma=1/(1+increment);
    end
    if(pmin>minPercentile) dmi=1/(1+increment);
    elseif (pmin<minPercentile) dmi=(1+increment);
    end;
    newTMin=tMin + (dmi-1) * (tMax-tMin);
    newTMax=tMax + (dma-1) * (tMax-tMin);
    if(newTMin<newTMax)
        tMin=newTMin;
        tMax=newTMax;
    else
        if(newTMax>tMax) tMax=newTMax; end; % continue changing if
        if(newTMin<tMin) tMin=newTMin; end; % they are moving apart
        % or if both moving towards each other, don't replace either;
    end;

    reward=exp(-max(0,(t-tMin))/(tMax-tMin));
else
    reward = 0.5; 
end;
tRecent=[t tRecent(1:n)]; %history of RTs



function varargout=items(X)
for(i=1:length(X)) varargout{i}=X(i); end;
