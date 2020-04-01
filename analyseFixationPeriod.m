function O = analyseFixationPeriod(A, info, condition, timerange, FFT)
% A ( trial, time ) = eye position as a function of time
% Provide stats on a fixation period, including
%  - microsaccades as a function of time
%  - fourier spectrum
%  - blink rate
% 
% info = structure saccade data in that period
%        read using snipSaccades with 'minsize',0  setting (permit
%        microsaccades.
% 
    
Exclude.blinks    = true; % blinks within the period of interest?
Exclude.saccades  = true; % macro-saccades within period?
Exclude.deviation = true ; % deviation from fixation?

sAmp    = info.sAmpl;
sRT     = info.sRT;
sBlink  = info.sBlink;
sMacro  = sAmp>1*33;            % exclude macrosaccades
sMicro  = sAmp<1*33;            % include microsaccades
sInTime = sRT>timerange(1) & sRT<timerange(2);  % times counted within the fixation period
% maximum deviaiton of eye position from mean position on each trial
% (sometimes missed by macrosaccades!)
maxdev = max(abs(bsxfun(@minus, A(:,timerange(1):timerange(2)),...
                        nanmean(A(:,timerange(1):timerange(2)),2) )),[],2);
% exclude any trial with a blink or a macroscaccade

badTr  = any(sBlink>0 & sInTime,2) | ...
         any(sMacro & sInTime,2)   | ...
         maxdev > 1.8*33;

if ~exist('FFT','var'); FFT   = true; end
PLOT  = false;
STAT  = true;

uc    = unique(condition); % list of conditions
uc(isnan(uc))=[];
NC    = length(uc); % num conditions

for i=1:NC % for each condition
  f = condition==uc(i);          % filter by condition
  pBad(i)=mean(badTr(f));      % how many bad trials in this condition
  f = f & ~badTr;                % remove bad trials from filter
  numOkTrials(i) = sum(f); % number of trials used for this condition
  % for every trial, how many microsaccades?
  nMicro{i} = sum( sMicro(f,:) & sInTime(f,:) ,2 );
  % mean amplitude of microsaccades in the window, for each trial
  ampf = sAmp(f,:);
  % mean amplitude of microsaccades
  aMicro{i} = nanmean( ampf + bool2nan(~sInTime(f,:) | sMacro(f,:)) ,2 ) ;
  % compile all microsaccades from all selected trials, by their RT 
  % ( useful for building a master time histogram of microsaccades )
  theseRTs = sRT( f,: );
  allMicroT{i} = theseRTs( sMicro(f,:) & sInTime(f,:) );
  if FFT
    fx   = nan(diff(timerange)+1,sum(f));
    x    = A(f,timerange(1):timerange(2))'; % position ( time, trial )
    % which time points are 'bad' i.e. high velocity or flat-line 
    % ( usually a sign of interpolation )
    badp = [ true(1,size(x,2)); abs(diff(x))>1.1 | abs(diff(x))==0 ];
    for tr=1:size(x,2) % for each trial
      % extended discrete fourier transform ( Liepin'sh 1996 algorithm that
      % accounts for missing data )
      fx(:,tr) = edft( x(:,tr)+bool2nan(badp(:,tr)) ); 
    end
    fft{i} = fx;
  else
    fft{i} = nan(diff(timerange)+1,sum(f));
  end
  try
    [D{1,i},xi{1,i}] = ksdensity(allMicroT{i},'support',timerange, 'kernel','normal','Bandwidth',0.15);
  catch mexcp
    mexcp
    D{1,i}=nan;
    xi{1,i} = NaN;
  end
end
O.pBad = pBad;
O.numOkTrials = numOkTrials;
O.nMicro = nMicro;
O.aMicro = aMicro;
O.allMicroT = allMicroT;
O.fft = fft;
O.D = nancat(1, D{:})';
O.xi = nancat(1, xi{:})';

O.microDensity = O.D * diff(timerange); % convert density to 'per second'
O.fft = first(nancat([2,3], O.fft)); % freq, trial, condition

if PLOT
  subplot(2,2,1)
  plot(O.microDensity,'xaxisvalues',xi);
  ylabel 'microsaccades / second'
  plotfreq=50:diff(timerange)/2; %  don't report ultra-high or low frequencies
  subplot(2,2,2); 
  errorBarPlot(permute(log(O.fft(plotfreq,:,:)),[2,1,3]),'area',1, 'xaxisvalues',plotfreq);
  ylabel 'log power'; xlabel 'frequency';
  drawnow
end
if STAT
  % compute effect of condition, for each time point
  b = ols( (reshape(permute(O.fft,[3,1,2]),NC,[])),[zscore([1:NC]'),  ones(NC,1) ],[1 0]);  
  O.fft_beta = b; 
end

%% do velocity

sDur = info.sDur; % duration
t = repmat(timerange(1):timerange(2), [size(sRT,1),1]);
for i = 1:size(sRT,2)
    t1(:,:,i) = t > sRT(:,i) & t < (sRT(:,i)+sDur(:,i));
end
saccInds = any(t1,3);
A2 = A;
A2(saccInds) = NaN;

speed = abs(diff(A2, [],2));

% estimate velocity over 40ms moving window
meanSpeed = movmean(speed, 40, 2);

% mean over trials
O.meanSpeed = meanSpeed;
O.meanSpeedCond = permute( groupMeans(O.meanSpeed,1,condition,'dim'), [1,3,2]);