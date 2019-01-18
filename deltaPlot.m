function [delta, bin, hplot, t_test_result, p_value, t_statistics, t_threshold] = ...
                     deltaPlot(DATA, varargin)
% function [delta, mean_bincentre, hplot,  
%           ttest_result, t_statistics, t_thresholds ]
%             = deltaPlots(DATA, [params...] )
%
% Create a delta plot (usually of reaction times) of the differences
% between conditions, at each quantile of the distribution.
% it takes Condition1 - Condition2, Condition2 - Condition3 etc.
% 
% DATA ( SUBJECT, CONDITION, TRIAL )
%  OR
% DATA { SUBJECT, CONDITION } ( TRIAL )
%   
% Creates delta plots averaged across subjects.
% 
% for a simple delta plot there should be 2 conditions, i.e. 
%   size(DATA) = subjects x 2 x trials
% if there are more than 2 conditions, each neighbouring pair is compared
% i.e. condition1-condition2, condition2-condition3, etc.
% 
% If the different subjects / conditions have different numbers of trials,
% you can either use nan-padding (e.g. using nancat) or use the cell-array
% version of DATA.
%
% RESULT:
% 
% delta          = the value of the difference between conditions, for
%                  each bin. DELTA ( SUBJECT, CONDITION, BIN )
% mean_bincentre = the value at the centre of each of the bins
% 
% if you request 'hplot' as ouput, it will plot the delta plot with dotted
% lines for error bars.
% 
% if you request 't_test_result', then a permutation-based test will
% calculate the threshold for t, corrected for FDR of multiple comparisons 
% across the bins, and will return 0 or 1 for each quantile. The null
% distribution is when condition is randomised - i.e. the
% permutation randomises the order of the two conditions for each subject.
% 
% Delta plots - what do they mean?
% (Ridderinkhof 2002, Ratcliffe 1979)
%
% SGM 2014

% Default 100 bins, gaussian smoothing=0.1 of range, p<0.05

NB     = 100;            % number of quantile bins
SMOOTH = floor(NB/5);    % smoothing window (as a number of bins). < 2 means no smoothing
SMOOTH_FILTER = 'loess'; % which smoothing function to use
EXTRASMOOTHING = 5;      % additional smoothing just before plotting. Does not apply to returned values
ALPHA  = 0.05;           % if t-test is requested, then 
TWO_TAILED = false;      % should I test if condition1 not equal to condition2?
WIDTH  = 0.2 ;           % width (in quantiles) of bins. Can be zero for a pure quantiled plot
PLOT   = true;           % draw a delta plot
PLOT_P = false;          % draw boxes for where p < alpha

if exist('parsepvpairs','file') % attempt to read in parameters
  [  NB,  WIDTH,  ALPHA,  SMOOTH,  TWO_TAILED,   PLOTARGS,  PLOT_P     PLOT] = parsepvpairs( ...
   {'NB','WIDTH','ALPHA','SMOOTH','TWO_TAILED', 'PLOTARGS', 'PLOT_P', 'PLOT'}, ...
   {100 , 0.2 ,    0.05,  1,       false      ,  {} ,       false  , true }, varargin{:});
elseif nargin>1, warning('parsepvpairs.m not found; parameters ignored!'); end


WIDTH=floor(WIDTH*NB); % convert width from a quantile fraction into a number of bins
if WIDTH>4 & mod(WIDTH,2), WIDTH=WIDTH-1; end % make it even - so that bin centre is even
if isnumeric(DATA)     % numeric array? convert to cells
  if ndims(DATA)<3, error('deltaplot:dimension','DATA should be an 3-dimensional array'); end
  for i=1:size(DATA,1) % each subject
    for j=1:size(DATA,2) % each condition
      DATA2{i,j}=squeeze(DATA(i,j,:)); 
    end % next condition
  end % next subject
  DATA=DATA2; % now it's DATA { SUBJECT, CONDITION } ( TRIAL )
end
NSubj = size(DATA,1); % number of subjects
NCond = size(DATA,2); % number of conditions

flattening = false;   % have we had to flatten any cells into columns (warn if so)
quantiles = linspace(0,1,NB+1);  % the actual quantiles to use
delta = nan(NSubj,NB-WIDTH, NCond-1); % create empty matrix for results: Delta between conditions
bin   = nan(NSubj,NB-WIDTH, NCond-1); % and this is for the abscissa.
for subject = 1:NSubj % for each subject
  for condition = 1:(NCond-1) % for each pair of neighbouring conditions
    c1 = DATA{subject,condition};       % get data for 1 subject for 2 conditions
    c2 = DATA{subject,condition+1};
    if ~isvector(c1) || ~isvector(c2),  % flatten to column vector if needed
      flattening=1; c1=flat(c1); c2=flat(c2); 
    end
    q1 = quantile(c1, quantiles);   % create a vector of quantiles for each condition
    q2 = quantile(c2, quantiles);   % quantile ignores nans.
    %if WIDTH == 0
    %  delta(subject,:,condition) = q1-q2;     % DELTA ( SUBJECT, QUANTILE, CONDITION_PAIR )
    %else
      for i=1:(NB-WIDTH) % for each quantile bin
        meanc1 = nanmean( c1( c1>=q1(i) & c1<q1(i+WIDTH+1) ) );  % c1 values within quantile range 
        meanc2 = nanmean( c2( c2>=q2(i) & c2<q2(i+WIDTH+1) ) );  % c2 
        delta(subject,i,condition) = meanc1-meanc2;
      end
    %end
    halfwidth=floor(WIDTH/2);
    % BIN   ( SUBJECT, QUANTILE, CONDITION_PAIR )
    bin(subject,:,condition)   = (q1(halfwidth+1:end-halfwidth-1)+q2(halfwidth+1:end-halfwidth-1))/2; 
  end % next condition pair
end % next subject
if flattening, warning('deltaplot:flatten','Flattening matrix data into columns'); end

mbin=squeeze(nanmean(bin));    % average across subjects so can be plotted on single graph
                               % and remove first dimension, to give 
                               % X ( QUANTILE, CONDITION_PAIR )
% take mean delta across subjects in each quantile bin
mdelta = permute(nanmean(delta),[2 3 1]);                 
% and calculate sd across subjects at each quantile bin
edelta = permute(( nanstd(delta)) / sqrt(NSubj),[2 3 1]); 


if SMOOTH>1 % SMOOTH
  for j=1:NCond-1 % for each pair of conditions (smooth requires single column)
    mdelta(:,j) = smooth(mdelta(:,j) ,SMOOTH,SMOOTH_FILTER);
    edelta(:,j) = smooth(edelta(:,j) ,SMOOTH,SMOOTH_FILTER);
  end
  delta=smoothn(2, delta, SMOOTH, SMOOTH_FILTER)
end

if nargout>2 && PLOT % requested plot handle? then PLOT GRAPH
  if NB<7     % this version with error bars is good for few bins
    hplot=errorBarPlot(delta, 'xaxisvalues',mbin, 'plotargs', PLOTARGS);
  else        % this version shows a curve - is good for many bins!
    washeld=ishold();
    plotArgs = PLOTARGS;
    if 0 % use dotted lines as error bars -- useful for exporting 
      hp1=plot(mbin, mdelta, plotArgs {:}, 'Marker','.'); % plot mean delta for each bin
      hold on
      hp2=plot(mbin, mdelta+edelta, plotArgs {:}, 'LineStyle',':'); % error lines above
      hp3=plot(mbin, mdelta-edelta, plotArgs {:}, 'LineStyle',':'); % and below
      hplot=[hp1 hp2 hp3];
    else % use Area errorBarPlot
      hplot=errorBarPlot(delta, 'area',1, 'xaxisvalues',mbin, 'plotargs', PLOTARGS, 'smooth',EXTRASMOOTHING );
      hold on;
    end
    plot(xlim,[0 0],':'); % dotted zero-line = no difference between conditions
    if ~washeld, hold off; end
    %ylabel(['\Delta']);
    xlabel('value');
    % create legend = "condition1-condition2" etc.
    labels = arrayfun(@(x) sprintf('condition%g-condition%g',x,x+1), [1:(NCond-1)]','uniform',0);
    % insert blanks for error-lines
    legend(flat([labels, repmat({'',''},NCond-1,1)]));
    % hplot=gca; % return axes handle
  end
end

if nargout>3 % request t statistic
  % perform a permutation test to calculate the FDR across all bins
  for i=1:5000 % iterate 5000 times
    % what if the order of cond1 and cond2 were randomly chosen?
    randbool = (rand(NSubj,1)>0.5)*2-1; % +1 or -1 for each subject
    permuted_delta  = bsxfun(@times,randbool,delta); % delta with random swapping of cond1/cond2
    % t statistic for this permutation = mean / stderr_of_mean
    % PERMUTED_TSTATS ( QUANTILE, CONDITION_PAIR )
    permuted_tstats = permute(mean(permuted_delta) ./ std(permuted_delta), [2 3 1]) / sqrt(NSubj); 
    % MAXT ( ITERATION, CONDITION_PAIR )
    maxt(i,:) = max(permuted_tstats, [], 1); % maximum t across all quantiles 
                                             % (for each condition-pair separately)
    
    % this following line is for if you wanted
    % a two-tailed test, but at the moment I am assuming the hypothesis is
    % that condition1 is always bigger than condition2, so it's a 1-tailed
    % test.
    if TWO_TAILED
      mint(i,:) = min(permuted_tstats, [], 2); 
    end
  end
  t_statistics  = squeeze(mean(delta) ./ std(delta)) / sqrt(NSubj); % actual t statistic!
  p_value_pos   = 1-mean( bsxfun(@gt, t_statistics, maxt ) ); % proportion of null distribution above max t
  if TWO_TAILED
    t_threshold_lower = quantile(mint, ALPHA/2);   % find upper 2.5% of max-t values
    t_threshold_upper = quantile(maxt, 1-ALPHA/2); % find lower 2.5% of min-t values
    % null hypothesis rejected at 5% if either threshold is passed.
    t_test_result = bsxfun(@gt, t_statistics, t_threshold_upper)  ...
                 |  bsxfun(@gt, t_statistics, t_threshold_lower); 
    p_value_neg   = 1-mean( bsxfun(@lt, t_statistics, mint ) ); % proportion of null below min above t
    p_value = (p_value_neg<0.5).*p_value_neg + (p_value_neg>=0.5).*p_value_pos;
  else
    t_threshold   = quantile(maxt,1-ALPHA);     % upper 5% of max-t values
    % null hypothesisi rejected at 5% if t_statistics > t_threshold.
    t_test_result = bsxfun(@gt, t_statistics, t_threshold); 
    p_value = p_value_pos;
  end
  if PLOT_P % PLOT the t test results? This shows the significant points where
    % delta is nonzero, below the graph as a bar. 
    % Note: this might upset adding any more subplots to the figure!
    set(gca,'color','none'); mainaxes = gca;
    rect=get(mainaxes,'position'); % select region in bottom 5% of axis
    axes('position',[rect(1) rect(2) rect(3),rect(4)*0.05]);
    imagesc(t_test_result); axis off; 
    axes(mainaxes); % revert to main axes
  end
end
