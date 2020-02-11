function [Xb,Yb, p,t, h, resid ]=conditionalPlot(X,Y, Nbins, varargin)
% [Xb,Yb ]                 = conditionalPlot(X,Y, [nbins], ...)
% [Xb,Yb, p, t, h, resid ] = conditionalPlot(X,Y, nbins,  ...)
%      X{ SUBJECT, CONDITION } ( TRIAL )
% or   X( TRIAL, SUBJECT )            -- for one condition
% or   X( TRIAL, SUBJECT, CONDITION ) -- for multiple conditions
% or   { Table, 'xvar', 'yvar', 'subject', 'condition' }
%
% Takes sliding quantile bins of X, and calculates the mean Y for 
% each bin. Nbins indicates the width of the window (5 means 20 percentiles)
% This is done for each subject and condition separately. 
% The bin values of Xb and Yb are plotted as mean and sem across subjects.
% Different conditions are plotted as separate lines.
%
% returns: Xb ( SUBJECT, CONDITION, QUANTILE ) = median of X for each bin
%          Yb ( SUBJECT, CONDITION, QUANTILE ) = mean   of Y for each bin
%          p ( QUANTILE ) = probability that conditions are significantly 
%             different at a quantile, using permutation test.
%             Within-subject comparison between conditions is performed by 
%             permutation, and multiple conditions are treated as a 
%             single continuous linear predictor (1,2,3 etc).
%          T ( QUANTILE )  = t statistic for linear effect of condition
%          H = handle of lines as returned by errorbarplot. 
%  'resid' - returns the residuals after the means of the X-bins surrounding 
%     each Y-data point have been subtracted. It is the Y-distance from the 
%     binned curve. 
%
% Parameters:
%  'meanfunc' = function handle to use instead of "nanmean". 
%  'xtransform', 'ytransform' - set transformations for the plot axes. Can
%     be a @function, or one of 'identity' (default), 'probit', 'recip',
%     'log', 'asin', 'gauss' (inverse gaussian from 0 to 1), or
%     'gauss2' (inverse gaussian from 0.5 to 1). Note that the output values
%     Xb and Yb will not be transformed.
%  'divisions' - number of subdivisions, default 100. If divisions == Nbins,
%      then there is no overlap 
%  'color' - set colour of plot. If omitted, the ColorOrder is used for
%     different conditions
%  'withinsubjecterror' - if 0, errorbars are sem across subjects (for comparing subjects). 
%     if 1, then a subject's values are mean-corrected before sem is taken (for comparing conditions).
%     if 2, then the conditions are individually mean-corrected
%     ('within-condition error', for comparing time points within condition)
%  'smooth' - 0 = no smoothing. n = smooth over n bins, using smoothn (defaults
%     to gaussian kernel.)
%  'plotIndividuals' - if true, then draw a line for each person.
%  'plotTrajectories' - it true, then draw lines linking the conditions. 
%  'area' - 0 = dotted lines, 1 = area, >3 = graded quantile fill
%  'plotstats' - if true, compare the lines using permutationOLS
%  remaining unrecognised arguments are passed on to the 'plot' command.
% 
% dependencies: smoothn, nanmean, nancat, removeIdenticals, permutationOLS
% sanjay g manohar 2013


%%%%%%%%%%%%%%%%%%%%%
% Setup Parameters

% handle tables as inputs
if iscell(X) & numel(X)<6 & istable(X{1}) % table form
  if numel(X)==5 % {T, 'xvar','yvar','sub','lines'}
    if isempty(X{4}), % if 'sub' is empty, then treat it as one subject
      X{4}='CONST';  X{1}.CONST=ones(numel(X{1}(:,1)),1);  
    end
    sub  = X{1}.(X{4}); 
    con  = X{1}.(X{5});
    usub = unique(sub);
    ucon = unique(con);
    for i=1:length(usub)
      for j=1:length(ucon)
        tX{i,j} = X{1}.(X{2})(sub==usub(i) & con==ucon(j));
        tY{i,j} = X{1}.(X{3})(sub==usub(i) & con==ucon(j));
      end
    end
  elseif numel(X)==4 % {T,X,Y,G}
    if isempty(X{4}), % if 'sub' is empty, then treat it as one subject
      X{4}='CONST';  X{1}.CONST=ones(length(X{1}),1);  
    end
    sub = X{1}.(X{4}); 
    usub =unique(sub);
    for i=1:length(usub)
      tX{i,1} = X{1}.(X{2})( sub==usub(i) );
      tY{i,1} = X{1}.(X{3})( sub==usub(i) );
    end    
  elseif numel(X)==3 % {T, X, Y}
    tX = X{1}.(X{2});
    tY = X{1}.(X{3});
  elseif numel(X)==2, error('please specify two variables')
  elseif numel(X)==1, error('please specify two variables')
  end
  if exist('Y','var') % shift the other parmeters along by one, since Y is already done
    if exist('Nbins','var')
      varargin=[{Nbins} varargin];
    end
    Nbins = Y;
  end
  labels=X(2:end); % keep values in table form
  X=tX;Y=tY;
end % table version

if(~exist('Nbins','var') || isempty(Nbins)) Nbins=5; end;


wascells = iscell(X); % keep track of what format input was
if ~iscell(X)
  if ndims(X)==2 % treat as X ( TRIAL, SUBJECT )
    if size(X,1)==1, X=X'; end
    if size(Y,1)==1, Y=Y'; end
    X  = mat2cell(X, size(X,1),   ones(1,size(X,2))   )';  % break into cells
    Y  = mat2cell(Y, size(Y,1),   ones(1,size(Y,2)) )';
  elseif ndims(X)==3 % treat as X ( TRIAL, SUBJECT, CONDITION ): convert to cell { SUBJECT, CONDITION } (TRIAL)
    X  = sq(mat2cell(X, size(X,1),   ones(1,size(X,2)),   ones(size(X,3),1) )); 
    Y  = sq(mat2cell(Y, size(Y,1),   ones(1,size(Y,2)),   ones(size(Y,3),1) ));
  end
end

Nq       = 100;     % number of quantiles 
binwidth = 1/Nbins; 
if exist('nanmean','file')  MEANFUNC = @nanmean;  % do we have 'nanmean' ?
else                        MEANFUNC = @mean;     % if not, use mean.
end
STDFUNC            = @(x)nanstd(x); % error bar size when only one subject given
doErrorBars        = true;  % on the non-overlapping binning plot, show error bars for both X and Y?
FILL_AREA          = 2;  % if false, show dotted lines for standard error, rather than shaded area
REMOVE_IDENTICAL_X = true;  % if there aren't enough unique X-values, then interpolate?
WITHINSUBJECTERROR = true;  % remove subject means before calculating standard error?
PLOT_INDIVIDUALS   = false; % draw a line for each subject, in addition to mean and error.
STANDARDISE_RESIDUALS = false;
PLOT_TRAJECTORIES  = false; % draw lines between conditions, for each quantile
% transform axes?
T.identity = @(x)x;
T.probit   = @(x) log(eps+x./(1-x+eps));
T.recip    = @(x) -1./x;     % reciprocal (eg for 'reciprobit'  
T.log      = @(x) log(x+eps);
T.asin     = @(x) asin(sqrt(x)); 
T.gauss    = @(x) erfinv( max(eps-1, min(1-eps, (x*2-1))) ); % gaussian error, range [0-1]
T.gauss2   = @(x) erfinv( max(eps-1, min(1-eps, (x-0.5)*4-1)) ); % gaussian in range [0.5-1]
% these two lines apply transforms to the data before plotting. 
% the function output values are unchanged by the transform.
YTRANSFORM = T.identity;  
XTRANSFORM = T.identity;
i=find(strcmpi(varargin, 'xtransform')); if i, if isstr(varargin{i+1}) XTRANSFORM=T.(varargin{i+1}); 
  else XTRANSFORM=varargin{i+1};  end ; varargin([i i+1])=[];end

i=find(strcmpi(varargin, 'xtransform')); if i, if isstr(varargin{i+1}) YTRANSFORM=T.(varargin{i+1}); 
  else YTRANSFORM=varargin{i+1};  end ; varargin([i i+1])=[];end

i=find(strcmpi(varargin, 'meanfunc')); if i, 
  MEANFUNC=varargin{i+1}; varargin([i,i+1])=[];
end

i=find(strcmpi(varargin, 'divisions')); if i, 
  Nq = varargin{i+1}; varargin([i,i+1])=[];
end
% use nonoverlapping bins? 
% advantage is the values are independent, and thus stats are simpler.
% if not, use overlapping bins, which must be compared with permutation test 
NO_OVERLAP = Nq*binwidth == 1; 
STANDARD_BINNING = false;
i=find(strcmpi(varargin, 'standardbinning')); if i, 
  STANDARD_BINNING = varargin{i+1}; varargin([i i+1])=[]; 
end



i=find(strcmpi(varargin, 'WithinSubjectError')); if i, 
  WITHINSUBJECTERROR = varargin{i+1}; varargin([i i+1])=[]; 
end
i=find(strcmpi(varargin, 'Smooth')); if i, 
  SMOOTH = varargin{i+1}; varargin([i i+1])=[]; 
end
i=find(strcmpi(varargin, 'PlotIndividuals')); if i, 
  PLOT_INDIVIDUALS = varargin{i+1}; varargin([i i+1])=[]; 
end
i=find(strcmpi(varargin, 'PlotTrajectories')); if i, 
  PLOT_TRAJECTORIES = varargin{i+1}; varargin([i i+1])=[]; 
end
i=find(strcmpi(varargin, 'Conditions')); if i, 
  C = varargin{i+1}; varargin([i i+1])=[]; 
  if ~isempty(C)
    if size(X,2)>1, error('cannot specify conditions since X already contains conditions'); end
    if isnumeric(C) && ndims(C)==2, C=mat2cell(C, size(C,1), ones(1,size(C,2)) )'; end
    uc=unique(nancat(2,C{:})); uc(isnan(uc))=[]; 
    for i=1:size(X,1) % subject
      for j=1:length(uc) % condition
        X_new{i,j}  =X{i}(  C{i}==uc(j) );
        Y_new{i,j}  =Y{i}(  C{i}==uc(j) );
      end
    end
    X=X_new; Y=Y_new;
  end
end
i=find(strcmpi(varargin, 'Area')); if i
  FILL_AREA = varargin{i+1}; varargin([i i+1])=[];
end
i=find(strcmpi(varargin, 'plotargs')); if i
  plotargs = varargin{i+1}; varargin([i i+1])=[];
else  plotargs={};
end

DO_STATS = nargout>2; % whether to do statistics?
DO_STATS_PLOT = nargout > 4; 
i=find(strcmpi(varargin,'plotstats')); if i
  DO_STATS_PLOT = varargin{i+1}; varargin([i i+1])=[]; 
end

if NO_OVERLAP , N=Nbins;   % N is the number of quantiles to sample, i.e.
else N=floor(Nq*(1-binwidth));  % number of points on a CAF curve, e.g. 100*0.8
end
Xb=nan(size(X,1), size(X,2), N); % create blank output variable (needed for parallel processing)
binwidth_n = floor( binwidth*Nq );  % points per window-bin
RESIDUALS = nargout>5; % should we calculate residuals of X?
if ~exist('SMOOTH','var'), SMOOTH   = binwidth_n;    end  % smoothing width = window size(used when overlap) 

i=find(strcmpi(varargin, 'standardise')); if i, 
  STANDARDISE_RESIDUALS= varargin{i+1}; varargin([i,i+1])=[];
end


%%%%%%%%%%%%%%%%%%%%%
% Calculate binning

FULL_RANGE = true;  % if false, exclude 1% of highest and lowest x-values.
%%%% change the outer 'for' to 'parfor' for faster execution.
for(i=1:size(X,1))       % for each subject
  x_i=nan(size(X,2),N);  % keep track of bin centres
  for(j=1:size(X,2))     % for each condition
    if ~STANDARD_BINNING      % VARIABLE WIDTH sliding window bins BY QUANTILE
      quant_levels = linspace(0,1,Nq);
      quantiles = quantile(X{i,j}, quant_levels); % create quantiles
      % go through bins, but there aren't 100 bins because of the width!
      for q=1:N   % e.g. 1:80 for 5 bins
        % get start, mid and end of the bin
        if FULL_RANGE % uses whole range of data 
          range=quantiles( q + [0,floor(Nq*binwidth/2), floor(Nq*binwidth)] );
        else % remove 1% highest and lowest x-data 
          range=quantiles( q + [1,floor(Nq*binwidth/2), floor(Nq*binwidth)-1] );
        end
        % filter data within range - keep symmetrical top and bottom edges
        % (otherwise gives biases)
        f = X{i,j} > range(1)  &  X{i,j} < range(3); 
        % Yb ( SUBJECT, CONDITION, QUANTILE )
        Yb(   i,j,q) = MEANFUNC( Y{i,j}(f) );      % mean of Y in each bin
        Yb_sd(i,j,q) = STDFUNC(  Y{i,j}(f) ) / sqrt(sum(f)) ; % Y standard error for each bin
        % Xb values will be the centre (median) RT of the bin.
        % Change this to range([1 3])/2 to get the midway-point between edges
        x_i(j,q)     = range(2);                       % centre (median) of the bin
        x_m(j,q)     = nanmean( X{i,j}(f) );             % centre (mean)   Xval of the bin
        x_sd(j,q)    = nanstd(  X{i,j}(f) ) / sqrt(sum(f)); % standard error of X for bin
      end
      % 8/7/14 added this to prevent nonunique bin centres
      % note that this will fail if X does not have enough unique values to
      % make bins. (and warning given)
      if REMOVE_IDENTICAL_X
        x_i(j,:)=removeIdenticals(flat(sq(x_i(j,:))));
      end
    else  % STANDARD_BINNING:  FIXED WIDTH BINS
      edges =  linspace(min(X{i,j}), max(X{i,j}),Nq+1);  % divide range equally
      edges = [-inf, edges(2:end-1), +inf];              % end bins are infinite
      for q=1:N % bins
        range        = [edges(q), edges(q+binwidth_n)];
        f            = X{i,j} > range(1)  &  X{i,j} <= range(2);
        Yb(   i,j,q) = MEANFUNC(  Y{i,j}(f) ) ;
        Yb_sd(i,j,q) = STDFUNC(   Y{i,j}(f) ) / sqrt(sum(f)) ; % Y error
        x_i(j,q)     = nanmedian( X{i,j}(f) );
        x_m(j,q)     = nanmean(   X{i,j}(f) );
        x_sd(j,q)    = nanstd(    X{i,j}(f) ) / sqrt(sum(f)); % standard error of X for bin
      end
    end
  end % next condition
  Xb(i,:,:) = x_i;     % store bins
end % next subject
if isempty(Xb) || ~exist('Yb','var')
  warning('no data for conditionalPlot');
  return;
end

%%%%%%%%%%%%%%%%%%%%%%
% Calculate residuals for each datapoint (26/2//16)

% subtracts Y from the mean Y in a window around the X coordinate.  The 
% window width is quantiled, i.e. each window has a given number of 
% datapoints in it, except at the extreme values of X, which are calculated 
% with smaller windows. 
% This effectively de-trends data, with a window size equal to the number 
% of data points divided by the number of bins. 
if RESIDUALS                  % calculating residuals? if so, create space to store them.
  for i=1:size(X,1);          % subjects
   for j=1:size(X,2)          % conditions
     NP = length(X{i,j});     % number of data points in this condition.
     window = floor(NP / Nbins / 2); % half-width (number of datapoints) of window around current X value
     [sx,ord] = sort(X{i,j}); % Sort data in this condition
     sy = Y{i,j}(ord);        % in order of ascending X.
     % this is a slow way of doing it, point by point, but I cant think of 
     % anything faster right now...
     for k=1:NP               % for each datapoint in the condition:
       left  = max(1,k-window);           % bin edges for window around 
       right = min(NP,k+window);          % the current data point
       % calculate residual as the difference from mean of Y in the window
       resid{i,j}(k,1) = Y{i,j}(k) - MEANFUNC(Y{i,j}(left:right)); 
       if STANDARDISE_RESIDUALS
         r_sd = nanstd( Y{i,j}(left:right) );
         resid{i,j}(k,1) = resid{i,j}(k,1) / r_sd ;
       end
     end
   end
  end
end
   

%%%%%%%%%%%%%%%%%%%%%%
% Now do the plotting

oldhold=ishold();             % keep track of whether "hold on"
t_x = XTRANSFORM( Xb );       % transform the data for plotting
t_y = YTRANSFORM( Yb );
h=[];                         % keep graphics / line handles
cols = get(gca,'ColorOrder'); % use axis colours for different lines

if NO_OVERLAP % for non-overlaping bins, we have just a couple of bins
  for j=1:size(X,2) % so for each bin, 
    if doErrorBars  % draw error-bars for both X and Y
      if size(X,1)>1  % num subjects > 1
        t_xn = t_x; % subtract off means to normalise error bars
        t_yn = t_y; % for within-subject effect
        if WITHINSUBJECTERROR == 1 % for whole subject's data
          t_yn = bsxfun( @minus, t_yn , nanmean(nanmean(t_yn,2),3) );
        elseif WITHINSUBJECTERROR == 2 % for each condition separately
          t_yn = bsxfun( @minus, t_yn , nanmean(t_yn,3) );
        end
        % calculate standard error in x and y
        ex = sq(nanstd(t_xn(:,j,:),[],1))/sqrt(size(Xb,1)) ;
        ey = sq(nanstd(t_yn(:,j,:),[],1))/sqrt(size(Yb,1)) ;
        errorbarxy(  XTRANSFORM(  sq(nanmean(t_x(:,j,:),1))  ), ...
                     YTRANSFORM(  sq(nanmean(t_y(:,j,:),1))  ), ...
          ex,ey ...
          ,[],[], cols(j,:)  , .2*[1 1 1]);
      else % only one subject
        if ~isequal(XTRANSFORM, T.identity) || ~isequal(YTRANSFORM, T.identity)
          warning('transforms not applied'); 
        end
        errorbarxy(  XTRANSFORM(  sq(t_x(:,j,:))  ), YTRANSFORM( sq(t_y(:,j,:)) ), ...
          sq(x_sd), sq(Yb_sd) ... 
          ,[],[], cols(j,:)  , .2*[1 1 1]);        
      end
    else            % or just plot the means
      h2=plot(sq(nanmean(t_x(:,j,:),1)), ...
              sq(nanmean(t_y(:,j,:),1)) ...
        , cols(j,:) );
      h=[h h2];
    end
    hold on
  end
else % SLIDING BIN = we have a whole curve to plot
  % Yb ( SUBJECT, CONDITION, TIME )
  for j=1:size(Xb,2)       % for each condition
    x = sq(nanmean(t_x(:,j,:),1)); line = sq(nanmean(t_y(:,j,:),1));
    if WITHINSUBJECTERROR==1  % 2016/1: subtract off subject average across conditions and time
        demeaned_t_y  = bsxfun(@minus, t_y, nanmean(nanmean(t_y,2),3));
    elseif WITHINSUBJECTERROR==2  % subtract mean for each condition and subject
        demeaned_t_y  = bsxfun(@minus, t_y, nanmean(t_y,3));
    else
      demeaned_t_y = t_y; % don't normalise - use between-subjects sem.
    end
    if size(demeaned_t_y,1)==1 % if there is only one 'subject' (i.e. dimension 1)
      erro = permute(Yb_sd(:,i,:),[3,1,2]);     % then use the within-bin stdev as error bar
    else% otherwise, multiple subjects, so use sem across subjects
      erro = sq(nanstd(demeaned_t_y(:,j,:)))/sqrt(size(Yb,1)); % std err at each time
    end
    if numel(erro)==1, erro=erro*ones(size(line)); end
    if SMOOTH, line=smoothn(line, SMOOTH); erro=smoothn(erro, SMOOTH); end
    args=varargin;                % which arguments to send to the plot command
    if any(strcmpi(args,'color')) % has color been provided as a parameter?
      ci = find(strcmpi(args,'color'));
      colj= args{ci+1};
    else                          % if not, cycle through list of colours
      % if there is more than one line, add in the color-order
      colj = cols(mod(j-1,size(cols,1))+1,:); 
      %if size(Xb,2)>1,           % if there's more than one condition,
        args=[args 'Color',colj]; % set colour for each condition separately
      %end
    end
    if FILL_AREA<3 % if not quantiled, plot the mean as a line with dots
      h2=plot(x, line, args{:}, 'Marker','.',plotargs{:}); h=[h h2];
      hold on
    end
    if ~FILL_AREA % show error bars as DOTTED LINES?
      h2=plot(x, line+erro, args{:}, 'LineStyle',':',plotargs{:}); h=[h h2];
      h2=plot(x, line-erro, args{:}, 'LineStyle',':',plotargs{:}); h=[h h2];
    elseif FILL_AREA<3 % show error bars as FILLED AREA?
      alpha = 0.5;
      bad = isnan(x) | isnan(line) | isnan(erro) | isinf(erro) | isinf(x) | isinf(line);
      h2 = fill(  [x(~bad); flipud(x(~bad))]', ...
                  [line(~bad)+erro(~bad); flipud(line(~bad)-erro(~bad))]', ...
                  colj, 'FaceAlpha', alpha,'linestyle','none');
      h=[h h2];
    else % if fill>=3, show QUANTILES? (see errorBarPlot code)
      errorBarPlot( sq(t_y(:,j,:)), 'xaxisvalues',nanmean( sq(t_x(:,j,:)) ) , ...
        'area',FILL_AREA,'smooth',SMOOTH,'plotargs',{'marker','.', plotargs{:}},'color',colj);
    end
    if PLOT_INDIVIDUALS
      plot( squeeze(t_x(:,j,:))',  squeeze(t_y(:,j,:))' ,'color',colj); 
    end
    hold on
  end
  if PLOT_TRAJECTORIES 
    plot( sq(nanmean(t_x)) , sq(nanmean(t_y)) ,'color',[.5 .5 .5] );
  end
  %%%%%%% STATISTICS on sliding bin data?
  if DO_STATS % DO STATS?  treats conditions as a continuous variable
    NS = size(Yb,1); NC = size(Yb,2);      % run mixed effects permutation test
    if nargout<5 || ~DO_STATS_PLOT, 
      ao={[],[],[],[]};        % did you request a figure handle? if so 
    else
      ao={[],[],[],[],[]};     % then request a significance bar from permutationOLS.
    end 
    % create a design matrix: column of 1s, then column of conditions (assume linear effect) 
    DES=[flat(ones(NS,NC))  zscore(flat(repmat(1:NC,[NS,1]))) ];
    if var(DES(:,2))>0  % assuming there are different conditions to compare:
      [ao{:}]=permutationOLS(reshape(t_y,NS*NC,[]),... % get list of traces over values of X
        DES, [0 1], ...            % statistics: contrast for linear effect of conditions
        repmat([1:NS]', NC,1)  );  % group by subject (i.e. permute conditions within subjects)
      p = ao{2}; t = ao{3};        % p-value and t-statistic
    else p=[]; t=[];    % only one condition - no stats possible
    end
  end
end

if exist('labels','var')
  if length(labels)>0,   xlabel(labels{1}); end
  if length(labels)>1,   ylabel(labels{2}); end
  if length(labels)>3,   legend(labels{4}); end
end


if ~oldhold % restore "graphics hold" state
  hold off;
end


% restructure the residuals to match the input format of X
if RESIDUALS && ~wascells
  resid = nancat([2,3],nancat([1,2],resid)); 
  resid = resid{1};
end


