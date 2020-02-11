function [acc,mrt, p,t, h]=conditionalAccuracyPlot(RT,Corr, Nbins, styles, varargin)
% [acc, mrt, p, t, h] = conditionalAccuracyPlot(RT,Corr, nbins, styles, ...)
%      RT{ SUBJECT, CONDITION } ( TRIAL )
% or   RT( TRIAL, SUBJECT )            -- for one condition
% or   RT( TRIAL, SUBJECT, CONDITION ) -- for multiple conditions
%
% columns are conditions, plotted as separate lines.
% styles is a cell with a string for each condition
%
% note that, since the function uses nanmean() to calculate the accuracy,
% instead of plotting a percentage correct, you could
% remove the transform and use this to plot the mean of any value when 
% binned by RT.
%
% extra params e.g. 'smooth', 'alpha' are passed to errorBarPlot. except:
%  'meanfunc' = function handle to use instead of "nanmean". 
%
% returns: ACC ( SUBJECT, CONDITION, QUANTILE ) = mean of Corr for each bin
%          MRT ( SUBJECT, CONDITION, QUANTILE ) = median RT for each bin
%          p ( QUANTILE ) = probability that conditions are significantly 
%             different at a quantile, using permutation test
%             (treats conditions as a single continuous variable)
%          T ( QUANTILE )  = t statistic for linear effect of condition
%          H = handle of lines as returned by errorbarplot
%
% 'xtransform', 'ytransform' - set transformations for the plot axes. Can
%    be a @function, or one of 'identity' (default), 'probit', 'recip',
%    'log', 'asin', 'gauss' (inverse gaussian from 0 to 1), or
%    'gauss2' (inverse gaussian from 0.5 to 1). Note that the output values
%    acc and mrt will not be transformed.
% 'StandardBinning' - if true, then use non-overlapping bins.
% 
% sanjay g manohar

if(~exist('Nbins','var') || isempty(Nbins)) Nbins=5; end;

if(~exist('styles','var') || isempty(styles)) 
  colstr='rgbcmyw';
  for(i=1:size(RT,2))
    styles{i}=colstr(mod(i-1,length(colstr))+1);
  end
end
if ~iscell(RT)
  if ndims(RT)==2 % treat as RT( TRIAL, SUBJECT )
    RT   = mat2cell(RT,   size(RT,  1),   ones(1,size(RT,2))   )';  % break into cells
    Corr = mat2cell(Corr, size(Corr,1),   ones(1,size(Corr,2)) )';
  elseif ndims(RT)==3 % treat as RT ( TRIAL, SUBJECT, CONDITION ): convert to cell { SUBJECT, CONDITION } (TRIAL)
    RT   = sq(mat2cell(RT,   size(RT,  1),   ones(1,size(RT,  2)),   ones(size(RT,  3),1) )); 
    Corr = sq(mat2cell(Corr, size(Corr,1),   ones(1,size(Corr,2)),   ones(size(Corr,3),1) ));
  end
end

% use nonoverlapping bins? 
% advantage is the values are independent, and thus stats are simpler.
% if not, use overlapping bins, which must be compared with permutation test 
STANDARD_BINNING=0; 
Nq=100; % number of quantiles (used when STANDARD_BINNING=0)
binwidth = 1/Nbins;
SMOOTH   = 15; % used for manual area plot at bottom. 
MEANFUNC = @nanmean; % can use nanstd here if we want variability

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

i=find(strcmpi(varargin, 'StandardBinning')); if i, 
  STANDARD_BINNING = varargin{i+1}; varargin([i i+1])=[]; 
end

i=find(strcmpi(varargin, 'Conditions')); if i, 
  C = varargin{i+1}; varargin([i i+1])=[]; 
  if ~isempty(C)
    if size(RT,2)>1, error('cannot specify conditions since RT already contains conditions'); end
    if isnumeric(C) & ndims(C)==2, C=mat2cell(C, size(C,1), ones(1,size(C,2)) )'; end
    uc=unique(nancat(2,C{:})); uc(isnan(uc))=[]; 
    for i=1:size(RT,1) % subject
      for j=1:length(uc) % condition
        rt{i,j}  =RT{i}(  C{i}==uc(j) );
        corr{i,j}=Corr{i}(C{i}==uc(j) );
      end
    end
    RT=rt; Corr=corr;
  end
end


if STANDARD_BINNING, N=Nbins, 
else N=floor(Nq*(1-binwidth));  % number of points on a CAF curve
end
mrt=nan(size(RT,1), size(RT,2), N); % create blank output variable (needed for parallel processing)
%%%%%%%% change the outer 'for' to 'parfor' for faster execution.
for(i=1:size(RT,1)) % for each subject
  mrt_i=nan(size(RT,2),N);
  for(j=1:size(RT,2)) % for each condition
    if STANDARD_BINNING % DIVIDE INTO BINS
      top=-inf;
      for q=1:Nbins   % bin the RTs into Nbins
        bottom=top;
        if(q<Nbins)   top=quantile(RT{i,j},q/Nbins);
        else          top=inf;
        end
        f=RT{i,j}>bottom & RT{i,j}<top;
        acc(i,j,q)=MEANFUNC( Corr{i,j}(f) );    % mean of Corr
        mrt_i(j,q)=quantile(RT{i,j},(q-0.5)/Nbins); % centre (median) RT of the bin
        % alternative: use mean of top and bottom of bin?
      end
    else % SLIDING WINDOW BINS
      quantiles = quantile(RT{i,j}, linspace(0,1,Nq)); % create quantiles
      % go through bins, but there aren't 100 bins because of the width!
      for q=1:floor(Nq*(1-binwidth));   % e.g. 1:80 for 5 bins
        % get start, mid and end of the bin
        range=quantiles( q + [0,floor(Nq*binwidth/2), floor(Nq*binwidth)] );
        f = RT{i,j}>range(1) & RT{i,j}<range(3); % filter data within range
        % ACCURACY ( SUBJECT, CONDITION, QUANTILE )
        acc(i,j,q)=MEANFUNC( Corr{i,j}(f) );      % mean of Corr 
        mrt_i(j,q)=range(2); % centre (median) RT of the bin        
      end
      % %% 8/7/14 added to prevent nonunique bin centres
      mrt_i(j,:)=removeIdenticals(flat(sq(mrt_i(j,:))));
    end
  end
  mrt(i,:,:) = mrt_i;
end
if exist('TRANSFORM') acc=TRANSFORM(acc); end % apply transform to mean

doErrorBars=1;

oldhold=ishold();
t_mrt = XTRANSFORM( mrt );
t_acc = YTRANSFORM( acc );
h=[]; 

if STANDARD_BINNING
  for(j=1:size(RT,2))
    if doErrorBars
      errorbarxy(  YTRANSFORM(  sq(mean(t_mrt(:,j,:),1))  ), ...
        sq(mean(t_acc(:,j,:),1)), ...
        sq(std(t_mrt(:,j,:),[],1))/sqrt(size(mrt,1)), ...
        sq(std(t_acc(:,j,:),[],1))/sqrt(size(acc,1)) ...
        ,[],[],styles{j}, .2*[1 1 1]);
    else
      h2=plot(sq(mean(t_mrt(:,j,:),1)), ...
        sq(mean(t_acc(:,j,:),1)) ...
        , styles{j});
      h=[h h2];
    end
    hold on
  end
else % SLIDING BIN = area plot
  % ACC ( SUBJECT, CONDITION, TIME )
  if 0 % Use errorbarplot to do areas?
    h2=errorBarPlot( permute(t_acc, [1 3 2]), ...
       'area', 1, ...
       'alpha',0.5, ...
       'xaxisvalues', sq(mean(mean(t_mrt(:,:,:),1),2)) ... 
       ,varargin{:}); % problem is, it doesn't get the x axis right (should be different for each condition)
    h=[h h2];
  else % so, do it manually?
    cols=get(gca,'ColorOrder');
    for j=1:size(mrt,2) % for each condition
      x = sq(nanmean(t_mrt(:,j,:),1)); line = sq(nanmean(t_acc(:,j,:),1)); 
      erro = sq(nanstd(t_acc(:,j,:)))/sqrt(size(acc,1)); % std err at each time
      if numel(erro)==1, erro=erro*ones(size(line)); end
      if SMOOTH, line=smooth_old(line, SMOOTH,'gauss'); erro=smooth(erro,SMOOTH,'sgolay'); end
      args=varargin;
      % if there is more than one line, add in the color-order
      colj = cols(mod(j-1,size(cols,1))+1,:);
      if size(mrt,2)>1, args=[args 'Color',colj]; end 
      h2=plot(x, line, args{:}, 'Marker','.'); h=[h h2]; 
      hold on
      if 0 % DOTTED LINES?
        h2=plot(x, line+error, args{:}, 'LineStyle',':'); h=[h h2];
        h2=plot(x, line-error, args{:}, 'LineStyle',':'); h=[h h2];
      else % FILLED AREA?
        alpha = 0.5;
        h2 = fill([x; flipud(x)]', [line+erro; flipud(line-erro)]', colj, 'FaceAlpha', alpha,'linestyle','none');
        h=[h h2];
      end
    end
  end
  if nargout>2 % DO STATS?  treats conditions as a continuous variable
    NS = size(acc,1); NC = size(acc,2); % run mixed effects permutation test
    if nargout<5, ao={[],[],[],[]};  % did you request a figure handle? if so 
    else          ao={[],[],[],[],[]}; end % then request a significance bar from permutationOLS.
    X=[flat(ones(NS,NC))  zscore(flat(repmat(1:NC,[NS,1]))) ];
    if var(X(:,2))>0
      [ao{:}]=permutationOLS(reshape(acc,NS*NC,[]),...
        X, [0 1], ...
        repmat([1:NS]', [1,NC])  ); 
      p = ao{2}; t = ao{3};
    else p=[]; t=[];
    end
  end
end

if ~oldhold
  hold off;
end


