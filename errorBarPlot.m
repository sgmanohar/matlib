function varargout = errorBarPlot(X, varargin)
% h = errorBarPlot(X, varargin)
%   X is an array
%     Plot with within-subject standard error bars. 
%     Means are taken over dimension 1.
%   X is a cell array { Table, 'yvar',  'grouping' , 'xvar' , 'linesvar'  }
%     grouping and linesvar are optional, eg: { T, Y, X } or { T, Y, G, X } 
%     If grouping is present, means are taken for each group, and the 
%       error bars will reflect the groupings (e.g. subjects). If not, error
%       bars are across all rows that match a level of xvar.
%     If linesvar is present, then multiple lines are plotted, for each
%       level of this variable
% options:
%           'area': 0 = error bars only
%                1 = shaded transparent area
%                2 = use dotted lines for min and max of std error. 
%                3 or greater: show Quantiles. 3=[.33 to .66], 4=[.25 to .75],
%                  5=[.2 to .8, and .4 to .6], etc. Note: median used, meanfun ignored. 
%               -1 = 'Violin' plot, one-sided; -2 = symmetrical area.
%           'alpha': after area - use an alpha for call to 'fill' when
%                plotting . Default 0.5
%           'xaxisvalues': specify the x-values of the graph. Otherwise
%                just uses 1,2,3...
%           'maineffect' dim 
%                dim should be 2 or 3 to indicate a main effect of
%                dimension 2 or 3 of the data.
%                subtract the mean for that dimension from all values
%                before averaging and taking std error. then add on the
%                mean before plotting.
%           'withinSubjectError' - subtract subjects' intercepts before
%                calculating error bars. Default 0.
%           'standardError': 1 = use SEM for error bars (default).
%                2 = use standard deviation
%                between 0 and 0.5: use percentile above and below median
%                otherwise : calculate bootstrapped confidence intervals at p=0.05.
%                    In that case, also specify 'nBoot' = num permutations.
%           'width':  adjust the width of the errorbars. needs 'errorbarT.m'
%               function (download from Matlab Central)
%           'type': 'line' (default) or 'bar'. Note, 'area' only works with
%                line. 'bar' needs 'barwitherr.m', download from Matlab Central.
%           'meanfun' : function to calculate mean. Default: @nanmean
%                if using this, and you want error bars to reflect the
%                deviation of this value, use 'standardError' false.
%           'doStats' : do repeated measures anova? default true for small datasets. 
%           'plotargs':  send the following items to 'plot'. Arguments
%                should be in a cell array.
%           'statargs':  send the following items to 'rmanova' or the
%                statistical function.
% Sanjay Manohar 2014

% 1. Decide upon parameters
DEFAULT_WITHIN_SUBJECT_ERROR = 1;
if exist('nanmean','file'),   fmean = @nanmean; % mean function
else                          fmean = @mean; end
Areas.VIOLIN = -1;
try
  % read parameters
  [medim, AREA, alpha, xaxisvalues, SMOOTH, plotargs,...
    WITHINSUBJECTERROR, STANDARD_ERROR, NBOOT, WIDTH, COLOR, ...
    TYPE, fmean, PLOT_INDIVIDUALS, DO_STATS, STAT_ARGS, LABELS] = parsepvpairs( ...
  {'MainEffect', 'Area', 'Alpha', 'xaxisvalues', 'smooth', 'plotargs', ...
    'withinSubjectError', 'standardError', 'nBoot' , 'width','color', ...
    'type' , 'meanfun','plotIndividuals','doStats','statargs', 'labels'}, ...
    {[], false, 0.5,     [],         1,       {} , ...
     DEFAULT_WITHIN_SUBJECT_ERROR, 1,     5000,   [],     [] ...
     'line', @nanmean ,false, [], {}, {}}, ... 
    varargin{:});
  varargin={};
catch me
  medim=[]; AREA=false; alpha=0.5; xaxisvalues=[]; SMOOTH=1; 
  WITHINSUBJECTERROR=DEFAULT_WITHIN_SUBJECT_ERROR;
  PLOT_INDIVIDUALS = false; 
  STANDARD_ERROR=1; WIDTH=[]; COLOR=[]; TYPE='line';
  DO_STATS= []; STAT_ARGS = {}; LABELS = {};
  warning('all varargins to errorbarplot will now be passed to Plot. Use plotArgs next time, please.');
  plotargs=varargin;
end
DOTTED_AREA = AREA==2; % use dotted line above and below, instead of area plot

if ~isempty(COLOR), plotargs=[plotargs {'color',COLOR}]; end

%%% HANDLE TABLES AS INPUT
if iscell(X) && isa( X{1}, 'table' ) % convert table to array
  % handle three possible formats of the variables:
  if numel(X)==3 % { T, Y, X }
    tmp = [ X{1}.(X{3}) X{1}.(X{2}) ]; % [X, Y]
    tmp = pivot(tmp)';
    LABELS = X([3,2]); % X and Y label
  elseif numel(X)==4 % { T, Y, G, X }
    tmp = [ double(X{1}.(X{3})) double(X{1}.(X{4})) double(X{1}.(X{2})) ]; % [G, X, Y]
    tmp = pivot(tmp);
    tmp = sq(nanmean(tmp, ndims(tmp))); % take means within groups
    % tmp = permute(tmp, [2,1]); % tmp( group, xvar )
    LABELS = X([4,2]); % X and Y label
  elseif numel(X)==5 % { T, Y, G, X, Z } 
    tmp = [ X{1}.(X{3}) X{1}.(X{4}) X{1}.(X{5}) X{1}.(X{2}) ]; % [G, X, Z,  Y]
    tmp = pivot(tmp); % tmp( group, xvar, zvar, samples_in_group )
    tmp = sq(nanmean(tmp, ndims(tmp))); % take means within groups
    LABELS = X([4,2,5]); % X, Y and line labels
    %tmp = permute(tmp, [
  else error('please provide names of two or three table columns.');
  end
  X=tmp;
end

%%% HANDLE LME OBJECT AS INPUT
if isa(X,'LinearMixedModel') 
  C=X.Coefficients;                    % get coefficients
  if strcmpi(C.Name(1),'(Intercept)'), % ignore intercept
    C = C(2:end,:);
  end
  barwitherr(C.SE,  C.Estimate);       % plot bar with error
  set(gca,'XTick',1:length(C.SE));     % all bars labelled
  set(gca,'XTickLabel',C.Name);        % label x axis
  tmp=gca; if isprop(tmp,'XTickLabelRotation')
    tmp.XTickLabelRotation = 45;
  end
  % y location to draw asterisks: 95% of height
  yl = ylim; ht=yl(2); dh = 0.05*(yl(2)-yl(1)); 
  for i=1:size(C,1)    % for each coefficient
    if C.pValue(i) < 0.05 % is it significant?
      text(i, ht-dh , '*'); 
      text(i, ht-2*dh+(dh/2)*(-1)^i, sprintf('p=%0.2f',C.pValue(i)));
    end
  end
  return;
end

% check data looks 'paired'
if size(X,1)>5 && any(any(all(isnan(X(end-5:end,:,:)),1)))
  warning('data looks unpaired. error bars will be between-subjects')
  WITHINSUBJECTERROR = false;
end

sx=size(X); 

% 2. Calculate main effects 
addon=0;
if ~isempty(medim) % main effect analysis? 
    % remove means across all dimensions except the specified one 
    % before doing standard deviation
    if(medim==3)
        addon=sq(repmat(fmean(fmean(X,1),2),[1,size(X,2),1]));
        X=X-repmat(fmean(X,3),[1,1,size(X,3)]);
    elseif medim==2
        addon=sq(repmat(fmean(fmean(X,1),3),[1,1,size(X,3)]));
        X=X-repmat(fmean(X,2),[1,size(X,2),1]);
    else error('Please could you use MainEffect on dimension 2 or 3');
    end
end

% 3. Calculate error bars
if WITHINSUBJECTERROR  % within subject error bars?
  % subtract between-subject means first
  M=fmean(fmean(fmean(fmean(X,1),2),3),4);
  m2 = fmean(X,2); if(size(m2,3)>1) m2=fmean(fmean(m2,3),4); end
  Xm=X-repmat(m2, [1, sx(2:end)] ); % subtract subject mean
  Xm=Xm+M;
  if all( sum(sum(~isnan(X),2),3) < 2 ), % are they all "one condition per subject"?
    Xm=X; 
    warning('one condition per subject: using between subjects error');
  end
else % Otherwise use raw values to compute standard error
  Xm=X;
end
NS = sum(~isnan(X),1); % num subjects

% use standard error? or bootstrapped confidence intervals?
if STANDARD_ERROR==1 % standard error calculation 
  yerror  = sq(bsxfun(@rdivide, sqrt(nanvar(Xm,[],1)), sqrt(NS)));
  yerrorm = yerror; % 'minus' errors same as 'plus' errors
elseif STANDARD_ERROR ==2 % use standard deviation?
  yerror  = sq(sqrt(nanvar(Xm,[],1)));
  yerrorm = yerror;
elseif STANDARD_ERROR < 0.5 && STANDARD_ERROR > 0 % treat fraction as percentile for bar
  yerror  =   quantile(Xm, 0.5+STANDARD_ERROR) - fmean(Xm);
  yerrorm = -(quantile(Xm, 0.5-STANDARD_ERROR) - fmean(Xm));
else              % bootstrap confidence interval at 5%
  yerror  = bootci(NBOOT, {fmean, Xm}, 'type', 'bca');
  yerrorm = -squeeze(yerror(1,:,:,:))+squeeze(fmean(Xm)); % minus errors
  yerror  = squeeze(yerror(2,:,:,:))-squeeze(fmean(Xm)); % plus errors
end

% 4. Draw the graph
ho=ishold(); % preserve hold-on status
h=[]; % store axes/graph handles when plotted.
if ~AREA % normal points/lines
  Y=sq(fmean(X,1))+addon; % get the mean across subjects
  if isrow(Y) Y=Y'; end;  % make sure it is a column
  if isempty(xaxisvalues) % X axis values not supplied?
    xaxisvalues = repmat([1:size(Y,1)]', [1,size(Y,2),size(Y,3)]);
  else % check the x-axis values shape is correct
    if (isrow(xaxisvalues) && iscolumn(Y)) || (iscolumn(Y)) && isrow(xaxisvalues)
      xaxisvalues=xaxisvalues';
    end;
    xaxisvalues=bsxfun(@(a,b)a,xaxisvalues, Y); % shape must match Y
  end

  switch TYPE
    case 'line',  
      % plot line with errorbars
      h=errorbar(xaxisvalues, Y,yerrorm,yerror, plotargs{:});
      if PLOT_INDIVIDUALS  % individuals as dotted lines
        hold on;
        plot(xaxisvalues,X',':',plotargs{:}); 
      end
    case 'bar', 
      % plot bar with errorbars
      h=barwitherr(cat(3, yerror, -yerrorm), xaxisvalues, Y, plotargs{:});
      if PLOT_INDIVIDUALS  % individuals as circles
        hold on;
        plot(xaxisvalues, X,'o',plotargs{:});  
      end
  end
  if WIDTH, errorbarT(h,WIDTH); end % adjust witdth of bars
else % AREA plot
  if DOTTED_AREA % use dotted lines above and below, instead of area?
    h1=plot(sq(fmean(X,1))+addon, plotargs{:}); % mean line
    hold on
    h2=plot(sq(fmean(X,1))+addon+yerror, plotargs{:},'linestyle',':');
    h3=plot(sq(fmean(X,1))+addon-yerrorm, plotargs{:},'linestyle',':');
    h=[h1 h2 h3];
  else % shade in area between errors
    if ~any(strcmpi(plotargs,'color')) % color not explicitly given?
      col=get(gca,'ColorOrder'); % use color order
    else % consume 'color' parameter
      ci=find(strcmp(plotargs,'color')); 
      col = plotargs{ci+1}; plotargs([ci ci+1])=[]; 
    end
    for i=1:size(X,3) % for each line
      coli = col(1+mod(i-1,size(col,1)),:); % color of line
      if ~isempty(xaxisvalues)
        % transpose to make it a column if needed
        if isvector(xaxisvalues) && size(xaxisvalues,1)==1, xaxisvalues=xaxisvalues'; end % make column
        % need to subselect columns of xaxisvalues?
        if size(xaxisvalues,2)>1, xav = xaxisvalues(:,i); else xav=xaxisvalues; end
      else % no x axis values specified; use 1 to N
        xav = [1:size(X,2)]'; % average x values
      end
      if AREA==1 % error region shading
        xm = fmean(Xm(:,:,i));   % calculate mean
        % fill in area as a polygon
        xx = [xav' fliplr(xav')]; % x coords of polygon for plotting
        yy = [xm+yerror(:,i)' fliplr(xm-yerrorm(:,i)')]; % y coords for plotting
        bad = isnan(yy);
        xx(bad)=[]; yy(bad)=[]; % simply remove nans from the polygon
        if isempty(yy) warning('line %g omitted because no data',i); continue; end
        h{i,2}=fill(xx ,yy, coli, 'linestyle','none','FaceAlpha',alpha);
        hold on;
        % plot mean line on top of area
        h{i,1}=plot( xav, xm, 'color', coli,'linewidth',1.5,plotargs{:} );
      elseif AREA>2 % AREA > 2 ==> QUANTILES
        % number of quantiles to calculate
        nquant = AREA; % AREA=3 gives 0.33,.67, AREA=4 gives .25,.75 etc.
        quantiles = [1:nquant-1]/nquant; % which quatiles to plot
        q = quantile(X(:,:,i), quantiles); % should give a [ nquant x size(X,2) ] array
        if SMOOTH % make the curves artificially smoother?
          q=smoothn(q',SMOOTH)'; % use a window convolution
        end
        for j=1:floor((nquant-1)/2); % for each pair of quantile lines, 
          % fill an area between them, using an alpha depending on the
          % number of quantiles to make the transparency look graded
          h(j)=fill( [xav', fliplr(xav')] ,  [ q(j,:) fliplr(q(length(quantiles)+1-j,:)) ],  coli, 'linestyle','none','FaceAlpha',2*alpha/nquant );
          hold on;
        end
        % plot the median line
        h(end+1)=plot( xav, smoothn(sq(nanmedian(X(:,:,i))),SMOOTH) , 'color', coli , 'linewidth',1.5, plotargs{:})
        hold off;
      elseif AREA<0 % violin plot?
        dxav = mean(diff(xav)); % difference in x axis values = scale of pdf
        xm = fmean(Xm(:,:,i));   % calculate mean
        yav=linspace(-2.5*yerror, 2.5*yerror, 30);
        for k=1:size(xav,1) % for each x-point, 
          %xx=xav(k)+dxav*normpdf(yav(k,:)',0,yerror(k));
          %yy=yav(k,:)'+xm(k);
          [xx,yy]=ksdensity(Xm(:,k,i)); % estimate the PDF using kernel smoothed density
          dyav = yy(end)-yy(1); % difference in y values over support of distribution
          % calculate the horizontal position & scaling. account for the
          % density being dependent on the y-scaling. area = 1/10th of
          % inter-xval rectangle of support.
          xx = dxav * xx * dyav / 10 + xav(k);  
          if AREA==-1   % one-side: flat left edge
            yy=[yy';yy(1)]; xx=[xx';xx(1)];
          elseif AREA==-2 % two-sides: mirror like a violin
            yy=[yy'; flipud(yy')]; xx=[xx';xav(k)-flipud(xx'-xav(k))];
          end % now draw the fill:
          h(i,k) =  fill(xx,yy, ...
                         coli, 'linestyle','none','FaceAlpha',alpha ); 
          hold on
        end
        h(i,1)=plot( xav, xm, 'color', coli,'linewidth',1.5,plotargs{:} );
      end
      hold on
    end % next line
  end % dotted vs fill 
end % errorbars/area
if(~ho) hold off; end % restore hold-on status

if length(LABELS)>0
  xlabel(LABELS{1});
end
if length(LABELS)>1
  ylabel(LABELS{2});
end
if length(LABELS)>2
  legend(LABELS{3});
end

if isempty(DO_STATS)  % do stats was not specified. 
  if all(sx(2:end)<5),   DO_STATS = true;  % do stats only if few datapoints
  else
    if AREA, DO_STATS = true; 
    else     DO_STATS = false;
    end
  end
  if any(sx==1),         DO_STATS = false; end
end    
if DO_STATS
  try
    if ~AREA % non-area? use rm-ANOVA for stats
      stat = rmanova(X, STAT_ARGS{:});   % use simple repeated measures anova
      txt='';
      for row=2:size(stat.pValue) % plot each p value
        txt_r  = sprintf('F(%g,%g)=%0.3g,p=%0.3g', stat.DF1(row), stat.DF2(row), ...
          stat.FStat(row), stat.pValue(row) );
        if stat.pValue(row)<0.05, txt_r=[txt_r ' *']; end
        if isempty(txt),  txt = txt_r;
        else              txt = [txt char(10) txt_r];
        end
      end
      text(0.5,0.9, txt, 'Units','normalized', ...
        'HorizontalAlignment','center', 'VerticalAlignment','Top', 'editing','on');
    else % area? use permutation OLS
      NS = size(X,1); Nx = size(X,2); NC = size(X,3);
      if NC>2, warning('testing linear effect of condition'); end
      [stat.htest, stat.p_vals, stat.t_stat, stat.t_thresh, hplot] = permutationOLS(...
        reshape( permute(X,[1,3,2]),NS*NC,[]),... put x-values on dim 2, and subj x cond on dim 1
        [ ones(NS*NC,1), kron(zscore(1:NC),ones(1,NS))' ], [0 1], ...
        repmat([1:NS]', NC,1)  );
      h={h hplot};
    end
  catch mx
    fprintf('unable to do statistics:\n');
    disp(mx);
  end
end
if nargout > 0       % do we need to return the figure object handles?
  varargout{1} = h;
end
if nargout > 1
  varargout{2} = stat;
end
