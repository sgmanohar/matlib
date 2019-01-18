function [b bint r rint stats handles]=scatterRegress(X,Y, varargin)
% [b bint r rint stats handles]=scatterRegress(X,Y, varargin)
%
% show a scatter of the two vectors X and Y
% perform a linear regression and return what would be expected from
% REGRESS.  The first regressor is X, the second is a constant term.
% the variable args get sent to SCATTER.
% 
% STATS is: 
% [ R-square, F statistic, p value for the full model, 
%   and estimate of the error variance. ]
%
% 2015: you can now pass X Y and Z to do a 3D scatter.
% in this case, Z is regressed against X, Y, X*Y, and constant, in that order.
% options: 
%   'pearson' : 1 = pearson correlation (default); 0 = spearman
%   'plot_ci' : plot regression confidence intervals (parabolic curves)
%   'fill'    : fill in the confidence interval as a shaded area, if plot_ci
%               is true
%   'plane'   : for 3D data, fit a plane to the data
%   'plot_2d' : only do 2D plot, even if 3D data supplied. Extra parameters get
%               passed to 'scatter', e.g. as Size and Color of points.
%   'DemeanXByColumn', 'DemeanYByColumn' : if the data is in multiple columns,
%               then align the means of each of the columns before plotting. 
% 
%   'Density' : For 2D data (X,Y), plot a 3D density surface over the
%               data points, using kernel density (gaussian convolution).
%               The specified value after 'density' determines the width
%               of the kernel.
%   'ShowZero': If the axes include zero, draw a 'zero' line intersecting 
%               the axis.
% 
PLOT_CI  = 1; % draw regression confidence intervals? 1 = if signif, 2= always
FILL     = 1; % fill in the confidence interval
FIT_SURF = 1; % use gridfit (for 3D data) to fit a z-surface
PLANE    = 1; % for 3D, draw the edges of the regression plane on the walls
PEARSON  = 1; % pearson or spearman correlation?
PLOT_2D  = 0; % plot only in 2D (extra dimensions go to 'scatter' function)
DEMEAN_X_BY_COLUMN = false; % align each column's X-values to the same global mean?
DEMEAN_Y_BY_COLUMN = true; % align each column's Y-values to the same mean?
DENSITY  = 0;
SHOW_ZERO = 1; % plot dotted lines at 'zero-zero', if within the axis
ALPHA    = 0.05; % significance threshold for plotting confidence intervals

i=find(strcmpi(varargin,'pearson'));
if i, PEARSON=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'plot_ci'));
if i, PLOT_CI=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'fill'));
if i, FILL=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'plane'));
if i, PLANE=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'plot_2D'));
if i, PLOT_2D=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'DemeanXByColumn'));
if i, DEMEAN_X_BY_COLUMN=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'DemeanYByColumn'));
if i, DEMEAN_Y_BY_COLUMN=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'Density'));
if i, DENSITY=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'ShowZero'));
if i, SHOW_ZERO=varargin{i+1}; varargin([i i+1])=[]; end

if ~exist('Y','var') || ~isnumeric(Y)
  if size(X,2)==2 && ndims(X==2)
    Y=X(:,2);X=X(:,1);
  elseif size(X,2)==3 && ndims(X==2)
    Y=X(:,2);X=X(:,1);varargin=[{X(:,3)} varargin];
  end
end

ohold=ishold();

if ~isempty(varargin) && isvector(varargin{1}) && length(varargin{1})==length(X) && ~PLOT_2D
  % provided 3D data!
  Z=varargin{1}; varargin(1)=[];
  % but we can deal with this...
  handle=scatter3(X,Y,Z,varargin{:});
  [b bint r rint stats] = regress( Z(:), [X(:), Y(:), X(:).*Y(:), ones(size(Y(:)))]);
  hold on
  xl=xlim; yl=ylim; zl=zlim;
  if PLANE 
    for i=1:4 % project regression plane onto the walls.
      switch i,
        case 1, xp = xl(1)*[1 1]; yp = yl;
        case 2, xp = xl(2)*[1 1]; yp = yl;
        case 3, xp = xl;          yp = yl(1)*[1 1];
        case 4, xp = xl;          yp = yl(2)*[1 1];
      end
      handles=plot3(xp,yp, [ xp; yp; xp.*yp; 1 1 ]'*b);
    end
  end
  if FIT_SURF % fit a smooth z-surface to the data
    [zsurf xsurf ysurf]=gridfit(X(:),Y(:),Z(:), ...
      linspace(min(X(:)),max(X(:)),50), linspace(min(Y(:)),max(Y(:)),50), 'smoothness',0.7); 
    surf(xsurf,ysurf,zsurf, 'edgecolor','none'); camlight; lighting phong; 
  end
  % annotate with p-values
  signif  = prod(bint,2)>0; 
  symbols = {'X* ','Y* ','XY*'}; 
  anno    = [symbols{signif(1:3)}]; 
  if any(anno)
    text(mean(xlim),max(ylim), anno);
  end
  axis vis3d
else % 2D data
  if isrow(Y), Y=Y'; end
  if isrow(X), X=X'; end
  if size(X,2)>1 && all(size(X)==size(Y)), % X is a matrix? loop several times
    if DEMEAN_X_BY_COLUMN % subtract column means, and add global mean
      mx=nanmean(X); X=bsxfun(@minus, X,mx) + mean(mx); 
    end
    if DEMEAN_Y_BY_COLUMN % subtract column means, and add global mean
      my=nanmean(Y); Y=bsxfun(@minus, Y,my) + mean(my); 
    end
    handles=[];
    for i=1:size(X,2)
      [b(i,:) bint(i,:,:) r(i,:) rint(i,:,:) stats(i,:) h] = scatterRegress(X(:,i),Y(:,i),...
        varargin{:},'plot_ci',0,'pearson',PEARSON);  % pass on parameters, but don't do CI
      handles=[handles h]; % keep figure handles
      set(h(1), 'CData', colourMap(i,size(X,2)) );  % set colours to match.
      set(h(2), 'Color', colourMap(i,size(X,2)) ) ;
      hold on
    end 
  else % X is a vector, Y might be a matrix though
    Yall=Y;
    if size(Yall,2)>3, PLOT_CI=false; end
    handles = [];
    for i=1:size(Yall,2) % do each column of Y separately
      Y=Yall(:,i); % select column
      h=scatter(X,Y, varargin{:});
      handles=[handles h];
      [b bint r rint stats] = regress(Y(:), [X(:), ones(size(Y(:)))]);
      if PEARSON
        anno=sprintf('r^2=%g\np=%g',stats(1),stats(3));
        p=stats(3);
      else % SPEARMAN
        bad = isnan(X) | isnan(Y);
        [rho p] = corr( X(~bad),Y(~bad)  , 'type','Spearman');
        anno=sprintf('\\rho=%g\np=%g',rho,p);
      end
      isSignif = p<ALPHA;
      hold on;
      h = plot(xlim(),b(1)*xlim()+b(2));
      handles=[handles h];
      % annotation('textbox',[0.5,0.5, 0.1 0.1],'String',sprintf('r^2=%g\np=%g',stats(1),stats(3)), 'LineStyle','none')
      text(mean(xlim),max(ylim),anno);
      if PLOT_CI && exist('regression_line_ci') && (isSignif || PLOT_CI>1)
        [top_int,bot_int,x_int]=regression_line_ci( 0.05, flipud(b), X(:),Y(:) , 100, min(xlim), max(xlim));
        colororder = get(gca,'colororder');
        if FILL
          fill([x_int fliplr(x_int)], [top_int fliplr(bot_int)], colororder(2,:), ...
            'FaceAlpha', 0.5,'linestyle','none');
        else
          plot(x_int,top_int,':','color',colororder(2,:));
          plot(x_int,bot_int,':','color',colororder(2,:));
        end
      end
      if DENSITY
        z = [X(:) Y(:)]; z(isnan(X) | isnan(Y),:)=[];
        gridx=linspace(nanmin(X(:)),nanmax(X(:)),40);
        gridy=linspace(nanmin(Y(:)),nanmax(Y(:)),40);
        if islogical(DENSITY), DENSITY=[]; end % if not numeric, don't pass it to ksdensity2d
        d = ksdensity2d(z,gridx,gridy ,DENSITY );
        surf(gridx,gridy,d, 'edgecolor','none');    camlight; lighting phong;
        hold on; plot3(z(:,1),z(:,2),zeros(size(z,1),1),'bo');
      end
      b_all{i}=b; bint_all{i}=bint; r_all{i}=r; rint_all{i}=rint; stats_all{i}=stats;
    end
  end
end
if SHOW_ZERO && prod(xlim)<0
  plot([0 0], ylim,':');
end
if SHOW_ZERO && prod(ylim)<0
  plot(xlim,[0 0],':');
end
if(~ohold) hold off; end;
