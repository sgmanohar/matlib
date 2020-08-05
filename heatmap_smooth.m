function varargout=heatmap_smooth(X,varargin)
%   heatmap( XY )
%   heatmap( X, Y )
%   heatmap( X, ... )
% 
% convenience function to call 
%   contourf( hist3( X ) )
% plot heatmap of paired X,Y data. 
%
% params:   'windows',  [ ngridX, ngridY ]  or  { centresX(:), centresX(:) }
%       or  'smooth',   gaussian width
%           'kernel',   'gauss' or 'flat'
%           'plotargs',  extra arguments passed to contourf.
%           'contour',   if false, just plot as squares
%           'transform', apply transform to the probability density (e.g. log) 
% sgm 2014

if isvector(X) && nargin>1 && all(size(varargin{1})==size(X))
  X=[X varargin{1}]; 
  varargin(1)=[];
end
  
[KERNEL D SM plotargs, CONTOUR, TRANSFORM] =parsepvpairs({
  'kernel','windows','smooth','plotargs','contour', 'transform' },{
  'gauss',  [] ,  [] , {}, true, [] }, varargin{:});

if isnumeric(D) && ~isempty(D)
  if isscalar(D), D=[D D]; end
elseif ~iscell(D)
  D = [1 1] * floor(sqrt(length(X))); 
end
if isempty(SM), SM=floor(D(1)/10);end;

% KERNEL = 'gauss'; % smoothing kernel
% WINDOWS = 15;  % number of windows per line, for smoothing
% SM = floor(D/WINDOWS); 

if strcmpi(TRANSFORM,'log')
  TRANSFORM = @(x)eps+log(x+eps); 
end

[h , c]= hist3(X,D); 
if ~isempty(TRANSFORM)
  h=TRANSFORM(h);
end

if SM
  switch KERNEL
    case 'gauss'
      lin    = linspace(-2.2,2.2,SM(1))';
      kernel = exp(-bsxfun(@plus, lin.^2, lin'.^2));
      kernel = kernel./sum(kernel(:));
    case 'flat'
      kernel = ones(SM); 
  end
  h = conv2(h,kernel,'same');
end
if CONTOUR
  if 0 % use axis of integer value bins
    contourf(h,plotargs{:},'edgecolor','none');
    set(gca,'xticklabel', ...
      arrayfun(@(x)sprintf('%0.2g',x), c{1}(get(gca,'xtick')),'uni',0), ...
      'yticklabel',...
      arrayfun(@(x)sprintf('%0.2g',x), c{2}(get(gca,'ytick')),'uni',0) ...
      );
  else % use actual values as axes
    contourf(c{1}, c{2}, h', plotargs{:}, 'edgecolor','none');
  end

else
  ht=imagesc(c{1},c{2},h');                           % show heatmap
  set(gca,'ydir','normal');                 % up is up!
  xt=get(gca,'xtick');
  yt=get(gca,'ytick');                      % set axes scales
  if all( xt==floor(xt) & xt>0 ) && all( yt==floor(yt) & yt>0 )
    makelabels = @(y) arrayfun( @(x)sprintf('%0.2g',x), y, 'uniformoutput',0);
    set(gca,'xticklabel',makelabels(c{1}(xt)), 'yticklabel',makelabels(c{2}(yt)) );
  end

end
if nargout>0,   varargout{1}=h;  end
if nargout>1,   varargout{2}=c;  end

