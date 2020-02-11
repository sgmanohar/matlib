function varargout = heatmap(x,y, SIZE, GAUSS)
% OLD  - use heatmap_smooth
% handle = heatmap(X,Y)     
%             where X and Y are vectors
%          heatmap(X,Y,SIZE)       
%             where SIZE  = [a,b] number of bins
%          heatmap(X,Y,SIZE, GAUSS)
%             where GAUSS = number from 1 to a indicating radius of 
%             gaussian kernel, in bins
%          heatmap( X )  or
%          heatmap( X, SIZE [, GAUSS ] )
%             where X is N x 2 matrix of X and Y values
% 
% plot a heatmap scatter plot of X and Y coordinates.
% uses HIST3 to produce the 2D-histogram of bin counts
% then IMAGESC to plot. The X values will be on the X-axis.
% thus, colour is determined by the number of points that fall in a given square.
% sgm
if size(x,2)==2                           % using N x 2 data?
  if exist('SIZE','var'), GAUSS=SIZE; end % shift other vars along by one
  if exist('y','var'), SIZE=y; end 
  y=x(:,2); x=x(:,1);                     % split into x and y
end

if exist('SIZE','var') && ~isempty(SIZE), % make sure N is 1-by-2
  N=SIZE; if isscalar(N), N=[1 1]*N; end  
else N=[10 10];                           % default value of N is a 10-by-10 grid
end

if ~exist('GAUSS','var') GAUSS=0; end     % default no gaussian kernel 

if ~isvector(x) || ~isvector(y)           % input must (by now) be vectors
  warning('heatmap:vector','x and y should be vectors of the same length. They will be flattened.');
  x=x(:); y=y(:); 
end

[n,c]= hist3([x,y], N);                   % Do the actual bin counts
if GAUSS                                  % smooth counts using gaussian window
  if GAUSS > 2, K=floor(2*GAUSS)+1; else K=5; end % size of kernel = 2K+1 by 2K+1 grid
  K=floor(2*GAUSS)+1;                     % default K is 2 x smoothing radius
  kern = exp(-bsxfun(@plus, [-K:K].^2, [-K:K]'.^2)/GAUSS.^2) ;
  kern = kern/sum(kern(:));               % normalise to 1
  n    = conv2(n, kern,'same');           % perform averaging
end
ht=imagesc(n');                           % show heatmap
xt=get(gca,'xtick'); 
yt=get(gca,'ytick');                      % set axes scales
if all( xt==floor(xt) & xt>0 ) && all( yt==floor(yt) & yt>0 )
  makelabels = @(y) arrayfun( @(x)sprintf('%0.2g',x), y, 'uniformoutput',0);
  set(gca,'xticklabel',makelabels(c{1}(xt)), 'yticklabel',makelabels(c{2}(yt)) );
end
if nargout>0, varargout{1}=ht; end        % return handle to the image (if requested)

