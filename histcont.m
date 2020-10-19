function varargout=histcont(Y, bins, width, varargin)
% Histogram of continuous data -- Sliding window.
%   - calculate a binned probability density
%   - works on n-dimensional data, using dimension 1.
% function [N, X] = histcont(X, bins, width, ...)
%   uses N = HIST(Y) to count elements of Y in sliding window of width 1/10th 
%               of the range of Y. Defaults to 200 window positions.
%               as with HIST, if Y=matrix, then work on each column.
%   BINS  = number of points to return. Like HIST, you can instead specify
%           the bin centres as a vector. 
%           default 200, or number of samples whichever is smaller.
%   WIDTH = the width of each bin as a proportion of the whole range.
%           default = 0.1
% 
% note that HISTCONT differs from HIST in that 
%   1) the returned N-values are the count per unit of X
%   2) the integral of N dX will equal the number of non-NaN values in the
%      sample; i.e. the denominator of the density counts NaNs
%
% params: 'Plot'              'Bar' (like original HIST), 'Line' (default), 
%                             'Area', 'None', 'CDF' (cumulative density)
%         'Count'  t/f        scale height to the raw count per bin, rather 
%                             than "density per unit X" (default false) 
%         'Proportion' t/f:   calculate proportion, rather than count (default false).
%         'Plotargs':         a cell array, whose elements get passed to "plot" 
%         'Support':          what range of x-values to use. 'range': min:max; (default)
%                             'sd': -2sd to +2sd around mean; 
%                             'quantile': 5th to 95th percentile  
% SGM 2014

Plot = enum({'LINE','BAR','CDF','AREA', 'NONE'});
PLOT = Plot.LINE; % do a bar histogram? otherwise a line.
[     COUNT,     PROP,        plotargs,    PLOT ,  SUPPORT  ] = parsepvpairs( ...
  {  'Count' ,  'Proportion', 'plotargs', 'plot', 'support' }, ...
  {   false ,    false ,      {},       nargout==0, 'range' }, ...
  varargin{:}); 

if isstr(PLOT)
  PLOT = Plot.getIndex(PLOT);
elseif PLOT==0
  PLOT = Plot.NONE;
end

if isrow(Y) Y=Y'; end % make into columns
if ~exist('width','var') || isempty(width), width = 0.1; end % width of a bin, as a proportion of range
if ~exist('bins','var') || isempty(bins),  % number of bins
  bins  = floor(size(Y,1)*(1-width)); % bin for every data point!
  bins=min(bins,200); % default: no more than 200 bins!
end 
if ~iscolumn(bins), bins=bins'; end % ensure column
if bins*width<1, warning(['the histogram bins are not wide enough to overlap!' ...
    'width should be greater than 1/bins']); end

% calculate the support of the distribution, based on the whole dataset
% across all columns.
switch SUPPORT
  case 'range'
    range = [min(Y(:)),max(Y(:))]; % work from max to min over all datapoints
  case 'sd'
    rangesize = [-2 +2];
    range = nanmean(Y(:)) + rangesize * nanstd(Y(:)); 
  case 'quantile'
    rangesize = [0.05 0.95];
    range = quantile(Y(:),rangesize);
end

binwidth = diff(range)*width;  % physical bin width
if length(bins)==1 % create edges and centres of bins
  binL = linspace(range(1), range(2)-binwidth, bins)'; % linear spaced bin starts
  binR = linspace(range(1)+binwidth, range(2), bins)'; % bin ends
  binC = (binL+binR)/2;                                % bin centres
else % user provided bin centres manually: ignore width parameter
  binL = bins-binwidth/2; 
  binR = bins+binwidth/2; % create edges between centres.
  binC = bins;
end
% calculate ( Y > binL & Y < binR )
% for every value of Y, and sum over Y's columns.
shdim1 = [2 3 1]; shdim2 = [2 1]; 
if ndims(Y)==3, shdim1=[2 3 4 1]; shdim2=[3 1 2]; 
elseif ndims(Y)==4, shdim1=[2,3,4,5,1]; shdim2=[4,1,2,3];
elseif ndims(Y)==5, shdim1=[2 3 4 5 6 1]; shdim2=[5 1 2 3 4]; end

N=permute( squeeze(sum( bsxfun(@ge, Y, permute(binL, shdim1 )) ...
                      & bsxfun(@lt, Y, permute(binR, shdim1 )),1) ) , shdim2); 

if isempty(plotargs) && ndims(N)>2 %% default area if more than 2 dimensions
  plotargs = {'area',1,'dostats',0};
end


if ~COUNT
  N=N/binwidth; % convert from "count per bin" into "Count per unit of X"
end
if PROP % proportions rather than counts? divide by number of samples
  N=bsxfun(@rdivide, N, sum(~isnan(Y)) ); % pray these are the same dimensions
end
X=binC;
if PLOT == Plot.CDF 
  N=cumsum(N,1); 
end

if PLOT ~= Plot.NONE,  % if doing plot: 
  if PLOT == Plot.BAR, h=bar(X,N);  % original bar histogram?
  elseif PLOT == Plot.AREA
    %h=area(X,N,'facealpha',0.8); % stacked areas?
    h=area(X, [N(:,1) diff(N,[],2) ] ,'facealpha',0.7); % unstacked!
  else % line plot
    if ndims(N)==2, 
      h=plot(X,N, plotargs{:});   % continuous line histogram
    else
      if ndims(N) < 4
        h=errorBarPlot(permute(N,[2 1 3]), 'xaxisvalues',X, plotargs{:});
      else
        plotn(permute(N,[2,1,3,4,5,6,]), ...
          'plotfun',@(x)errorBarPlot(x,'xaxisvalues',X,plotargs{:}) , ...
          'threed',1,'fixshape',1) ;
      end
    end
  end
end
if nargout==1, varargout={N}; 
elseif nargout==2, varargout={N X};
end

