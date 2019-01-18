function varargout=histcont(Y, bins, width, varargin)
% HISTCONT Histogram of continuous data -- Sliding window.
%   - calculate a binned probability density
% function [N, X] = histcont(X, bins, width, ...)
%   N = HIST(Y) count elements of Y in sliding window of width 1/10th 
%               of the range of Y. Defaults to 10 positions.
%               as with HIST, if Y=matrix, then work on each column.
%   BINS  = number of points to return. Like HIST, you can instead specify
%           the bin centres as a vector. 
%           default 10.
%   WIDTH = the width of each bin as a proportion of the whole range.
%           default = 0.1
% 
% note that HISTCONT differs from HIST in that 
%   1) the returned N-values are the count per unit of X
%   2) the integral of N dX will equal the number of non-NaN values in the
%      sample; i.e. the denominator of the density counts NaNs
%   3) currently configured to plot a line histogram. If you prefer the
%      original bar histograms, edit the line at the top of the code!
%
% params: 'Bar'   true/false: draw bars instead of line
%         'Count'             scale height to the raw count per bin, rather 
%                             than "per unit X"
%         'Proportion' t/f:   calculate proportion, rather than count.
%         'Plotargs':         a cell array, whose elements get passed to "plot" 
% SGM 2014

BAR = false; % do a bar histogram? otherwise a line.
[BAR, COUNT, AREA, PROP, plotargs ] = parsepvpairs( ...
  { 'Bar', 'Count' , 'Area', 'Proportion', 'plotargs'}, ...
  {  false,  false , false , false , {}      }, ...
  varargin{:}); 



if isrow(Y) Y=Y'; end % make into columns
if ~exist('width','var'), width = 0.1; end % width of a bin, as a proportion of range
if ~exist('bins','var'),  % number of bins
  bins  = floor(size(Y,1)*(1-width)); % bin for every data point!
  bins=min(bins,200); % default: no more than 200 bins!
end 
if ~iscolumn(bins), bins=bins'; end % ensure column
if bins*width<1, warning(['the histogram bins are not wide enough to overlap!' ...
    'width should be greater than 1/bins']); end

range = [min(Y(:)),max(Y(:))]; % work from max to min over all datapoints
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
  plotargs = {'area',1};
end


if ~COUNT
  N=N/binwidth; % convert from "count per bin" into "Count per unit of X"
end
if PROP % proportions rather than counts? divide by number of samples
  N=bsxfun(@rdivide, N, sum(~isnan(Y)) ); % pray these are the same dimensions
end
X=binC;
if     nargout==0, 
  if BAR, h=bar(X,N);  % original bar histogram?
  elseif AREA
    %h=area(X,N,'facealpha',0.8); % stacked areas?
    h=area(X, [N(:,1) diff(N,[],2) ] ,'facealpha',0.7); % unstacked!
  else
    if ndims(N)==2, 
      h=plot(X,N, plotargs{:});   % continuous line histogram
    else
      h=errorBarPlot(permute(N,[2 1 3]), 'xaxisvalues',X, plotargs{:}); 
    end
  end
elseif nargout==1, varargout={N}; 
elseif nargout==2, varargout={N X};
end