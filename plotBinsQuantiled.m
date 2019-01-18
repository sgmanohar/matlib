function varargout=plotBinsQuantiled(X,Y, N, varargin)
%  h = plotBinsQuantiled ( X, Y, N ... )
%         plot the quantile lines of Y, for each quantile of X.
% 
%  [h, xbin, ybin] = plotBinsQuantiled ( X, Y, N ... )
%         return the quantiles of x and y 
%
% N is number of quantiles, can be a vector [NX NY] to specify number of
% quantiles for X and Y separately
% 
% 'Bins' specifies the how big a bin should be. 
%  E.g. if N=100 and Bins = 5, then the bin width will be 1/5 of the 100
%  quantiles, i.e. 20% of the X quantiles. Default = 5
%
% 'DoPlot' - if 0, then don't plot anything
%
% returns 
%    xbin:   a vector length N(1)
%    ybin:   a matrix N(1) x N(2) matrix with quantiles of y for each x-bin
% 
% if X/Y are matrices, then quantiles are calculated for each column, 
% and then the values for each quantile bin are averaged across all 
% rows and columns that match the binning.
% 

if(~exist('N','var'))
  N=100;
end
if(length(N)==1) NX=N; NY=N; 
else NX=N(1); NY=N(2);
end

i=find(strcmpi(varargin,'Bins')); if any(i)
  NBINS = varargin{i(1)+1};
  varargin( [ i(1) i(1)+1 ] ) = []; 
else NBINS = 5;
end

qx=quantile(X,linspace(0,1,NX+1)); % create quantiles of X
WIDTHX = floor(NX/NBINS); % number of quantiles in one bin
WIDTHY = floor(NY/NBINS); % number of quantiles in one bin

% go through bins of x, and select y-values in each bin.
% then bin these y values according to their quantiles.
% the x-quantiles are stored in xbin(i), and the corresponding y-quantiles are
% stored in ybin(i,j).
ybin=nan(NX-WIDTHX+1, NY-WIDTHY+1);    % y values in each x-bin
for(i=1:NX-WIDTHX+1)                   % for each x-bin
  filterx=X>=qx(i) & X<qx(i+WIDTHX);   % select trials in this x-range
  yi = Y(filterx);                     % get their y values
  qy=quantile(yi, linspace(0,1,NY+1)); % and find their quantiles 
  for(j=1:NY-WIDTHY+1)                    % for each of theese y quantile bins
    filtery=yi>=qy(j) & yi<qy(j+WIDTHY);  % select the values in this bin.
    if ~isempty(yi(filtery))              % if there are values, 
      ybin(i,j) = nanmean( yi(filtery) ); % calculate y bin centre value
    else ybin(i,j)=nan; end               % otherwise give nan.
  end
  xbin(i)=qx(i+floor(WIDTHX/2));       % store x bin centre
end

if ~any(strcmpi(varargin,'doplot'))
  h=plot(xbin, ybin, varargin{:});
else h=nan;
end

if nargout==1,     varargout={h};
elseif nargout==2, varargout={h xbin}; 
elseif nargout==3, varargout={h xbin ybin}; 
end