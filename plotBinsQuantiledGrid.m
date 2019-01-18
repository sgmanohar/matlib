function varargout=plotBinsQuantiledGrid( X,Y, NB, varargin)
% plotBinsQuantiled2D( X,Y, NB, varargin)
% calls plotBinsQuantiled
% 
%  NB      = number of bins (determines width of sliding window, default 5)
% 'subdiv' = number of subdividsions of each bin, for gridlines
%            (determines step of sliding window, default 10)
% 'plotrargs' = a cell array of parameters to pass to 'plot' when plotting
%               grid lines.
if ~exist('NB','var')
  NB=5;    % have the number of bins been provided? if not, use 5 bins
end        % i.e. each window contains 1/5th of the dataset

[SUBDIV, plotargs] = parsepvpairs({'Subdiv', 'plotargs'},...
                                  {10, {} }, varargin{:});
% Has the number of subdivisions been provided? if not use 10 subdivisions.
% i.e. each quintile window will have 10 overlap points, giving 50 lines
                                
SM=NB/2;      % smoothing kernel width (number of subdivisions to average over)
NQ=SUBDIV*NB; % number of quantiles for the q plot = number of grid lines

% calculate gridlines by first finding y-quantiles when binning by X, 
[~, ~,bx1]=plotBinsQuantiled(Y,X, NQ, 'doplot',0, 'Bins',NB); 
[~, ~,by2]=plotBinsQuantiled(X,Y, NQ, 'doplot',0, 'Bins',NB);

if NQ>10
  gausskern = exp(-bsxfun(@plus, [-SM:SM].^2/SM, [-SM:SM]'.^2/SM));
  gausskern = gausskern/sum(gausskern(:));
  bx1=conv2(bx1, gausskern,'valid');
  by2=conv2(by2, gausskern,'valid');
end

ish=ishold();
NB=size(bx1,2);
for i=1:NB; 
  plot(bx1(:,i),by2(i,:), plotargs{:}); hold on; 
end; 
for i=1:NB; 
  plot(bx1(i,:),by2(:,i), plotargs{:}); 
end
if ~ish
  hold off;
end

if nargout==2
  varargout={bx1, by2};
end
  