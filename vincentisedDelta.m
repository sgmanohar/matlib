function [x,y,h]=vincentisedDelta(RT1, RT2, N, varargin)
% [x,y]=vincentisedDelta(RT1, RT2, N, varargin)
%
% draw the vincentised delta plot of the two RT distributions
% with N quantile bins. RT1 is the longer group or RTs
% The delta plot is
% X coordinate = mean of i-th quantile of RT1 and i-th quantile of RT2
% Y coordinate = RT difference between mean RT in i-th quantile bin of RT1 
%                and mean RT in i-th quantile of RT2
%                (RT1 - RT2)
%
% default N = 5
% extra parameters are passed to 'plot'.
% 
% if you specify 2 outputs x and y, nothing is plotted, the bins are just
% calculated. If a third output arg is specified, the figure handles for
% the line, x-error, and y-error, are returned.
%
% (Ridderinkhof 2002, Ratcliffe 1979)
% sanjay manohar

if(~exist('N','var')) 
  N=5;
  warning('quantile','using %g quantiles',N);
end
fprintf('%d samples per point in group 1, %d in group 2', ...
  floor(length(RT1)/N), floor(length(RT2)/N));

q1=quantile(RT1,([1:N+1]-1)/N);
q2=quantile(RT2,([1:N+1]-1)/N);
for(i=1:N)
  f1=RT1>q1(i) & RT1<q1(i+1); % filter for current quantile bin in RT1
  f2=RT2>q2(i) & RT2<q2(i+1); % filter for current quantile bin in RT2
  % take n as min number of trials in the corresponding bin for RT1 and RT2
  n=min(sum(f1) + sum(f2));
  y(i)=nanmean(RT1(f1)) - nanmean(RT2(f2));
  alldiff=bsxfun(@minus, RT1(f1), RT2(f2)');
  ys(i)=nanstd(alldiff(:)) / sqrt(n);
  x(i)=nanmean([q1(i),q1(i+1),q2(i),q2(i+1)]);
  xs(i)=nanstd([q1(i),q1(i+1),q2(i),q2(i+1)]) / 2;
end

if(nargout~=2)
  h(1)=plot(x,y,varargin{:});
  held=ishold();
  for(i=1:N)
    hold on;
    h(2)=plot(x(i)+[-1 1]*xs(i), [1 1]*y(i), '-'); % x errorbar
    h(3)=plot([1 1]*x(i), y(i)+[-1 1]*ys(i), '-'); % y errorbar
  end
  xlabel('RT');
  ylabel('delta RT between RT1 and RT2');
  if(~held) hold off;end;
end