function [x,y,hh]=plotBins(X,Y, N, varargin)
% [x,y,hh]=plotBins(X,Y, N, varargin)
%
% Plot a the value of Y for each value of X,
% but bin X according to quantiles.
% uses nanmean to calculate the Y value for each X-bin,
% and  nanmean to calculate the X bin-centres.
% 
%   varargin: is passed to 'errorbar'.
%   returns: x = bin centres, y = mean of Y in each bin
%            hh = figure handle
% sgm
if(~exist('N','var'))  || isempty(N)
  N=5;
  warning('quantile','using %g quantiles',N);
end

if any(strcmpi(varargin, 'xerror')), 
  XERROR=1; varargin(strcmpi(varargin, 'xerror'))=[]; 
else XERROR = 0; end

fprintf('%d samples per point\n', ...
  floor(length(X)/N));

SE=1; % stdev = 0, stderr = 1;

q=quantile(X,([0:N])/N);
qcentre=[0.5:(N+0.5)]/N;
q(end)=inf;
for(i=1:N)
  f=X>=q(i) & X<q(i+1);
  % take n as min number of trials in the corresponding bin for RT1 and RT2
  n=sum(f);
  if(n==0)
    % a quantile is repeated - too many points for one bin
    warning('too many quantiles')
  end
  y(i)=nanmean(Y(f));
  x(i)=nanmean([q(i),q(i+1)]);
  x(i)=quantile(X,qcentre(i));
  ys(i)=nanstd(Y(f));
  xs(i)=nanstd(X(f));
  if(SE) ys(i)=ys(i)/sqrt(n); end % standard error
end
%x(end)=quantile(X,1);

held=ishold();
if XERROR % manual with X and Y error bars
  plot(x,y,varargin{:});
  for(i=1:N)
    hold on;
    plot(x(i)+[-1 1]*xs(i), [1 1]*y(i), '-'); % x errorbar
    plot([1 1]*x(i), y(i)+[-1 1]*ys(i), '-'); % y errorbar
  end
else % automatic- use errorbar
  hh=errorbar(x,y,ys, varargin{:});
end

if(~held), hold off;end;
