function u=ehazard(X,XI)
% function ehazard(X,XI)
% calculate empirical hazard function
% using ksdensity
% supply the data, and a set of evaluation points.
% returns Pi( Xi < X < Xi+dt | X > Xi) / dt
% treat nans as events that never occur
% if no ouput specified, plot the hazard function
if(~exist('XI','var'))
  XI=X;
end
if(~any(isnan(X)))
  u=ksdensity(X,XI,'function','pdf')./(1-ksdensity(X,XI,'function','cdf'));
else
  X(isnan(X))=bitmax;
  u=ksdensity(X,XI,'function','pdf', 'censoring',isnan(X))./(1-ksdensity(X,XI,'function','cdf', 'censoring',isnan(X)));
end
if(nargout==0)
  [s,s]=sort(XI);
  plot(XI(s),u(s));
  ylabel('Hazard');
end