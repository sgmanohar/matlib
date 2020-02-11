function [chain, s2chain, res, olsfit] = mcmc_lm(y,X, varargin)
% [chain, s2chain, res, ols] = mcmc_lm( y, X , ... )
% [...]                      = mcmc_lm( Y , ... )
% convenience function to call MCMC toolbox with a linear model. 
% if X specified, then fit y=Xb. If X not specified, Y should be a
% n-dimensional array, where each dimension is a predictor. 
% (like for 'anovanTable').
% 
%   'modelfun': function @(X, beta) to give estimate of y.
%               default: y = X * beta'
%   'ones':     insert a column of ones (intercept) if it doesn't exist 
%               default: true
%   'plot': 


if ~isnumeric(X) && ~isvector(y) % accept N-dimensional arrays for Y, if no X specified
  y=depivot(y);  % convert the dimensions of y into predictors
  X=y(:,1:end-1); % take first columns as X (predictors)
  y=y(:,end); 
end

[ modelfun    INSERT_ONES , PLOT   ] = parsepvpairs(...
{ 'modelfun' ,'ones' ,    'plot'   }, ...
{  [],         true         false   }, varargin{:} );

if size(X,2)>size(X,1), X=X'; end
if isrow(y), y=y'; end;
[nD, nP] =size(X); 
varx=var(X);
if varx(1)==0 % is it a column of ones on the left?
  X=X(:,[2:end 1]); % put it at the end.
end
if all(varx==0) && INSERT_ONES
  X=[X ones(nD,1)]; 
end
data.xdata   = X;
data.ydata   = y;
if isempty(modelfun)
  modelfun     = @(X,b) X*b';
end
model.ssfun  = @(b,data) (data.ydata - modelfun(data.xdata, b')).^2;
model.modelfun = modelfun;
% apply ordinary least squares to get estitletimate
[olsfit.b, olsfit.var, olsfit.t, olsfit.res] = ols(y,X,eye(nP));
% this is how good we expect our model to fit!
model.sigma2 = mean(olsfit.res.^2); % estimate error variance
% parameter array { { 'p1', ols.b(1) }, 
%                   { 'p2', ols.b(2) } ... }
params = arrayfun(@(x){ sprintf('p%g',x), olsfit.b(x) },[1:nP],'uni',0)';
options.nsimu = 5000;
options.updatesigma = 1;
% options.qcov = ??
  
[res, chain, s2chain] = mcmcrun(model,data,params,options);

if PLOT % plot the data against each predictor, along with predictions
  out=mcmcpred(res, chain, [], X, modelfun);
  for i=1:nP; subplot(1,NP,i)
    o=out; o.data=X(:,i); 
    mcmcpredplot(o);
    hold on; hold on; plot(X(:,i), y,'s'); hold off;
    xlabel(params{i}{1});
  end
end
  
