function [params, values, hplot, means]=fitSigmoid(Y, varargin)
% [params, values] = fitSigmoid(Y, varargin)
% Fit a sigmoid curve to each column of Y.
% Fitting is done by least squares. 
% 
% 'X' specifies X values corresponding to each Y. If X is not specified, 
%         X increases linearly, i.e. X = [ -(N-1)/2 : (N-1)/2 ]' for each column.
% 'model' should either a number: one of 
%     1. sig( p1 * (x - p2) ) * dY + Y0    : basic sigmoid with slope and intercept 
%     2. sig( p1 * (x - p2) ) + p3         : model pedestal 
%     3. sig( p1 * (x - p2) ) * p3 + Y0    : model scaling
%     4. sig( p1 * (x - p2) ) * p3 + p4    : model both pedestal and scaling
% 'Iter' number of iterations when searching for optimal fit. Default=10
% 'Xmodel' is the values of X that will be used to return modelled values. for
%         example if you wanted to know the modelled Y value when X = 0, you can
%         put 0 here.  Note that these X values are used for plotting if
%         specified, and they default to a range of 50 points between the
%         minimum and maximum X.
% 'Plot' could be 0 = No plot; 1 = overlaid plots, 2 = individual subplots,
%         3 = show group means (mean taken on dimension 2 of data). This plot
%         will show standard error of mean values, model according to the mean
%         of parameter estimates, and separate model fitted to the mean data.
% 'optimset' is a structure created with the OPTIMSET function, specifying any
%         parameters for fminsearch
% 'function' is a different function to use for the fitting. The default
%         function is the logistic sigmoid: Y ~ 1/(1+exp(-x)). 
%         The first parameter to the function is X, the subsequent parameters
%         are fitted to the data. The parameter 'model' is discarded.
%         e.g. 'function', @(x,p1,p2) p1./(1+exp(-p2*x))
%
% returns:
%     parameter estimates for p1...p4 
%     values corersponding to Xmodel
%     hplot is the handles of the plotted objects, column 1 = model, column 2 =
%         data, column 3 = dotted line for bias, column 4 = horizontal 50%
%         meridian.
% 
% sanjay manohar 2015

if size(Y,1)==1, Y=Y'; end % row should be column!
ND     = size(Y,1); % num data points

% parse parameters
[  MODEL,   X,          ITER,      xmodel,   PLOT ,  OPTIM,   FITFUNC] = parsepvpairs(...
  {'model', 'X',          'Iter' , 'XModel', 'Plot', 'optimset', 'function'} ...
 ,{4, [1:ND]'-(ND+1)/2 ,    10,      [],      2,     [],          [] } ...
 ,varargin{:});

% make sure X is same size as Y by expanding singleton dimensions
X = bsxfun(@times, ones(size(Y)), X); 

% collapse across dimensions 2 upwards, to give a matrix
oshape = size(Y); 
oY     = Y; % original shape
oX     = X; 
Y      = reshape(Y, ND, []); 
X      = reshape(X, ND, []);
NS     = size(Y,2); % num columns (after reshaping)

% basic sigmoid function from 0 to 1, taking a value of 0.5 when x=0
sigma  = @(x) 1./(1+exp(-x)); 
minY   = min(Y(:));
Yrange = max(Y(:))-minY; 
Xrange = max(X(:))-min(X(:)); 
minX   = min(X(:)); maxX = max(X(:)); 
if isempty(xmodel), 
  xmodel = linspace(minX, maxX, 50);
end
if isempty(OPTIM), OPTIM = optimset('display','off'); end 
nPlot = 1+floor(sqrt(NS));   % how many subplots required, to tile all values
if PLOT>2 && nPlot.^2==NS, nPlot=nPlot+1; end % space for the grand average?

prevhold=ishold();
for i=1:NS % for each column
  % function to optimise
  if isempty(FITFUNC)
    if MODEL==1 % simple 2 parameter fit, slope and intercept
      parnames = {'slope','bias'};
      predict  = @(x,p) minY+Yrange.* sigma( -p(1) * (x-p(2)) ) ;
      p0e      = [ Yrange/Xrange mean(X(:)) ];
      p0 = [ -1 -1 ; 1 1 ]; % initial parameters?
    elseif MODEL==2 % include pedestal
      parnames = {'slope','bias','pedestal'};
      predict = @(x,p) Yrange*sigma( -p(1) * (x-p(2)) )  + p(3);
      p0 = [ -1 -1 -1; 1 1 1 ]; % initial parameters
    elseif MODEL==3
      parnames = {'slope','bias','sigheight'};
      predict = @(x,p) minY + sigma( -p(1) * (x-p(2)) ) * p(3);
      p0 = [ -Yrange/Xrange -Xrange  0
        -Yrange/Xrange  Xrange  Yrange*2  ]; % initial parameters
    elseif MODEL==4 % similar but with scaling of height as well as pedestal
      parnames = {'slope','bias','sigheight','pedestal'};
      predict = @(x,p)  sigma( -p(1) * (x-p(2)) )  * p(3) + p(4);
      p0 = [ -Yrange/Xrange -Xrange 0         0;
        Yrange/Xrange  Xrange Yrange*2  minY*2]; % initial parameters
    end
  else % provided a custom function 
    parnames = arrayfun(@(x)sprintf('par%g',x),[1:nargin(FITFUNC)-1],'unif',0);
    predict  = @(x,p) call(FITFUNC,[{x} num2cell(p)]);
    p0       = repmat([-1;1],1,nargin(FITFUNC)-1);
  end
  % squared errors of prediction to actual SIP
  sse     = @(p)  nansum( ( Y(:,i) - predict(X(:,i),p) ).^2 );
  [p1(i,:) tsse(i)] = fminsearchs( sse, p0 , ITER , OPTIM);
  ymodel = predict(xmodel, p1(i,:)) ;
  ymodel_a(i,:) = ymodel;
  if PLOT
    if PLOT>1,     subplot(nPlot,nPlot,i); end
    h(1,i)=plot( xmodel, ymodel );  % model
    hold on;
    h(2,i)=plot(   X(:,i) , Y(:,i) ,'o'); % data
    if p1(i,1)>minX & p1(i,1)<maxX
      h(3,i)=plot( p1(i,2) * [ 1 1 ], ylim ,':'); %  bias line
    end
    if i==1 || PLOT>1 % plot baseline
      h(4,i)=plot( xlim, (minY+Yrange/2)*[1 1] ,':');
    end
    % set(gca,'xtick',X(:,i), 'xticklabel',X(:,i));
    if PLOT>1, 
      hold off; 
      title(['params = [' sprintf('%0.2g ',p1(i,:)) ']']);
    end
    drawnow
  end
end
params=reshape(p1', [size(p1,2), oshape(2:end) ]); % restore to original size
values = reshape(ymodel_a', [size(ymodel,2), oshape(2:end)]); 

perm = [2,1,3,4,5,6]; % permutation of dimensions so errorBarPlot can take the mean on dimension 2
if PLOT>2 && size(oY,2)>2 % are there multiple columns, and are we plotting grand mean?
  subplot(nPlot,nPlot,nPlot*nPlot);
  errorBarPlot( permute(oY,perm)  , 'xaxisvalues',sq(mean(permute(oX,perm))) );
  means.ymodel_meanp = predict( xmodel, mean(p1) );
  hold on;  plot(xmodel, means.ymodel_meanp,':'); 
  % now try fitting the grand mean curve?
  sse     = @(p)  sum( ( mean(Y,2) - predict(X(:,i),p) ).^2 );
  [means.p1 means.tsse] = fminsearchs( sse, p0 , ITER , OPTIM);
  means.ymodel_meany = predict(xmodel, means.p1) ;
  plot(xmodel, means.ymodel_meany,'r:'); 
  legend({'mean','mean(params)','model mean(y)'});
end
if PLOT>1, makeSubplotScalesEqual(nPlot,nPlot, [1:NS, nPlot*nPlot]); end
if ~prevhold, hold off; end

function y=call(fun, params)
% shortcut to call a function with an array as parameters
% only needed because you can't do num2cell( ... ){:} 
y=fun(params{:});