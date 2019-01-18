function varargout = corrplot(X,varargin)
%     corrplot(X, ... )
%        correlate all pairs of columns of X
%     corrplot(X, Y, ... )
%        correlate columns of X with columns of Y
%     [r, p] = corrplot( ... )
%         return matrix of correlation coefficients and p-values using Matlab's 'CORR'
% 
% plot the correlations in a grid of subplots.
%
%   'histfun':     Histogram function, used for the diagonal subplots. 
%                  This will be called with each column of X, as a vector.
%                  This is only used if Y is not specified.
%   'scatterfun':  Scatter function, used for the off-diagonal subplots. 
%                  This will be called with each pair of columns (X1,X2).
%   'colnames':    Names of columns for title.
%   'colnames2':   if Y is specified, gives column names for Y for title.
%   'list':        if true, display a list of significant pairs of columns
%   'corrargs':    a cell containing arguments passed to 'corr', e.g. {'type','spearman'} 
%   'doPlot':      deafult true; if false, don't actually draw anything
% sgm 2015


if nargin>1 && isnumeric(varargin{1}) % allow correlating columns of two different matrices
  Y=varargin{1}; varargin(1)=[]; 
  yisx=false;
else
  Y=X; yisx=true;
end


[  SCATTERFUN,    HISTFUN,   COLNAMES , COLNAMES2,   TRANSFORM,    LIST, ALPHA,  CORRARGS    DOPLOT     TYPE] = parsepvpairs( ...
  {'scatterfun', 'histfun', 'colnames', 'colnames2', 'transform', 'list', 'alpha','corrargs' ,'doplot','type'}, ...
  {@scatterRegress, @hist ,  []  ,       [],           @(x)x,      false ,  0.05 , {} ,        1,      'pearson'}, varargin{:});

if yisx && ~isempty(COLNAMES2), error('colnames2 should only be used if correlating two matrices.'); end
if yisx, COLNAMES2=COLNAMES; end; 
if isequal(SCATTERFUN,@scatterRegress) && strcmpi(TYPE,'spearman') % pass appropriate correlation type to scatterRegress
  SCATTERFUN = @(varargin)scatterRegress(varargin{:},'pearson',0);
end
if ~strcmp(TYPE,'pearson') % add correlation type to the "corr" parameters if necesary
  CORRARGS=[CORRARGS {'type',TYPE}];
end

Nx=size(X,2);
Ny=size(Y,2);
if DOPLOT
  for i=1:Nx
    if yisx % if autocorrelating, then only use below-diagonals
      j0=i;
    else
      j0=1;
    end
    for j=j0:Ny
      subplot(Ny,Nx, i+(j-1)*Nx);
      if j~=i || ~yisx
        SCATTERFUN(X(:,i),Y(:,j));
        if ~yisx
          if ~isempty(COLNAMES)
            xlabel(COLNAMES{i});
          end
          if ~isempty(COLNAMES2)
            ylabel(COLNAMES2{j});
          end
        end
      else % j==i
        if yisx
          HISTFUN(X(:,i));
          if ~isempty(COLNAMES),
            title(COLNAMES{i});
          end
        end
      end
    end
  end
end
% run correlations
if nargout>0 || LIST
  % calculate correlation coefficients and p values using CORR
  [varargout{1:2}] = corr(X,Y,CORRARGS{:});
end


if LIST 
  if isempty(COLNAMES2), COLNAMES2=arrayfun(@(x)sprintf('C%g',x),1:size(Y,2),'uni',0); end
  sig = varargout{2}<ALPHA; % which pairs are significant?
  [i j] = find(sig); 
  table = {};
  for ii=1:length(i) % for each significant pair
    if yisx                        % when only one matrix is used,
      if j(ii)>=i(ii); continue; end % only check nondiagonal correlations, and ignore one triangle.
    end
    table = [table; { COLNAMES{i(ii)}, COLNAMES2{j(ii)}, varargout{1}(i(ii),j(ii)), varargout{2}(i(ii),j(ii)) } ];
  end
  displaytable(table, {'x','y','r','p'}, [35 25 11 11]);
end