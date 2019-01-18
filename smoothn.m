function X = smoothn(X, n, method, varargin)
% X = smoothn( X, n, method)
% X = smoothn( DIM, X, n, method )
% 
% runs SMOOTH on each column of the data X.
% syntax - same as SMOOTH.
% also allows SMOOTHN( DIMENSION, X, ...)
% default N=5.
% if 'method' is 'nan', or if there's no smooth file, then uses nanconv
% with a gaussian filter.
% sgm 2014

%%%%%%%% parase params
if ~exist('n','var') || isempty(n)   % if smoothing range not specified
  n=5;                               % select default value
end
if(n==0) warning('smooth called with zero width - nothing done!'); return; end
if isscalar(X) && X-floor(X)==0      % is first arg a dimension? (integer)
  DIM = X; 
  X=n;                               % bump along parameters
  if exist('method','var'), n=method; method=[]; else n=[]; end
  if ~isempty(varargin), method=varargin{1}; varargin=varargin(2:end); end
else DIM=1; end

if ~exist('method','var') || isempty(method)  % if kernel not specified, select a default
  if regexp(help('smooth'), 'gauss') % if we have the newer SMOOTH function
    method='gauss';                  % use gaussian kernel
  else                               % otherwise use moving average
    method='moving';
  end
end

if strcmpi(method,'nan') || ~exist('smooth','file') % create a gaussian kernel but use nanconv
  smoother = @(x,n,m) nanconv( x,exp(-(3*[-n:n]'/n).^2)/sum(exp(-(3*[-n:n]'/n).^2)) ,'nonanout','edge','1d'); 
else smoother = @smooth;
end
if ~exist('n','var') || isempty(n)   % if smoothing range not specified
  n=5;                               % select default value
end

%%%%%%%% do smooth
sz0=size(X);                          % remember original shape of data
ND=length(sz0);                       % number of dimensions
if DIM~=1,                            % if not going down columns, 
  otherdims = [1:ND]; otherdims(DIM)=[]; 
  X=permute(X, [DIM otherdims]);      % bring DIM to columns
  unpermute = [ [2:DIM] 1 [DIM+1:ND] ]; % and keep track of how to 
end                                   % restore the position of DIM.
sz=size(X);
X=reshape(X,sz(1),[]);                % collapse later dimensions
for i=1:size(X,2)                     % for each column of X
  X(:,i) = smoother(X(:,i), n, method); % smooth it
end
X=reshape(X, sz);                     % return data to original shape
if DIM~=1
  X=permute(X, unpermute);
end
if sz0 ~= size(X), error('dimensional mishap in smoothn'); end