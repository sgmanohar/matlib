function L = logOR(X,varargin)
% log odds ratio (this is just a shortcut)
%  logOR(x,y)
%  logOR( [X,Y] )
% sgm 2016
if nargin==2
  Y = varargin{1};
else
  Y = X(:,2);
  X = X(:,1);
end

L = log( X./(1-X) ./ ( Y./(1-Y) ) );
