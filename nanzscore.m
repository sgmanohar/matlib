function Z=nanzscore(X, SDTYPE,DIM)
%  Z = nanzscore (X)
%  Z = ( X - nanmean(X) ) / nanstd( X )
%      with singleton expansion
% i.e. 
%    = bsxfun(@rdivide, bsxfun(@minus, X, nanmean(X,1)), 
%                       nanstd(X, [], 1));
% optional args SDTYPE and DIM.
% sgm
if ~exist('SDTYPE','var'), SDTYPE=[]; end
if ~exist('DIM','var'),    DIM=1;     end
Z = bsxfun(@rdivide, bsxfun(@minus, X, nanmean(X,DIM)), nanstd(X, SDTYPE, DIM));

if isempty(Z) return; end
% allow some constant columns
badcolumns = all(isnan(Z));
if any(badcolumns)
  warning('input to nanzscore contains constant columns - ignored');
  Z(:,badcolumns) = X(:,badcolumns);
end
