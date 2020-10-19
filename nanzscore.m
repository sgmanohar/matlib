function Z=nanzscore(X, SDTYPE,DIM, varargin)
%  Z = nanzscore (X)
%  Z = ( X - nanmean(X) ) / nanstd( X )
%      with singleton expansion
% i.e. 
%    = bsxfun(@rdivide, bsxfun(@minus, X, nanmean(X,1)), 
%                       nanstd(X, [], 1));
% optional args SDTYPE and DIM.
% 'group', grouping variable.
% if X is a table, then apply to each numeric column
%                  group can be a column name
% sgm
if ~exist('SDTYPE','var'), SDTYPE=[]; end
if ~exist('DIM','var'),    DIM=1;     end
i=strcmp(varargin,'within');
if i>0
  GROUP = varargin{i+1};
  varargin([i i+1])=[];
else
  GROUP = [];
end

if isa(X,'table') % handle tables by recursing nanzscore on each column
  fn = fieldnames(X);
  for i=1:length(fn)
    if isa( X.(fn{i}), 'numeric' )
      if GROUP
        if isstr(GROUP) 
          if strcmp(fn{i}, GROUP) % don't zscore the grouping var!
            continue
          else
            X.(fn{i}) = nanzscore( X.(fn{i}), SDTYPE, DIM, 'group', X.(GROUP), varargin{:} );
          end
        else % numeric group
          X.(fn{i}) = nanzscore( X.(fn{i}), SDTYPE, DIM, 'group',GROUP );
        end
      else % no grouping
        X.(fn{i}) = nanzscore( X.(fn{i}), SDTYPE, DIM, varargin{:} );
      end
    end
  end
  Z = X;
  return;
end

if GROUP % recursively call nanzscore on each group
  % group is a vector of same length as size(X,1)
  ug = unque(GROUP);
  for i=1:length(ug)
    idx = GROUP==ug;
    Z(idx,:,:,:,:,:) = nanzscore( X(idx,:,:,:,:,:), SDTYPE, DIM, varargin{:} );
  end
  return
end
  


Z = bsxfun(@rdivide, bsxfun(@minus, X, nanmean(X,DIM)), nanstd(X, SDTYPE, DIM));

if isempty(Z) return; end
% allow some constant columns
badcolumns = all(isnan(Z));
if any(badcolumns)
  warning('input to nanzscore contains constant columns - ignored');
  Z(:,badcolumns) = X(:,badcolumns);
end
