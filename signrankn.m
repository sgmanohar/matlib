function [p,h, stats]=signrank(x,y,varargin)
% same as RANKSUM but handles matrices as inputs

sz=size(x);
if length(sz)>2
  x=x(:,:); y=y(:,:);
end
for i=1:size(x,2)
  xi=x(:,i); yi=y(:,i);
  bad = isnan(xi) | isnan(yi);
  xi(bad)=[]; yi(bad)=[];
  if all(isnan(xi) | isnan(yi))
    p(:,i)=nan; 
  else
    [p(:,i) h(:,i), stats{i}] = signrank(xi,yi,varargin{:});
  end
end
if length(sz)>2
  p=reshape(p,sz(2:end)); 
  h=reshape(p,sz(2:end));
  stats=reshape(stats,sz(2:end));
end
