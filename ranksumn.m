function [p,h, stats]=ranksumn(x,y,varargin)
% same as RANKSUM but handles matrices as inputs
szx=size(x); szy = size(y);
if size(szx,1)~=size(szy,1) %|| any(szx~=szy)
  error('X and Y must have same size');
end
x=reshape(x, szx(1), []);
y=reshape(y, szy(1), []);
for i=1:size(x,2)
  xi=x(:,i); yi=y(:,i);
  xi(isnan(xi))=[];
  yi(isnan(yi))=[];
  [p(:,i) h(:,i), stats(i)] = ranksum(xi,yi,varargin{:});
end
if length(szx)>2 % more than 2 dimensions? 
  p=reshape(p, szx(2:end) );
  h=reshape(h, szx(2:end) );
end
