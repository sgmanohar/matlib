function [U P] = uniqueProps(X)
% [U P] = uniqueProps(X)
% like unique, but also return the proportion of items in X that are of
% each type.
% sgm
nanplace=[];
if any(isnan(X))
  if ~any(X==-inf)
    nanplace=-inf;
  else
    nanplace=max(X)+1;
  end
  X(isnan(X))=nanplace;
end
[U,i,j]=unique(X);

N=length(U);
if isrow(j) j=j'; end
if ~iscolumn(j) 
  throw('matrices not yet supported');
end
P = mean( bsxfun(@eq, j, 1:N) );

if nanplace
  U(U==nanplace)=nan;
end