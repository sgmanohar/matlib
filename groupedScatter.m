function ym=groupedScatter(X,Y,N,q,r)
% groupedScatter(X,Y,N)
% take x and y coords, and plot a kind of
% scatter where points are grouped according to quantiles
% and mean Y values plotted.
% N is number of quantile bins, if X is continuous

if(exist('q','var')) 
  if(exist('r','var'))
    ym=groupedScatter2(X,Y,N,q,r); return;
  else
    ym=groupedScatter2(X,Y,N,q); return;
  end
end;

if ~exist('N','var') 
  u=unique(X);
  if(length(u)<=length(X)/2) % at least 2 points per x-coord?
    N=length(u);
    xbins=[u, u+eps]; % treat X as discrete.
  else
    N=length(x)/10
  end
end
%calculate percentiles if X not discrete or if N specified
if ~exist('xbins','var')
  for(i=1:N)
    xbins(i,1)=prctile(X,100*(i-1)/(N-1));
    xbins(i,2)=prctile(X,100*i/(N-1));
  end;
end;
for(i=1:N)
  xm(i)=mean(xbins(i,:));
  yi=Y( (X>=xbins(i,1)) & (X<xbins(i,2)) );
  ym(i)=mean(yi);
  ys(i)=std(yi);
  yn(i)=length(yi);
end;
errorbar(xm,ym,ys./sqrt(yn));



function ym=groupedScatter2(X,Y,criteriaX,criteriaZ, individualX)
% groupedScatter(X,Y,N)
% take x and y coords, and plot a kind of
% scatter where points are grouped according to quantiles
% and mean Y values plotted.
% criteriaX specify how the data is grouped into x-values;
% criteriaZ specify how to split the data into series lines.
%
% if criteriaX is a single integer, this specifies the number 
% of quantile bins to split X into. If individualX is true, 
% then each series is independently binned, otherwise the data
% is quantile-binned en bloc (resulting in different numbers of 
% datapoints in each bin, though).

crX=unique(criteriaX);
crZ=unique(criteriaZ);
if(~exist('individualX','var')) individualX=0;end;

if(size(crX)==size(X)) % criteriaX bins the data.
  for(i=1:length(crZ))
    for(j=1:length(crX))
      xm(i,j)=crX(j);
      yi=Y(criteriaZ==crZ(i) & criteriaX==crX(j));
      ym(i,j)=nanmean(yi);
      ys(i,j)=nanstd (yi);
      yn(i,j)=length (yi);
    end
  end
elseif length(crX)==1 % criteriaX specified number of quantile bins
  N=crX(1);
  for(i=1:length(crZ))
    for(j=1:N)
      if(individualX)
        xbins(i,j,1)=prctile(X(criteriaZ==crZ(i)), 100*(j-1)/(N));
        xbins(i,j,2)=prctile(X(criteriaZ==crZ(i)), 100*(j)/(N));
      else
        xbins(i,j,1)=prctile(X, 100*(j-1)/(N));
        xbins(i,j,2)=prctile(X, 100*(j)/(N));
      end
      xm(i,j)=nanmean(xbins(i,j,:));
      yi=Y( (X>=xbins(i,j,1)) & (X<xbins(i,j,2)) & (criteriaZ==crZ(i)) );
      ym(i,j)=nanmean(yi);
      ys(i,j)=nanstd(yi);
      yn(i,j)=length(yi);
    end;
  end
else
  error('criteriaX must be of the same size as data, or a single value');
end

errorbar(xm',ym',ys'./sqrt(yn'));


