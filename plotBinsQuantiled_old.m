function h=plotBinsQuantiled(X,Y, N, varargin)
% plotBinsQuantiled ( X, Y, N ... )
% plot the quantile lines of Y, for each quantile of X.
% 
% N is number of quantiles, can be a vector [NX NY] to specify number of
% quantiles for X and Y separately

if(~exist('N','var'))
  N=5;
end
if(length(N)==1) NX=N; NY=N; 
else NX=N(1); NY=N(2);
end
qxO=quantile(X,0);
for(i=1:NX)
  qx=quantile(X,i/NX);
  filterx=X>qxO & X<qx;
  qyO=quantile(Y(filterx),0);
  for(j=1:NY)
    qy=quantile(Y(filterx),j/NY);
    filtery=Y>qyO & Y<qy;
    % calculate bin value
    ybin(i,j)=quantile(Y(filterx),(j-0.5)/NY);
    
    qyO=qy;
  end
  xbin(i)=quantile(X,(i-0.5)/NX);
  qxO=qx;
end

h=plot(xbin, ybin, varargin{:});
