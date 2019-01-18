function [cope,varcope,tstat,res, covb]=ols(data,des,tc,varargin)
% [COPE,VARCOPE TSTAT RES]=OLS(DATA,DES,TC)
% DATA IS T x V
% DES IS T x EV (design matrix)
% TC IS NCONTRASTS x EV  (contrast matrix)
%
% returns COPE = contrasts of parameter estimates
%      VARCOPE = variance of each contrast
%        TSTAT = t statistic for each contrast
%          RES = residuals (T x V) 
% TB 2004
% modified SGM 2012-2018

if ndims(data)>2
  warning('operating on dimension 1 of %g',ndims(data));
  dataSize = size(data); 
  data = reshape(data, size(data,1),[]);
  RESHAPE_AFTER =true;
else RESHAPE_AFTER = false;
end

if ~exist('des','var')
  warning('des = [ [1:%g]'', ones(%g,1) ]',size(data,1),size(data,1));
  des=[zscore([1:size(data,1)]'), ones(size(data,1),1)];
end

if ~exist('tc','var'), tc=eye(size(des,2)); end

if(size(data,1)~=size(des,1))
  error('OLS::DATA and DES have different number of time points');
elseif(size(des,2)~=size(tc,2))
  error('OLS:: DES and EV have different number of evs')
end




    
%%%%% ORIGINAL CODE:

pdes=pinv(des);
pe=pdes*data;
cope=tc*pe; 

%%%%% 


if RESHAPE_AFTER
  newsize = dataSize; newsize(1) = size(cope,1); 
  newsize = num2cell(newsize);
  cope = reshape(cope,[], newsize{:}); 
end

if(nargout>1)
  prevar=diag(tc*pdes*pdes'*tc');
  R=eye(size(des,1))-des*pdes;
  tR=trace(R);
  res=data-des*pe;
  sigsq=sum(res.*res/tR);       
  varcope=prevar*sigsq;
  if(nargout>2)
    tstat=cope./sqrt(varcope);
  end
end

if nargout>4
  covb=(tc*des'*des*tc')*sigsq;
end