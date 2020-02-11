function X=previous(M,n, DIM)
% X=previous( M [,n [,DIM] ] )
% return the matrix shifted to the right n spaces (default n=1),
% wrapping round across rows.
% if DIM is specified, don't wrap.
if exist('DIM','var')
  sz=size(M);
  for j=1:length(sz) % for each dim
    if(j==DIM) 
      ix{j}=1:(sz(j)-n); % indices are whole matrix except 1st plane along dim DIM
    else
      ix{j}=1:sz(j);
    end
  end 
  d(DIM)=n;d(d==0)=1;
  X=nancat(DIM, repmat(nan,d) , M(ix{:}));
else
  
  if exist('n')~=1
    X=reshape([NaN M(1:end-1)], size(M));
  else
    X=reshape([NaN*ones(1:n) M(1:end-n)], size(M));
  end;
  
end