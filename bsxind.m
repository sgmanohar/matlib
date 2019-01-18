function Y = bsxind(X,I)
% Y = bsxindex ( X , I )
% 
% For singleton indices, I is expanded to the same size as X
% then the operation X(I) is performed.
% 
% I must be a scalar / vector / matrix of integers or logicals.
% Also, size(I) must match size(X) for all non-singleton values of size(I) 
% 


% perform singleton expansion using BSXFUN
if isnumeric(I)  
  Ifull = bsxfun(@times, I, ones(size(X)) );
  if all(Ifull(:)==0 | Ifull(:)==1)
    Ifull=logical(Ifull);
    warning('bsxfun','Assuming logical subscripting');
  end
  Xfull=X;
elseif islogical(I)
  Ifull = bsxfun(@and,  I, true(size(X))  );
  Xfull = bsxfun(@plus, X, zeros(size(I)) );
else
  error('bsxindex','Unrecognised index type, %s', class(x));
end
% and apply the indexing.
Y = Xfull(Ifull);
Y=reshape(Y,size(X));
