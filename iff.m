function RESULT = iff(CONDITION, TRUEVAL, FALSEVAL)
% RESULT = iff(CONDITION, TRUEVAL, FALSEVAL)
% implements a ternary operator in matlab
% 
% i.e. iff ( condition,  trueval,  falseval )
%        ==  condition ? trueval : falseval    [in C or Java]
%
% CONDITION can be scalar or array. 
% If it's a scalar, TRUEVAL and FALSEVAL can be anything you like.
% If it's an array, 
%   TRUEVAL must either be scalar or same size as CONDITION
%   Same goes for FALSEVAL. 
% 
% sgm

C=logical(CONDITION);

if length(C)==0
  RESULT = [];
  return;
end

if isscalar(C)
  if C, RESULT=TRUEVAL;
  else  RESULT=FALSEVAL;
  end
  return
end

if length(C)>1 % Vector condition
  if isscalar(TRUEVAL)
    RESULT(C) = TRUEVAL;
  else
    RESULT(C) = TRUEVAL(C);
  end
  if isscalar(FALSEVAL)
    RESULT(~C) = FALSEVAL;
  else
    RESULT(~C) = FALSEVAL(~C);
  end
end