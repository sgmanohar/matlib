function S = emptyStructLike(S)
% create an empty struct, with the same fields as S.
% it's what Matlab fills the struct with, if you extend it automatically.
S=S(1);
S(3)=S(1);
S=S(2);