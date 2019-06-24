function p = inversePermutation(from,to)
% p = inversePermutation(from,to)
%  return the indices of the items of 'from', in 'to'.
%
%  this is useful if you have used a permutation to get from one list to a
%  new ordering of that list, and you want to know what that permutation
%  was.
%  
%  in other words, this returns the vector p, such that
%  
%    to(p) == from
% 
% 'from' and 'to' must be vectors.
% 
% basically, if the inputs are 
%  from   [1 2 3 4 5],
%  to     [1 5 4 2 3]
% then the result p is
%         [1 4 5 3 2]
%
%  warnings will occur if 'from' contains items that are not in 'to',
%  or if items in 'from' are contained multiple times in 'to'
%
% sgm



if(length(from)~=length(to)) 
  error('''from'' and ''to'' should be same length');
end;
fe=0; te=0;
for(i=1:length(from))
  q=find(to==from(i));
  if(isempty(q)) p(i)=nan; te=1;
  else p(i)=q(1);
    if(length(q)>1) fe=1;end;
  end;
end;
if(fe) 
  warning('permute:duplication','some items were duplicated in ''to''');
end
if(te)
  warning('permute:duplication','some items were not found in ''to''');
end