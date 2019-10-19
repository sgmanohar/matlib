function [extra1, extra2] = compareLists(l1, l2)
% [extra1, extra2] = compareLsts(l1, l2)
% compare two cell arrays
% the terms in l1 that aren't in l2 are given in extra1

extra1={};
extra2={};
VERBOSE = true;

for i=1:length(l1)
  c = l1{i};
  found=false;
  for j=1:length(l2)
    if strcmp(c,l2{j})
      found=true;
      break
    end
  end
  if ~found
    extra1=[extra1; {c}];
    if VERBOSE
      disp(c)
    end
  end
end
      
if VERBOSE, fprintf('------\n'); end

for i=1:length(l2)
  c = l2{i};
  found=false;
  for j=1:length(l1)
    if strcmp(c,l1{j})
      found=true;
      break
    end
  end
  if ~found
    extra2=[extra2; {c}];
    if VERBOSE
      disp(c)
    end
  end
end
