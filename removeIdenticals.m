function V=removeIdenticals(V)
% function V2=removeIdenticals(V)
% vector V - replaces identicals with forward-interpolated versions
rep = [0; diff(V)==0];
N=length(V); 
iter=0;
while any(rep)
  i=find(rep,1); % find first repeated value
  j=i;
  while rep(j) && j<N % find end of repeats (j)
    j=j+1;
  end
  if rep(j) % reached the end
    if i>2
      V(j)=V(j)+(V(j)-V(i-2))/2; % add on a bit
    elseif i>1
      V(j)=V(j)+(V(j)-V(i-1))/2;
    else
      V(j)=V(j)+eps; %% problem
    end
  end
  % interpolate values
  V(i:j-1) = interp1([i-1 j], [V(i-1) V(j)], i:j-1); 
  rep = [0; diff(V)==0];
  iter=iter+1;
  if(iter>1000)
    warning('unable to remove identicals in (%g-%g)=(%g-%g)',i,j,V(i),V(j));
    return
  end
end

