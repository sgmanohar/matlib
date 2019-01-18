function s=medianSplit(X)
% median split that gracefully deals with badly behaved samples
% sgm 2015
m = nanmedian(X);
s = X>m;
if mean(s) < 0.25 && mean(X>=m)>mean(s) && (mean(X>=m) < 1-mean(s))
  s=X>=m;
end
if mean(s) < 0.25 || mean(s) > 0.75
  warning('unable to evenly split data; %g of sample is extreme', mean(s));
end


