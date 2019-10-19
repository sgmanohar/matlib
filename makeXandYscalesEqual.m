function makeXandYscalesEqual()
% makeXandYscalesEqual()
% makes the X and Y axes the same limits (square plot)
xl=xlim(); yl=ylim();
lo=min(xl(1),yl(1));
hi=max(xl(2),yl(2));
xlim([lo hi]); ylim([lo hi]);