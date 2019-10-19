
function sc=isFixated(p, thresholdv)
% quick helper function to check if a set of coordinates counts as
% "fixated". called by WaitForFixation.m.
% p: coordinates of eye ( time , xy )
% sgm 2011
    dt=(p(2:end,1) - p(1:end-1,1));
    vx=(p(2:end,2) - p(1:end-1,2)) ./ dt;
    vy=(p(2:end,3) - p(1:end-1,3)) ./ dt;
    v2=vx.*vx + vy.*vy; %individual v-squared
    vv=((p(end,2)-p(1,2))^2+(p(end,3)-p(1,3))^2) / (p(end,1)-p(1,1))^2; %nett v-squared
    sc=vv<thresholdv & all(v2<2*thresholdv);
