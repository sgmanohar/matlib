
function sc=isSaccade(p, thresholdv)
% helper function to calculate if a sequence of eye
% samples should count as a saccade.
% Useful for online tracking.
% sgm 2011
    dt=(p(2:end,1) - p(1:end-1,1));
    vx=(p(2:end,2) - p(1:end-1,2)) ./ dt;
    vy=(p(2:end,3) - p(1:end-1,3)) ./ dt;
    v2=vx.*vx + vy.*vy; %individual v-squared
    vv=((p(end,2)-p(1,2))^2+(p(end,3)-p(1,3))^2) / (p(end,1)-p(1,1))^2; %nett v-squared
    sc=vv>2*thresholdv & all(v2>thresholdv);
