function j = jet2(N)
% jet2 palette - sgm 2019
% removes the dark bits at the top and bottom of the palette
% to make it brighter
if ~exist('N','var'), N=size(colormap,1); end
front = floor(0.24*N);  % discard the dark blue
back  = floor(0.09*N);  % discard the dark red
j=jet(N+front+back);    % get original jet, with larger size than needed
j=j(front+1:end-back,:);% cut out unwanted bits

