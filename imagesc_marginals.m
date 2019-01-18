function [A H X2]=imagesc_marginals(X, varargin)
% function A=imagesc_marginals(X)
% calls imagesc(X)
% and also plots the marginal means of the image in panels at the edges.
% 
% arguments: 
%  'Ratio' = fraction of current axes that are used for the 
%            main image (the remainder is used for the marginal plots)
%            default 0.8
%  'RemoveMarginals' = before plotting X, subtract out the marginal means
%            to give the deviation of the joint distribution from the 
%            product of the marginals. This is very useful for revealing 
%            a correlations  
% returns :
%  H = a list of 3 axis handles
%      [ main_axis, right_axis, top_axis ]
%  A = a list of 3 graphics objects handles
%      [ image_handle, right_line, top_line ]
%  X2 = the actual values plotted (i.e. not X if RemoveMarginals is used)

a=gca; r=get(a,'position'); % left, bottom, width, height
delete(a);
F=0.8;
REMOVE_MARGINALS=false;
i=find(strcmpi(varargin,'ratio')); if i, F=varargin{i+1}; varargin([i i+1])=[]; end
i=find(strcmpi(varargin,'removemarginals')); if i, REMOVE_MARGINALS=varargin{i+1}; varargin([i i+1])=[]; end

if 0 % NORTHWEST?
  a1=axes('position',[r(1), r(2)+r(4)*(1-F), r(3)*F, r(4)*F]); % main
  a2=axes('position',[r(1)+r(3)*F, r(2)+r(4)*(1-F), r(3)*(1-F), r(4)*F]); % right
  a3=axes('position',[r(1), r(2), r(3)*F, r(4)*(1-F)]); % bottom
else
  a1=axes('position',[r(1), r(2), r(3)*F, r(4)*F]); % main
  a2=axes('position',[r(1)+r(3)*F, r(2), r(3)*(1-F), r(4)*F]); %right
  a3=axes('position',[r(1), r(2)+r(4)*F, r(3)*F, r(4)*(1-F)]); % top
end
A=[a1 a2 a3];

if REMOVE_MARGINALS
  X2 = bsxfun(@times,mean(X,1),mean(X,2)) ; % product of marginals
  X2 = X - X2*sum(X(:))/sum(X2(:)); % deviation of joint 
else X2=X; end

axes(a1);
h1=imagesc(X2,varargin{:});
axes(a2);
h2=plot(sum(X,2), fliplr([1:size(X,1)]) ); % flip because image has 1 at top!
set(a2,'ytick',[]);
ylim([0.5 size(X,1)+0.5]);
axes(a3);
h3=plot([1:size(X,2)], sum(X,1));
set(a3,'xtick',[]);
xlim([0.5 size(X,2)+0.5]);
axes(a1);

H=[h1 h2 h3];

