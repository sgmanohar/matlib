function [e]=entropy(x, limit, circular)
% [E] = entropy (x, limits)
%
% estimate entropy of probability distribution of x
% calculated as sum(p log(p)) between the limits
% assumes that p(limit(0)) = p(limit(1)) and p is zero outside.
%
% the result contains [entropy, fisher_information]
% 
% sanjay manohar

if(length(limit)~=2) error('limit must have 2 elements'); end;
if ~exist('circular','var')
    circular=1;
end

% use binary logarithms? 1=true
BITS = 1;


% problem distributions: try shifting the distribution - 
% helps problematic values by getting rid of singularities...
% not sure why, but it helps!
limit=sort(limit);
if prod(limit)<=0 % limit crosses zero
    x     = x - limit(1) + 1;
    limit = limit - limit(1) + 1;
end


% calculate probability density function
delta=diff(limit)/1000;
i=[limit(1):delta:limit(2)];
p=ksdensity(x,i,'support',limit);
if(circular)
    p=p(2:end-1); % truncate endpoints
end

% Calculate integral p log p
ent=sum(delta*p.*log(p+eps));
% Try with fisher information = E [ d/dS( d/dS( log p ) ) ]
logp = conv(log(p+eps), ones(50,1)/50, 'full'); % smoothing
logp = logp(60:end-60);                         % truncate 10 points off each end
if(BITS)
    logp=logp/log(2);
end
gr   = diff(logp);                              %/delta ?
fisher=sum(p(36:end-36).*gr.*gr);


% convert to bits if needed?
if(BITS)
    ent=ent/log(2);
end
e = [ent fisher];
return










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% OLD VERSION!!! %%%%%%%%%%%
%%%% BEFORE I DISCOVERED KSDENSITY %%%%
% i should stop reinventing the wheel %


N=length(x);
x=sort(x);
if(size(x,1)>1)
    if(size(x,2)>1) warning('data matrix will be linearised'); x=x(:);end;
    x=x';
end;
if(x(1)<limit(1))   warning('lower limit breached'); limit(1)=x(1)-1;  end;
if(x(end)>limit(2)) warning('upper limit breached'); limit(2)=x(end)+1;end;
midpoints=(x(1:end-1)+x(2:end))/2;
dx=diff([limit(1) midpoints limit(2)]);


% repeated values give singularities
singularities=find(dx==0);
while length(singularities)>0
    i=singularities(1);
    j=1; 
    while (j+1<=length(singularities) & singularities(j+1)==singularities(1))
        j=j+1;
    end
    % divide delta of previous point
    dx(i-1)=dx(i-1)/numsing;
    % delete zeros
    dx(i:i+j-1)=[];
    singularities=find(dx==0);
end


% probabilities at the midpoints between measurements
%  estimate probability of each point from intervals
%  as 1/dx as a proportion of what is expected with a flat distr.
p=(limit(2)-limit(1)) / ((N-1)) ./ diff(x);

% what happens at the extremes?
%  take weighted mean between p(1) and p(end), depending on distance of 
%  most extreme data point and the limits.
dxleft=x(1)-limit(1);
dxright=limit(2)-x(end);
if circular
    limp=p(end)+(p(1)-p(end))*(dxright)/(dxright+dxleft); 
else
    limp=eps;
end


xx =[limit(1) x limit(2)];
dxx=diff(xx);
pp =[limp p limp];
aucP = sum((pp(1:end-1)+pp(2:end))/2.*dx);
pp=pp/aucP;
dpdx0 = diff(pp)./diff(xx(1:end-1));
dpdx = diff(pp)./dx;
ddxx=(xx(3:end)-xx(1:end-2))/2;

% Integral[ p log p ] for p = ax+b     approximates to
% [ 2(ax+b)^2 Log(ax+b) - ax(2b+ax) ] /4a

hfn=@(x,a,b)(2*(a*x+b)^2*log(a*x+b)-a*x*(a*x+2*b))/(4*a);
for i=1:length(pp)-1
    h(i) = hfn(dx(i), dpdx(i), pp(i)) - hfn(0, dpdx(i), pp(i));
end;
h2= p.*log(p).*ddxx(1:end-1);
I = (pp(1:end-1)+pp(2:end)).*(dxx(2:end))/2/(limit(2)-limit(1));
h3=dx.* (  pp(2:end).^2 ...
          - pp(1:end-1).^2 ...
          + 2*(pp(1:end-1).^2).*log(pp(1:end-1)) ...
          - 2*(pp(2:end)  .^2).*log(pp(2:end)) ...
         )./(4*(pp(1:end-1)-pp(2:end)));
h4=dx./(4*(pp(2:end)-pp(1:end-1))).* ...
     (   pp(1:end-1).^2.*(1-2*log(pp(1:end-1))) ...
       - pp(2:end)  .^2.*(1-2*log(pp(2:end)  ))  );
  
e= [ks_e sum(h) sum(h2) sum(h3) sum(h4)];
