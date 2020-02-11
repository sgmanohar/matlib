function p=my_recinormpdf(x,m,s)
% p=my_recinormpdf(x,m,s)
%    reci-normal probability density function
%
% probability density = my_recinormpdf( RT, MU, SIGMA )
%  
%                                2
%                       _  (mx-1)
%                         --------
%         1                   2  2
%  --------------          2 s  x
%     2            *  e
%  s x  sqrt(2pi) 
% 
%
%   sgm

 p = 1 / sqrt(2*pi*s^2) ./ (x.^2) .* ...
     exp( -((m*x-1).^2) / ...
         2/(s^2) ./ (x.^2) );
 