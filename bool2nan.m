function Y=bool2nan(X)
% Y = bool2nan(X)
% convenience method to convert a matrix of booleans to a matrix of 0 or nan
% just replaces all the 'trues' with 'nans'.
%   i.e.  X(X) = nan
%
% example of how it is useful: 
%    X = rand(10)
%    nanmean( X .* bool2nan(X<0.5), 1)
% Takes the means of the columns, ignoring items that are less than 0.5.

Y=X*1; % convert to numeeric first
Y(X>0)=nan;