function S = fieldfun(fun, S, varargin)
% S = fieldfun(fun, S)
% apply a function to each field of structure S
% (and also to each element of the array if S is a struct-array)
% additional arguments are passed to fun.

% for each single structure-array element, apply the function to it
S = arrayfun( @(s)fieldfun_single(fun,s, varargin{:}),   S ); 



function S = fieldfun_single(fun,S, varargin)
% apply function to each field of a single structure S
fn = fieldnames(S); 
for i=1:length(fn) % for each field
  S.(fn{i}) = fun(S.(fn{i}), varargin{:}); 
end




