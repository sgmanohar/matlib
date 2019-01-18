function O = getOutput(N, func, varargin)
% O = getOutput( N, function, varargin )
% calls function( varargin{:} ) and returns the Nth output of the function
% sgm 2014
O = repmat({[]}, 1, N); % create a cell array of N empty arrays, to hold the outputs
[O{:}] = func(varargin{:}); % capture the first N outputs
O = O{N}; % return just the Nth output