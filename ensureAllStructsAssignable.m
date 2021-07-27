function new = ensureAllStructsAssignable(old, DEBUG)
% checks each structure in an array to ensure they all have the same
% fields, will also reorder all the fields into the same order
% new = ensureAllStructsAssignable(old, DEBUG)
% old is a structure array, DEBUG is [0 1] on whether to print messages
% when mismatches are found (default is 1)
%% check inputs

if ~iscell(old)
    error('old should be a cell array of structures');
end

if ~exist('DEBUG','var') || isempty(DEBUG)
    DEBUG=true;
end
%% loop through each struct, updating all previous if one changes

n = numel(old);
new = old; % copy

i = 2;
while i <= n
    
    [new{i}, new1] = ensureStructsAssignable(new{i}, new{1}, DEBUG);
    
    if ~equals(new1, new{1}) % if new has been updated
        new{1} = new1; % store it
        % and redo all others
        i = 2;
    else
        i = i + 1;
    end
end



