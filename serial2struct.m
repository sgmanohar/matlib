function structure = serial2struct(serialObj)
% function structure = serial2struct(serialObj)
% convert from a serial into a structure
% JPG.

if ~isa(serialObj, 'serial')
    error('input is not a serial');
end

% get names in serialObj
fn = fieldnames(serialObj);

for i = 1:length(fn)
    
    % store each into a new structure
    structure.(fn{i}) = serialObj.(fn{i});
    
end

end