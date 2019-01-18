function E=enum(names)
% E = enum(names)
% an enumeration is a list of possible categorical values, represented by
% strings.
% 
% names: a cell array of strings, containing the names of the constant 
%        variables in the enum.
% E:     returns a structure containing fields with the names specified, 
%        taking unique integer values.
%        E.names stores the original list, so you can easily look up the 
%        name for a given integer.

for i=1:length(names)
  E.(names{i}) = i;
end
E.names  = names;
E.values = 1:length(names);

