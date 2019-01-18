function struct2workspace(s, CHECK)
% STRUCT2WORKSPACE( struct )
% STRUCT2WORKSPACE( struct, CHECK )
%
%  Dump the fields of a structure into caller's workspace as variables.
% 
% CHECK is true by default, and protects vars in current workspace from
% being overwritten.
% WARNING 
% if 'check' is false, then don't check for overwriting.
%   - this will overwrite any variables in your current workspace 
%   that correspond to fields in the structure.
%
% Sanjay Manohar 2014


if ~exist('CHECK','var'), CHECK=true;  end;

f=fieldnames(s); for i=1:length(f)
  if CHECK
    try
      tmp=evalin('caller',f{i});
      warning('workspace:clash',sprintf('%s not assigned',f{i}));
    catch me
      assignin('caller',f{i}, s.(f{i}) );
    end
  else
    assignin('caller',f{i}, s.(f{i}) );
  end
end

