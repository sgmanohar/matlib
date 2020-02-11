function x=currentWorkspace()
% return a struct containing the variables of the current workspace
% - a simple routine that should be part of the API

z=evalin('caller','who');
for(i=1:length(z))
    x.(z{i})=evalin('caller',z{i});
end;
