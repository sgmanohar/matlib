function x=workspace2struct(exclude)
% X = WORKSPACE2STRUCT( [EXCLUDE] )
% 
%  Return a struct containing the variables of the current workspace
%  Wach element of the struct X is a current workspace variable.
%
%  note that a similar effect could be achieved using 
%     save('tmp'); X=load('tmp');
% 
%  to copy struct fields back to workspace, use struct2workspace.
%  if EXCLUDE is specified, it is a regexp string for variable names to be
%  excluded from the structure. 
%
% Sanjay Manohar 2014

z=evalin('caller','who');
for(i=1:length(z))
  if ~exist('exclude','var') || ~any(regexp(z{i},exclude))
    x.(z{i})=evalin('caller',z{i});
  end
end;
