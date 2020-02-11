function [d]=makedstruct()
% create a structure with all the current workspace variables in it.
vars = evalin('base', 'who');
N=prod(size(vars));
for x=1:N
    if(vars{x}=='d') 
        error('d already exists');
    end
end
data=cell(1,N);
for x=1:N
    data{x}=evalin('base', vars{x});
end
d=cell2struct(data,vars,2);