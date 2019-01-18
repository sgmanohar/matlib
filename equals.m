function r=equals(a,b)
% returns whether a and b are identical
% i.e. false if they are of different types or sizes.
% true if they are identical. structs must have identical
% fieldnames and contents but order doesnt't matter.
%
% e.g. a=[], b=[] returns true
%      a=struct('x',1,'y',2), b=struct('y',2,'x',1) returns true
%      
% sanjay manohar 2010
if(exist('b')~=1) error('syntax: equals(a,b)');end;
r=0;
if(class(a)~=class(b)) return;end;
if(isnumeric(a))
        if(size(size(a))~=size(size(b))) return;end;
        if(size(a)~=size(b)) ;return;end;
        if(length(a)==0) r=1;return ;end;
        if any(a~=b) ;return;end;
        r=1;return;
end;
if(iscell(a))
    if(size(size(a))~=size(size(b))) return;end;
    if(size(a)~=size(b)) return;end;
    for i=1:prod(size(a))
        if(~equals(a{i},b{i})) return;end;
    end;
    r=1;return;
end;
if(isstruct(a))
    fa=fieldnames(a);
    if(size(fa)~=size(fieldnames(b))) return;end;
    for(i=1:length(fa))
        if(~isfield(b,fa{i})) return;end;
        if(~equals(a.(fa{i}),b.(fa{i}))) return;end;
    end;
    r=1;return;
end;
if(islogical(a))
    r=all(a==b);
    return;
end;
if(ischar(a))
    r=strcmp(a,b);
    return;
end;
