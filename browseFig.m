function [selection]=browseFig(variable, command, valueSet, selection)
% browseFig(variable, command)
% allows you to go back and forward through figures 
% with different values of a parameter
% if valueSet is not given, values of variable start at 1 and
% range up and down the integers.
% keys: p,n

idx=1;
if(exist('valueSet')~=1)
    valueSet=[];
end
if(~exist('selection','var'))
    selection=[];
end;
helptext=['\nH:this screen, P:previous, N:next'...
    '\nS:toggle select, G:goto trial, Q/esc:exit'];

global key;
key=0; oval=0;
while key~=27
    if(length(valueSet)==0) value=idx;
    else value=valueSet(idx); end;
    if(oval~=value)
        oval=value;
        evalin('caller', [variable '=' num2str(value) ';' command ';']);
        set(gcf,'KeyPressFcn',@keypress);
        title ([variable '=' num2str(value)]);
    end
    key=0; redraw=0;
    while ~redraw
        pause(0.1);
        if exist('kbcheck','file')
          [z,z,kcode]=kbcheck;
          key=find(kcode);
        else
          [key key key]=ginput(1);
        end
        if(key=='N') idx=idx+1; redraw=1; end;
        if(key=='P') idx=idx-1; redraw=1; end;
        if(valueSet~=[] & idx<1) idx=1; redraw=1; end;
        if(valueSet~=[] & idx>length(valueSet)) idx=length(valueSet); redraw=1; end;
        if(key=='Q') key=27;end;
        if(key=='G') idx=input('Goto index:'); redraw=1;end;
        if(key=='S') 
            if(any(selection==idx)) selection(selection==idx)=[];
            else selection=[selection idx]; end
            selection
        end
        if(key=='H')fprintf(helptext);end;
        if(key==27) redraw=1;end
        if exist('kbcheck','file') % wait for key to be released
          while any(kcode) [z,z,kcode]=kbcheck; end;
        end
    end;
end

function keypress(s,k)
global key;
key=k.Character;