function txt = as_text(var)
%   txt = as_text(var)
% Convert a variable to a displayed string.
% how it works:
%  - Display the variable as if matlab had output it using 'disp',
%  - Capture the text in a temp file, and return it in a string
% sgm 2020

d0 = get(0,'Diary');     % get diary state
f0 = get(0,'DiaryFile'); % get diary file
f1 = [tempname '.txt'];  % get a temporary file
diary(f1);               % start diary to that file
disp(var)      % output the variable to the screen, whcih should now be logged to diary
if strcmp(d0,'off')      % restore the original state of the diary
  diary('off');
else
  diary(f0);    
end
txt = fileread(f1);   % grab the text into a variable
r0  = recycle('off'); % turn off recycler
delete(f1);           % delete temporary file
recycle(r0);          % set recycler back to original state


