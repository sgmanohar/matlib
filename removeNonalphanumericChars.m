function str=removeNonalphabeticChars(str)
% Remove all characters that are not alphanumeric, . or _,  and remove 
% any initial digits.
% useful for making field names or variable names.
    i=1;
    while i<=length(str)
        if any(str(i)=='_.') i=i+1;continue; end;
        if str(i)<65 | (str(i)>90 & str(i)<97) | str(i)>122 ...
                | (i==1 & str(i)>47 & str(i)<58)
            str=[ str(1:i-1)  str(i+1:end) ];
        else i=i+1;
        end
    end;
    