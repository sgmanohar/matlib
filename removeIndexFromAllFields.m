function struc=removeIndexFromAllFields(struc,idx)
% each field must be an array/cell. This removes the item with index (idx) 
% from each field of the structure.
f=fieldnames(struc);
for(i=1:length(f))
    struc.(f{i})(idx)=[];
end

