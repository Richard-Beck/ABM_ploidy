function output=readCurrentPop(filename)
pop=load(filename);
index= pop(:,1) ~= 0 | pop(:,2) ~= 0; % Extract all those locations with cells
pop=pop(index,1:2);

output=pop;

end