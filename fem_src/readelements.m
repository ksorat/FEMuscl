function [ElList,numels] = readelements(filename,elementline,elementend)

fid = fopen(filename);
numels = 0;

for i = 1:elementline % skip to element section
    
    s = fgetl(fid);
    
end
    
for i = elementline + 1 : elementend
    
    s = fgetl(fid);
    comtest = strfind(s,'$'); % Test for comment lines
    
    if isempty(comtest) 
        
        numels = numels + 1;
        ElList(numels,:) = str2num(s); %#ok<*AGROW,*ST2NM>
        
    end
    
end

global Elements M GDOF;

M = zeros(GDOF);
Elements = cell(length(ElList),1);

for i = 1:length(Elements)
    
    temp = ElList(i,3:end);
    temp = temp(temp>0);
    
    Elements{i} = element(i,temp,ElList(i,2));
    
end
    

fclose(fid);