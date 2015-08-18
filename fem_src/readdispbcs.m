function DispBC = readdispbcs(filename,dispbcline,dispbcend)

global Nodes;
fid = fopen(filename);
count = 0;

for i = 1:dispbcline % skip to element section
    
    s = fgetl(fid);
    
end
    
for i = dispbcline + 1 : dispbcend
    
    s = fgetl(fid);
    comtest = strfind(s,'$'); % Test for comment lines
    
    if isempty(comtest) 
        
        count = count + 1;
        DispBC(count,:) = str2num(s); %#ok<*AGROW,*ST2NM>
        
        if DispBC(count,3)
            Nodes{DispBC(count,1)}=Nodes{DispBC(count,1)}.setSPCs(1);
        end
        if DispBC(count,4)
            Nodes{DispBC(count,1)}=Nodes{DispBC(count,1)}.setSPCs(2);
        end
        if DispBC(count,5)
            Nodes{DispBC(count,1)}=Nodes{DispBC(count,1)}.setSPCs(3);
        end
        
    end
    
    
    
end

for i = 1:length(Nodes)
    
    Nodes{i} = Nodes{i}.setDOFs;
    
end


fclose(fid);