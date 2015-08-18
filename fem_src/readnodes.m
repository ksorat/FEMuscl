function [NodeList,numnodes] = readnodes(filename,nodeline,nodeend)

fid = fopen(filename);
numnodes = 0;

for i = 1:nodeline % skip to node section
    
    s = fgetl(fid);
    
end
    
for i = nodeline + 1 : nodeend
    
    s = fgetl(fid);
    comtest = strfind(s,'$'); % Test for comment lines
    
    if isempty(comtest) 
        
        numnodes = numnodes + 1;
        NodeList(numnodes,:) = str2num(s); %#ok<*AGROW,*ST2NM>
        
    end
    
end

NodeList = NodeList(:,1:4);

global Nodes dim;

Nodes = cell(length(NodeList),1);

for i = 1:length(Nodes)
    
    if dim == 2
        
        Nodes{i} = node(i,NodeList(i,2),NodeList(i,3));
        
    else
        
        Nodes{i} = node(i,NodeList(i,2),NodeList(i,3),NodeList(i,4));
        
    end
    
end

fclose(fid);