function [dispbcline,elementline,nodeline,dispbcend,elementend,nodeend] = LineFind(filename)
fid = fopen(filename);
endlines = zeros(1,3);
endflags = zeros(1,3);
starflag = 0;
count = 1;
while ~feof(fid)
    
    s = fgetl(fid); % load each line
    
    if ischar(s) % if it's a character
        
        if starflag % if the star flag is on
        
            star = strfind(s,'*'); % search for a *
            
            if star >= 1 % if there is a star
                
                endlines(find(endflags)) = count - 1; % mark the end of corresponding section
                endflags(find(endflags)) = 0; % reset the endflag marker
                starflag = 0; % reset the starflag marker
                
            end
            
        end
                
        
        switch s

            case '*BOUNDARY_SPC_NODE'
                dispbcline = count;
                endflags(1) = 1;
                starflag = 1;
            case '*ELEMENT_SHELL'
                elementline = count;
                endflags(2) = 1;
                starflag = 1;
            case '*NODE'
                nodeline = count;
                endflags(3) = 1;
                starflag = 1;
        end
        
    end
    
    count = count + 1;
    
end

dispbcend = endlines(1);
elementend = endlines(2);
nodeend = endlines(3);

fclose(fid);
end
            
    