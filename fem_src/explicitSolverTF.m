classdef explicitSolverTF
    
    properties
        
        tstepOld; % last timestep (FLOAT)
        tstep = 0; % timestep for solver (FLOAT)
        
    end
    
    methods
        
        function obj = explicitSolverTF()
            
            obj = timestep(obj);
            
        end % constructor
        
        function obj = run(obj,tsScale,endtime)
            %Inputs: tsScale - FLOAT
            %        endtime - FLOAT
            % It could be advantageous to keep this function operating on
            % global variables like it's written. The u, v, a, and Fext
            % vectors could be very long depending upon the number of DOFs
            % in the system and making copies of them could be annoying.
            % I'm sure you can use fancy stuff like pointers to do all this
            % stuff.
  
            global u v a t GDOF Fext Fmax Elements M Damage;
            
            t = 0; % start time (FLOAT)
            F = zeros(GDOF,1); % total force (GDOF x 1 VECTOR of FLOATS)
            for i = 1:length(Elements)
                
                [Elements{i},Fe] = Elements{i}.getforce; % Fe is a GDOF x 1 VECTOR of FLOATS
                F = F + Fe;
                
            end % pull initial forces from elements
            a = F./M; % calculate initial accelerations            
            
            while t < endtime
                
                obj = obj.timestep; % calculate timestep
                
                %Fext(79) = interp1([0,endtime],[0,Fmax],t); % get rid of this
                v = v + 0.5*tsScale*(obj.tstep + obj.tstepOld)*a; % velocity update
                u = u + tsScale*obj.tstep*v; % displacments update
                F = zeros(GDOF,1); % total force (a GDOF x 1 VECTOR of FLOATS)
                
                for i = 1:length(Elements)
                
                    [Elements{i},Fe] = Elements{i}.getforce;
                    F = F + Fe;
                
                end % pull forces from elements
                
                a = F./M; % calculate accelerations
                t = t + tsScale*obj.tstep; % update time
                
                if ~isempty(Damage)                
                    node.damageCheck(); % check for damage
                end
                
            end
            
            
        end % run the time integrator
        
        function obj = timestep(obj)
            
            global Lengths Speeds Elements;
            
            obj.tstepOld = obj.tstep; % log last timestep

            for i = 1:length(Elements)
                
                Elements{i} = Elements{i}.charLength;
                
            end % calculate new element characteristic lengths
            
            obj.tstep = min(Lengths)/max(Speeds);
            
            if isempty(obj.tstepOld)
                obj.tstepOld = 0;
            end
            if obj.tstepOld == 0
                obj.tstepOld = obj.tstep;
            end
            
        end % calculate timestep       
        
    end
    
end
            
            
            
            
            
            
            
            
            
            
            
            