classdef implicitSolver
    
    properties
        
        tstep; % timestep for solver (FLOAT)
        beta = 0.25; % beta parameter
        gamma = 0.5; % gamma parameter
        
    end
    
    methods
        
        function obj = implicitSolver(tstep,beta,gamma)
            
            obj.tstep = tstep;
            obj.beta = beta;
            obj.gamma = gamma;
            
        end % constructor
        
        function obj = run(obj,endtime)
            %Inputs: endtime - FLOAT
  
            global u v a t GDOF Fext Fmax Elements M;
            
            t = 0; % start time (FLOAT)
            F = zeros(GDOF,1); % total force (GDOF x 1 VECTOR of FLOATS)
            for i = 1:length(Elements)
                
                [Elements{i},Fe] = Elements{i}.getforce; % Fe is a GDOF x 1 VECTOR of FLOATS
                F = F + Fe;
                
            end % pull initial forces from elements
            Minv = inv(M);
            a = Minv*F; % calculate initial accelerations
            uNew = u; % initialize displacments vector (a GDOF x 1 VECTOR of FLOATS)
            
            while t < endtime
                
                Fext(79) = interp1([0,endtime],[0,Fmax],t); % get rid of this
                
                while conv > criteria
                    
                    for i = 1:length(Elements)

                        [Elements{i},Fe] = Elements{i}.getforce;
                        F = F + Fe;

                    end % pull forces from elements                    
                    uPre = u + obj.tstep*v + .5*obj.tstep^2*(1 - 2*obj.beta)*a; % displacment predictor (a GDOF x 1 VECTOR of FLOATS)
                    vPre = v + (1 - obj.gamma)*obj.tstep*a; % velocity predictor step (a GDOF x 1 VECTOR of FLOATS)
                    aNext = (uNew - uPre)/(obj.beta*obj.tstep^2); % calculate acceleration (a GDOF x 1 VECTOR of FLOATS)
                    vNext = vPre + obj.gamma*obj.tstep*aNext; % velocity corrector (a GDOF x 1 VECTOR of FLOATS)
                    r = M*aNext - F; % calculate residual (a GDOF x 1 VECTOR of FLOATS)
                    
                    
                end

                
            end
            
            
        end % run the time integrator
        
        function obj = timestep(obj)
            
            global Lengths Speeds Elements;

            for i = 1:length(Elements)
                
                Elements{i} = Elements{i}.charLength;
                
            end % calculate new element characteristic lengths
            
            obj.tstep = min(Lengths)/max(Speeds);
            
        end % calculate timestep       
        
    end
    
end
    