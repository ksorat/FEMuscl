classdef MaterialModel
% Superclass for material models. Contains the constitutive C matrix, 
% material properties, and the internal force solver for the specific
% material model.    
    
    properties
        
        ID; % Material model ID (INT)
        C; % stress-strain relation matrix (3 x 3 or 6 x 6 MATRIX of FLOATS)
        rho; % mass density (FLOAT)
        c; % sound speed (FLOAT)
        
    end
    
    methods
        
        function obj = MaterialModel(ID,rho)
            % Inputs: ID - INT
            %         rho - FLOAT

            obj.rho = rho; % mass density
            obj.ID = ID; % material model ID        
            
        end % Constructor
        
           
        
    end
    
end
            
            
            
            
            
        

        