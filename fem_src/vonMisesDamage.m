classdef vonMisesDamage
    
    properties
        
        threshold; % value for damage threshold
                
    end
    
    methods
        
        function obj = vonMisesDamage(threshold)
            
            obj.threshold = threshold;
            
        end
        
        function D = damageCalc(obj,sigma)
            
            global dim           
            
            if dim == 2                
                vm = sqrt(sigma(1)^2 + sigma(2)^2 - sigma(1)*sigma(2) + 3*sigma(3)^2); % von Mises equivalent stress
            elseif dim == 3
                vm = sqrt(.5*((sigma(1) - sigma(2))^2 + (sigma(2) - sigma(3))^2 + (sigma(3) - sigma(1))^2 + 6*(sigma(6)^2 + sigma(4)^2 + sigma(5)^2))); % von Mises equivalent stress
            end
            
            D = vm/obj.threshold; % calculate damage value
            
        end            
        
    end
    
end