classdef Part
    
    properties
        
        elements; % element list (ARRAY of INTs or LONGs full of indices to the global Element Array
        material; % material model object
        elForm; % parent element object
        damage; % damage object reference
        
    end
    
    methods
        
        function obj = Part(material,elForm)
            
            obj.material = material;
            obj.elForm = elForm;
            
        end
        
        
        
    end
    
end
        