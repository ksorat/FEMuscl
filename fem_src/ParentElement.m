classdef ParentElement
% Parent element superclass. This contains most of the properties that the
% subclasses (such as Q4 and H8) populate or modify for their own purposes)

    properties
        
        N; % shape functions (#int pnts x 1 ARRAY)
        NMatrix; % matrices of shape functions for mass matrix construction 
        Ns; % matrices of shape functions for nodal stress L2 projection
        StressMatrix; % "mass" matrices for nodal stress L2 projection
        dN; % shape function derivatives
        intpnts; % Gauss integration points
        Lintpnts; % Lobatto integration points for mass matrix construction
        w; % Gauss integration weights
        Lw; % Lobatto integration weights
        
    end
        
   
    methods
        
        function obj = ParentElement()

        end % Constructor
        
        function N = getN(obj)
            
            N = obj.N;
            
        end  % Get shape functions
        
        function dN = getdN(obj)
            
            dN = obj.dN;
            
        end  % Get shape function derivatives
        
        function intpnts = getIntpts(obj)
            
            intpnts = obj.intpnts;
            
        end  % Get integration points
        
        function w = getW(obj)
            
            w = obj.w;
            
        end  % Get integration weights           
        
    end
end