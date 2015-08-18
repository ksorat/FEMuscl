classdef Q4 < ParentElement
% Create a Q4 parent element object.

    properties
        
        connections = [2 4;...
                       1 3;...
                       2 4;...
                       1 3]; % array mapping connected nodes (4 x 2 ARRAY of INTs)
                   
    end

    methods
        
        function obj = Q4()
            
            obj@ParentElement(); % identify with the parent element superclass
            tmp = sqrt(1/3); % integration point temp constant (FLOAT)
            obj.intpnts = [-tmp,-tmp;...
                            tmp,-tmp;...
                            tmp, tmp;...
                           -tmp, tmp]; % Gauss integration points (4 x 2 ARRAY of FLOATS)
            obj.Lintpnts = [-1,-1;...
                             1,-1;...
                             1, 1;...
                            -1, 1]; % Lobatto integration points (4 x 2 ARRAY of FLOATS)
            obj.w = ones(4,1); % Gauss integration weights (4 x 1 ARRAY of FLOATS)
            obj.Lw = ones(4,1); % Lobatto integration weights (4 x 1 ARRAY of FLOATS)
            obj = obj.shapeFunctions; % call the shapeFunctions function to create the shape functions
            
        end % Constructor
        
        function obj = shapeFunctions(obj)
            
            obj.N = cell(4,1); % initialize shap function cell (4 x 1 ARRAY of 4 x 1 ARRAYS of FLOATS)
            Nm = cell(4,1); % intialize shape function cell for mass matrix (4 x 1 ARRAY of 4 x 1 ARRAYS of FLOATS)
            obj.dN = cell(4,1); % initalize shape function derivatives cell (4 x 1 ARRAY of 2 x 4 ARRAYS of FLOATS)            
            obj.NMatrix = cell(4,1); % initialize shape function cell for mass matrix creation (4 x 1 ARRAY of 2 x 8 MATRICES of FLOATS)
            obj.StressMatrix = cell(4,1); % initialize shape function cell for stress matrix creation (4 x 1 ARRAY of 12 x 12 MATRICES of FLOATS)
            
            for i = 1:4
            
                obj.N{i} = .25*[(1-obj.intpnts(i,1))*(1-obj.intpnts(i,2)),(1+obj.intpnts(i,1))*(1-obj.intpnts(i,2)),(1+obj.intpnts(i,1))*(1+obj.intpnts(i,2)), ... % Shape functions for bilinear quads
                        (1-obj.intpnts(i,1))*(1+obj.intpnts(i,2))];
                
                % Create "mass" matrix for L2 stress projection to nodes    
                obj.Ns{i} = zeros(3,12);    
                    
                obj.Ns{i}(:,1:3) = obj.N{i}(1)*eye(3);
                obj.Ns{i}(:,4:6) = obj.N{i}(2)*eye(3);
                obj.Ns{i}(:,7:9) = obj.N{i}(3)*eye(3);
                obj.Ns{i}(:,10:12) = obj.N{i}(4)*eye(3);
                
                obj.StressMatrix{i} = obj.Ns{i}'*obj.Ns{i}*obj.w(i);
                obj.Ns{i} = obj.Ns{i}*obj.w(i);
                
                % Create mass matrix for dynamic solver
                Nm{i} = .25*[(1-obj.Lintpnts(i,1))*(1-obj.Lintpnts(i,2)),(1+obj.Lintpnts(i,1))*(1-obj.Lintpnts(i,2)),(1+obj.Lintpnts(i,1))*(1+obj.Lintpnts(i,2)), ... % Shape functions for bilinear quads
                        (1-obj.Lintpnts(i,1))*(1+obj.Lintpnts(i,2))];                    

                obj.NMatrix{i} = zeros(2,8);    
                    
                obj.NMatrix{i}(:,1:2) = Nm{i}(1)*eye(2);
                obj.NMatrix{i}(:,3:4) = Nm{i}(2)*eye(2);
                obj.NMatrix{i}(:,5:6) = Nm{i}(3)*eye(2);
                obj.NMatrix{i}(:,7:8) = Nm{i}(4)*eye(2);                    

                obj.dN{i} = .25*[-(1-obj.intpnts(i,2)), (1-obj.intpnts(i,2)), (1+obj.intpnts(i,2)), -(1+obj.intpnts(i,2)); ... % Derivative wrt xi
                        -(1-obj.intpnts(i,1)), -(1+obj.intpnts(i,1)),  (1+obj.intpnts(i,1)),   (1-obj.intpnts(i,1))]; % Derivative wrt eta
                    
            end
            
        end % Solve shape functions at integration points for each node        

        
    end
    
end