classdef BMatrix2D
    
    properties
        
        numNodes; % number of nodes (INT)
        
    end
    
    methods
        
        function obj = BMatrix2D(numNodes)
            % Inputs: numNodes - INT
            
            obj.numNodes = numNodes;
                       
        end % Constructor
        
        function [Bv,B,J] = BMatrix(obj,xyz,parent,num)
            % Inputs: xyz - 3 x 1 ARRAY of FLOATS
            %         parent - parent element object
            %         num - current integration point (INT)
            
            Na_x = zeros(1,obj.numNodes); % 1 x numNodes ARRAY of FLOATS
            Na_y = zeros(1,obj.numNodes); % 1 x numNodes ARRAY of FLOATS
            Bv = zeros(3,2*obj.numNodes);% 3 x 2*numNodes MATRIX of FLOATS
            B = zeros(2,obj.numNodes); % 2 x numNodes MATRIX of FLOATS
            x_xi = 0; % FLOAT
            x_eta = 0; % FLOAT
            y_xi = 0; % FLOAT
            y_eta = 0; % FLOAT           
            
            for i = 1:obj.numNodes
                
                x_xi = x_xi + parent.getdN{num}(1,i)*xyz(i,1); % derivative of x wrt xi
                x_eta = x_eta + parent.getdN{num}(2,i)*xyz(i,1); % derivative of x wrt eta
                y_xi = y_xi + parent.getdN{num}(1,i)*xyz(i,2); % derivative of y wrt xi
                y_eta = y_eta + parent.getdN{num}(2,i)*xyz(i,2); % derivative of y wrt eta
                
            end % components of element to parent jacobian matrix
            
            J = x_xi*y_eta - x_eta*y_xi; % jacobian determinent (FLOAT)

            for i = 1:obj.numNodes
                
                Na_x(i) = (parent.getdN{num}(1,i)*y_eta - parent.getdN{num}(2,i)*y_xi)/J; % derivative of shape function wrt x
                Na_y(i) = -(parent.getdN{num}(1,i)*x_eta - parent.getdN{num}(2,i)*x_xi)/J; % derivative of shape function wrt y
                
            end
            
            for i = 1:obj.numNodes
                
                Bv(:,2*i-1:2*i) = [Na_x(i),     0; ...  % B matrix in Voigt form
                                     0,     Na_y(i); ...
                                  Na_y(i),  Na_x(i)];
                              
                B(1,i) = Na_x(i); % B matrix in non-Voigt form
                B(2,i) = Na_y(i);
                
            end
            
        end % Create a B Matrix
        
        function Bv0 = B0Voigt(obj,element,num)
            
            [~,F] = element.defGrad(num); % Call for deformation gradient
            B = element.B0{num}; % Choose template B
            Bv0 = zeros(3,2*obj.numNodes);
            for i = 1:obj.numNodes
                
                Bv0(:,2*i-1:2*i) = [B(1,i)*F(1,1),                 B(1,i)*F(1,2); ...  % B matrix in Voigt form
                                    B(2,i)*F(2,1),                 B(2,i)*F(2,2); ...
                            (B(1,i)*F(2,1)+B(2,i)*F(1,1)), (B(1,i)*F(2,2)+B(2,i)*F(1,2))];
                              
            end
            
        end % Create B0 Matrix in Voigt form for TL calculations
        
        function J = Jacobian(obj,xyz,parent,num)
            % Inputs: xyz - 3 x 1 ARRAY of FLOATS
            %         parent - parent element object
            %         num - current integration point (INT)
            
            x_xi = 0; % FLOAT
            x_eta = 0; % FLOAT
            y_xi = 0; % FLOAT
            y_eta = 0; % FLOAT           
            
            for i = 1:obj.numNodes
                
                x_xi = x_xi + parent.getdN{num}(1,i)*xyz(i,1); % derivative of x wrt xi
                x_eta = x_eta + parent.getdN{num}(2,i)*xyz(i,1); % derivative of x wrt eta
                y_xi = y_xi + parent.getdN{num}(1,i)*xyz(i,2); % derivative of y wrt xi
                y_eta = y_eta + parent.getdN{num}(2,i)*xyz(i,2); % derivative of y wrt eta
                
            end % components of element to parent jacobian matrix
            
            J = x_xi*y_eta - x_eta*y_xi; % jacobian determinent (FLOAT)
        end % calculate only the Jacobian determinent
        
    end
    
end
    