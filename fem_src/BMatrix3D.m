classdef BMatrix3D
    
    properties
        
        numNodes; % number of nodes (INT)
        
    end
    
    methods
        
        function obj = BMatrix3D(numNodes)
            % Inputs: numNodes - INT            
            obj.numNodes = numNodes;
                       
        end % Constructor
        
        function [Bv,B,J] = BMatrix(obj,xyz,parent,num)
            % Inputs: xyz - 3 x 1 ARRAY of FLOATS
            %         parent - parent element object
            %         num - current integration point (INT) 
            
            Na_x = zeros(1,obj.numNodes); % 1 x numNodes ARRAY of FLOATS
            Na_y = zeros(1,obj.numNodes); % 1 x numNodes ARRAY of FLOATS
            Na_z = zeros(1,obj.numNodes); % 1 x numNodes ARRAY of FLOATS
            Bv = zeros(6,3*obj.numNodes); % 6 x 3*numNodes MATRIX of FLOATS
            B = zeros(3,obj.numNodes); % 3 x numNodes MATRIX of FLOATS
            x_xi = 0; % FLOAT
            x_eta = 0; % FLOAT
            x_zeta = 0; % FLOAT
            y_xi = 0; % FLOAT
            y_eta = 0; % FLOAT
            y_zeta = 0; % FLOAT
            z_xi = 0; % FLOAT
            z_eta = 0; % FLOAT
            z_zeta = 0; % FLOAT
            
            for i = 1:obj.numNodes
                
                x_xi = x_xi + parent.getdN{num}(1,i)*xyz(i,1); % derivative of x wrt xi
                x_eta = x_eta + parent.getdN{num}(2,i)*xyz(i,1); % derivative of x wrt eta
                x_zeta = x_zeta + parent.getdN{num}(3,i)*xyz(i,1); % derivative of x wrt zeta
                y_xi = y_xi + parent.getdN{num}(1,i)*xyz(i,2); % and so on
                y_eta = y_eta + parent.getdN{num}(2,i)*xyz(i,2);
                y_zeta = y_zeta + parent.getdN{num}(3,i)*xyz(i,2);
                z_xi = z_xi + parent.getdN{num}(1,i)*xyz(i,3);
                z_eta = z_eta + parent.getdN{num}(2,i)*xyz(i,3);
                z_zeta = z_zeta + parent.getdN{num}(3,i)*xyz(i,3);
                
            end % components of element to parent jacobian matrix
            
            % Jacobian cofactors
            cof11 = y_eta*z_zeta - y_zeta*z_eta; % FLOAT
            cof12 = y_zeta*z_xi - y_xi*z_zeta; % FLOAT
            cof13 = y_xi*z_eta - y_eta*z_xi; % FLOAT
            cof21 = z_eta*x_zeta - z_zeta*x_eta; % FLOAT
            cof22 = z_zeta*x_xi - z_xi*x_zeta; % FLOAT
            cof23 = z_xi*x_eta - z_eta*x_xi; % FLOAT
            cof31 = x_eta*y_zeta - x_zeta*y_eta; % FLOAT
            cof32 = x_zeta*y_xi - x_xi*y_zeta; % FLOAT
            cof33 = x_xi*y_eta - x_eta*y_xi; % FLOAT
            
            J = x_xi*cof11 + x_eta*cof12 + x_zeta*cof13; % Jacobian determinent (FLOAT)

            for i = 1:obj.numNodes
                
                Na_x(i) = (parent.getdN{num}(1,i)*cof11 + parent.getdN{num}(2,i)*cof12 + parent.getdN{num}(3,i)*cof13)/J; % derivative of shape function wrt x
                Na_y(i) = (parent.getdN{num}(1,i)*cof21 + parent.getdN{num}(2,i)*cof22 + parent.getdN{num}(3,i)*cof23)/J; % derivative of shape function wrt y
                Na_z(i) = (parent.getdN{num}(1,i)*cof31 + parent.getdN{num}(2,i)*cof32 + parent.getdN{num}(3,i)*cof33)/J; % derivative of shape function wrt z
                
            end
            
            for i = 1:obj.numNodes
                
                Bv(:,3*i-2:3*i) = [Na_x(i)      0          0; ...  % B matrix in Voigt form
                                     0        Na_y(i)      0; ...
                                     0          0        Na_z(i); ...
                                     0        Na_z(i)    Na_y(i); ...
                                   Na_z(i)      0        Na_x(i); ...  
                                   Na_y(i)    Na_x(i)      0];
                              
                B(1,i) = Na_x(i); % B matrix in non-Voigt form
                B(2,i) = Na_y(i);
                B(3,i) = Na_z(i);
                
            end
            
        end % Create a B Matrix 

        function Bv0 = B0Voigt(obj,element,num)
            
            [~,F] = element.defGrad(num); % Call for deformation gradient
            B = element.B0{num}; % Choose template B
            Bv0 = zeros(6,3*obj.numNodes);
            for i = 1:obj.numNodes
                
                Bv0(:,3*i-2:3*i) = [B(1,i)*F(1,1),                 B(1,i)*F(1,2),                 B(1,i)*F(1,3); ...  % B matrix in Voigt form
                                    B(2,i)*F(2,1),                 B(2,i)*F(2,2),                 B(2,i)*F(2,3); ...
                                    B(3,i)*F(3,1),                 B(3,i)*F(3,2),                 B(3,i)*F(3,3); ...
                            (B(2,i)*F(3,1)+B(3,i)*F(2,1)), (B(2,i)*F(3,2)+B(3,i)*F(2,2)), (B(2,i)*F(3,3)+B(3,i)*F(2,3)); ...
                            (B(1,i)*F(3,1)+B(3,i)*F(1,1)), (B(1,i)*F(3,2)+B(3,i)*F(1,2)), (B(1,i)*F(3,3)+B(3,i)*F(1,3)); ...
                            (B(1,i)*F(2,1)+B(2,i)*F(1,1)), (B(1,i)*F(2,2)+B(2,i)*F(1,2)), (B(1,i)*F(2,3)+B(2,i)*F(1,3))];
                            
                              
            end
            
        end % Create B0 Matrix in Voigt form for TL calculations
        
        function J = Jacobian(obj,xyz,parent,num)
            % Inputs: xyz - 3 x 1 ARRAY of FLOATS
            %         parent - parent element object
            %         num - current integration point (INT)           

            x_xi = 0; % FLOAT
            x_eta = 0; % FLOAT
            x_zeta = 0; % FLOAT
            y_xi = 0; % FLOAT
            y_eta = 0; % FLOAT
            y_zeta = 0; % FLOAT
            z_xi = 0; % FLOAT
            z_eta = 0; % FLOAT
            z_zeta = 0; % FLOAT
            
            for i = 1:obj.numNodes
                
                x_xi = x_xi + parent.getdN{num}(1,i)*xyz(i,1); % derivative of x wrt xi
                x_eta = x_eta + parent.getdN{num}(2,i)*xyz(i,1); % derivative of x wrt eta
                x_zeta = x_zeta + parent.getdN{num}(3,i)*xyz(i,1); % derivative of x wrt zeta
                y_xi = y_xi + parent.getdN{num}(1,i)*xyz(i,2); % and so on
                y_eta = y_eta + parent.getdN{num}(2,i)*xyz(i,2);
                y_zeta = y_zeta + parent.getdN{num}(3,i)*xyz(i,2);
                z_xi = z_xi + parent.getdN{num}(1,i)*xyz(i,3);
                z_eta = z_eta + parent.getdN{num}(2,i)*xyz(i,3);
                z_zeta = z_zeta + parent.getdN{num}(3,i)*xyz(i,3);
                
            end % components of element to parent jacobian matrix
            
            % Jacobian cofactors
            cof11 = y_eta*z_zeta - y_zeta*z_eta; % FLOAT
            cof12 = y_zeta*z_xi - y_xi*z_zeta; % FLOAT
            cof13 = y_xi*z_eta - y_eta*z_xi; % FLOAT
            
            J = x_xi*cof11 + x_eta*cof12 + x_zeta*cof13; % Jacobian determinent (FLOAT)
            
        end % just calculates the Jacobian determinent
    end
    
end
    