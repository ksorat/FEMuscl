classdef LinearElasticTL < MaterialModel
    
    properties
        
            lambda; % Lame's First Parameter (FLOAT)
            mu; % Lame's Second Parameter (shear modulus) (FLOAT)
            
    end
    
    methods
        
        function obj = LinearElasticTL(rho,ID,E,nu)
            % Inputs: ID - INT
            %         rho - FLOAT
            %         E - FLOAT
            %         nu - FLOAT
            
            obj@MaterialModel(rho,ID); % call superclass
            obj.lambda = (E*nu)/((1+nu)*(1-2*nu)); % Lame's First Parameter
            obj.mu = E/(2*(1+nu)); % Lame's Second Parameter (shear modulus)
            obj.c = sqrt(E*(1-nu)/((1+nu)*(1-2*nu)*rho)); % sound speed (FLOAT)
            obj.C = [(obj.lambda+2*obj.mu)  obj.lambda        obj.lambda          0       0       0; ...
                          obj.lambda   (obj.lambda+2*obj.mu)  obj.lambda          0       0       0; ...
                          obj.lambda        obj.lambda   (obj.lambda+2*obj.mu)    0       0       0; ...
                               0                 0                 0            obj.mu    0       0; ...
                               0                 0                 0              0     obj.mu    0; ...
                               0                 0                 0              0       0     obj.mu]; % Constitutive Matrix (6 x 6 MATRIX of FLOATS)           
            
        end % Constructor
        
        function [f,sig] = fInternal(obj,element,num)            
            
            [H,F] = element.defGrad(num); % gradient of displacement field and deformation gradient (both are 3 x 3 MATRICES of FLOATS)
            J = det(F); % Jacobian determinent (FLOAT)
            E = element.greenStrain(H); % Green Strain Tensor
            Ev = [E(1,1);E(2,2);E(3,3);2*E(2,3);2*E(1,3);2*E(1,2)]; % Green strain tensor in voigt notation (6 x 1 MATRIX of FLOATS)
            Sv = obj.C*Ev; % PK2 tensor in voigt notation (6 x 1 MATRIX of FLOATS)
            S = [Sv(1),Sv(4),Sv(5);...
                 Sv(4),Sv(2),Sv(6);...
                 Sv(5),Sv(6),Sv(3)]; % PK2 tensor (3 x 3 MATRIX of FLOATS)
            sig = (1/J)*F*S*F'; % stress tensor (3 x 3 MATRIX of FLOATS)
            P = S*F'; % Nominal Stress (analagous to Engineering Stress) (3 x 3 MATRIX of FLOATS)
            f = element.B0{num}'*P*element.J0(num); % internal force vector (#nodes x 3 MATRIX of FLOATS)
            f = reshape(f',size(element.DOFids)); % reshape into a column vector (3*#nodes x 1 VECTOR of FLOATS)           
            
        end
            
            
        
    end
    
end
            
            