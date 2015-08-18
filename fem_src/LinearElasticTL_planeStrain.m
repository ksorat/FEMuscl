classdef LinearElasticTL_planeStrain < LinearElasticTL
% Version for 
    
    methods
        
        function obj = LinearElasticTL_planeStrain(ID,rho,E,nu)
            % Inputs: ID - INT
            %         rho - FLOAT
            %         E - FLOAT
            %         nu - FLOAT
            
            obj@LinearElasticTL(ID,rho,E,nu); % call superclass
            obj.c = sqrt(E/((1-nu^2)*rho)); % sound speed (FLOAT)
            obj.C = [(obj.lambda+2*obj.mu)      obj.lambda           0; ...
                          obj.lambda       (obj.lambda+2*obj.mu)     0; ...
                               0                     0            obj.mu]; % Constitutive Matrix (3 x 3 MATRIX of FLOATS)
            
        end % Constructor
        
        function [f,sig] = fInternal(obj,element,num)            
            % Input: element - element object
            %        num - current integration point (INT)
            
            [H,F] = element.defGrad(num); % gradient of displacement field and deformation gradient (both are 2 x 2 MATRICES of FLOATS)
            J = det(F); % Jacobian determinent (FLOAT)
            E = element.greenStrain(H); % Green Strain Tensor (2 x 2 MATRIX of FLOATS)
            Ev = [E(1,1);E(2,2);2*E(1,2)]; % Green strain tensor in voigt notation (3 x 1 VECTOR of FLOATS)
            Sv = obj.C*Ev; % PK2 tensor in voigt notation (3 x 1 VECTOR of FLOATS)
            S = [Sv(1),Sv(3);...
                 Sv(3),Sv(2)]; % PK2 tensor (2 x 2 MATRIX of FLOATS)
            sig = (1/J)*F*S*F'; % stress tensor (2 x 2 MATRIX of FLOATS)
            P = S*F'; % Nominal Stress (analagous to Engineering Stress) (2 x 2 MATRIX of FLOATS)
            f = element.B0{num}'*P*element.J0(num); % internal force vector (#nodes x 2 MATRIX of FLOATS)

            f = reshape(f',size(element.DOFids)); % reshape into a column vector (2*#nodes x 1 VECTOR of FLOATS)           
            
        end
        
    end
    
end