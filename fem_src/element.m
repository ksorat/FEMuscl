classdef element
% Create an Element object. This is the lifeblood of the simulation. The 
% element object holds state data, and pointers to node objects, and a
% pointer to the part it resides in. It holds a copy of its node
% coordinates and node degrees of freedom but that probably isn't
% completely necessary. The B matrix formulation is variable upon startup
% for 2D and 3D implementation; this decision is made upon construction, 
% not in an if statement during the solve. Element forces are solved here 
% and deformation and strain tensors are also computed.

    properties
        
        ID; % Element ID (INT OR LONG)
        nodes; % array of node IDs (VECTOR of INTs or LONGs)
        part; % part ID (INT)
        rho; % mass density (FLOAT)
        sig; % stress values (num integration pnt x 1 ARRAY of 2x2 or 3x3 FLOAT ARRAYS)
        m; % mass matrix (dim*#nodes x dim*#nodes FLOAT ARRAY)
        B; % B Matrix object for creating correct B Matrix for this element
        DOFids; % DOF IDs of the contained nodes (dim*#nodes x 1 FLOAT VECTOR)
        damage; % damage values for element (num integration pnt x 1 ARRAY FLOATS)
        xyz; % spatial coords ( 3 x 1 FLOAT VECTOR)
        XYZ; % material coords ( 3 x 1 FLOAT VECTOR)
        l = 10000; % characteristic length (FLOAT)
        Bv0; % B Matrix (Voigt form) evaluated at material coords (only needs to be done once!) (#int pnts x 1 ARRAY of either 3 x 2*#nodes MATRICES of FLOATS or 6 x 3*#nodes MATRICES of FLOATS)
        B0; % B Matrix evaluated at material coords (only needs to be done once!) (#int pnts x 1 ARRAY of dim x #nodes MATRICES of FLOATS)
        J0; % Jacobian determinant at material coords (only needs to be done once!) (#int pnts x 1 ARRAY of FLOATS)
        ue; % Working element displacement (originally created as a dim*#nodes x 1 VECTOR of FLOATS; reshaped via matlab magic into a dim x #nodes MATRIX of FLOATS)
        ve; % Working element velocity (originally created as a dim*#nodes x 1 VECTOR of FLOATS; reshaped via matlab magic into a dim x #nodes MATRIX of FLOATS)
        
    end
    
    methods
        
        function obj = element(ID,nodes,part)
        % Inputs: ID - INT or LONG
        %         nodes - INT or LONG ARRAY of node IDs
        %         part - INT
        % Assigns ID, nodes, and part. Automatically pulls out initial
        % coordinates, degrees of freedom, mass density, assigns a B Matrix
        % solver, and solves for the unchanging B0 matrices and element
        % Jacobian determinents.
            global Parts dim Nodes;
            obj.ID = ID;
            numNodes = length(nodes);
            obj.nodes = nodes;
            for i = 1:length(nodes)
               
                Nodes{nodes(i)} = Nodes{nodes(i)}.setOwners(ID); % set this element as an owner of this node
                
            end
            obj.part = part;
            obj = obj.pullMCoords; % pull out XYZ coords from nodes
            obj.xyz = obj.XYZ; % Upon initialization spatial coords = material coords
            obj = obj.pullDOFs; % pull out DOFs from nodes
            obj.rho = Parts{part}.material.rho; % element mass density
            if dim==2
                obj.B = BMatrix2D(numNodes); % assign 2D B Matrix solver to the element               
            elseif dim==3
                obj.B = BMatrix3D(numNodes); % assign 3D B Matrix solver to the element
            end
            parent = Parts{part}.elForm; % query part for element form
            dN = parent.getdN; % shape function derivatives for creating the B0 Matrices (#int pnts x 1 ARRAY of dim x #nodes ARRAY of FLOATS)
            obj.B0 = cell(length(dN),1);
            obj.Bv0 = obj.B0;
            obj.J0 = zeros(length(dN),1);
            obj.sig = cell(length(dN),1);
            obj.damage = zeros(length(dN),1);
            
            for i = 1:length(dN)
                
                [obj.Bv0{i},obj.B0{i},obj.J0(i)] = obj.B.BMatrix(obj.XYZ,parent,i);
                
            end % populate B0 Matrices
            
            obj = obj.massMatrix(parent); % create element mass matrix            
            
        end % Constructor
        
        function obj = pullSCoords(obj)
            
            global Nodes;
            
            obj.xyz = zeros(length(obj.nodes),3);
            
            for i = 1:length(obj.nodes)
                
                obj.xyz(i,:) = Nodes{obj.nodes(i)}.getSCoords;
                
            end
        
        end % Pull spatial coords
        
        function obj = pullMCoords(obj)
            
            global Nodes;
            
            obj.XYZ = zeros(length(obj.nodes),3);
            
            for i = 1:length(obj.nodes)
                
                obj.XYZ(i,:) = Nodes{obj.nodes(i)}.getMCoords;
                
            end
        
        end % Pull material coords
        
        function obj = pullDOFs(obj)
            
            global dim Nodes;
            
            obj.DOFids = zeros(dim*length(obj.nodes),1);
            count = 0;
            
            for i = 1:length(obj.nodes)
                
                temp = Nodes{obj.nodes(i)}.getDOF;
                
                for j = 1:dim
                    
                    count = count + 1;
                    obj.DOFids(count) = temp(j);
                    
                end
                
            end
            
        end % Pull DOF IDs
        
        function obj = massMatrix(obj,parent)
            
            global M;
            
            obj.m = zeros(length(obj.DOFids));
            
            for i = 1:length(parent.N)
                
                mi = obj.rho*parent.NMatrix{i}'*parent.NMatrix{i}*parent.Lw(i)*obj.J0(i); % evaluate mass matrix at integration point
                obj.m = obj.m + mi; % sum integration point values together
            
            end % Lobatto quadrature
            
            for i = 1:length(obj.m)
                
                if obj.DOFids(i) == 0                    
                    continue                    
                end
                
                for j = 1:length(obj.m)
                    
                    if obj.DOFids(j) == 0                    
                        continue                    
                    end
                
                    M(obj.DOFids(i),obj.DOFids(j)) = M(obj.DOFids(i),obj.DOFids(j)) + obj.m(i,j);
                    
                end
                                
            end % scatter nodal forces into global vector
            
        end % Create element mass matrix and scatter into global M            
        
        function [Bv,B,J] = BMatrix(obj,num)
            global Parts;
            
            [Bv,B,J] = obj.B.BMatrix(obj.xyz,Parts{obj.part}.elForm,num);            
            
        end % Call B Matrix
        
        function sig = VoigtStress(obj,num)
            
            global dim
            
            switch dim
                
                case 2
                    
                    sig = [obj.sig{num}(1,1);obj.sig{num}(2,2);obj.sig{num}(1,2)];
                    
                case 3
                    
                    sig = [obj.sig{num}(1,1);obj.sig{num}(2,2);obj.sig{num}(3,3);obj.sig{num}(2,3);obj.sig{num}(1,3);obj.sig{num}(1,2)];
                    
            end 
            
        end % put stress in Voigt form            
            
        function obj = charLength(obj)
            
            global Lengths;

            for i = 1:length(obj.nodes)
                
                for j = 1:length(obj.nodes)
                    
                    if j ~= i
                        
                        d = sqrt(sum((obj.xyz(i,:) - obj.xyz(j,:)).^2));
                        
                        if d < obj.l
                            
                            obj.l = d;
                            
                        end
                        
                    end
                    
                end
                
            end
            
            Lengths(obj.ID) = obj.l;
            
        end % calculate the characteristic length of the element
        
        function [H,F] = defGrad(obj,num)
            
            global dim
                                    
            H = obj.ue*obj.B0{num}'; % gradient of displacement field
            F = eye(dim) + H; % deformation gradient
            
        end % calculate the deformation gradient        
        
        function [B,L] = velGrad(obj,num)
                                    
            B = obj.Bmatrix(obj,num);
            L = obj.ve*B';
            
        end % calculate B Matrix and Velocity Gradient        
        
        function [obj,Fe] = getforce(obj)

            global Parts u v Fext GDOF dim Nodes;

            obj.ue = zeros(size(obj.DOFids));
            obj.ve = obj.ue;
            fext = obj.ue; % element external force vector (dim*#nodes x 1 VECTOR of FLOATS)
            fint = obj.ue; % element internal force vector (dim*#nodes x 1 VECTOR of FLOATS)
            Fe = zeros(GDOF,1); % element contribution to global force vector (GDOF x 1 VECTOR of FLOATS)
                        
            for i = 1:length(obj.DOFids)
                
                if obj.DOFids(i) == 0 % DOFid of 0 indicates that node cannot move in that direction. Skip it.
                    continue
                end
                % using the DOFids of the nodes in this element as indices
                % of the global disp and velocity vectors
                obj.ue(i) = u(obj.DOFids(i));
                obj.ve(i) = v(obj.DOFids(i));
                fext(i) = Fext(obj.DOFids(i));
                
            end % gather displacements, velocities, external forces

            count = 0; % counter for indexing inside of coordinate update loop (INT)

            for i = 1:length(obj.nodes)
                
                for j = 1:dim
                    
                    count = count + 1;                    
                    Nodes{obj.nodes(i)}.xyz(j) = obj.XYZ(i,j) + obj.ue(count);
                    
                end
                                
            end % update node spatial coordinates
            obj = pullSCoords(obj); % update element coordinates
            
            obj.ue = reshape(obj.ue,dim,length(obj.nodes)); % rearrange displacements for force update
            obj.ve = reshape(obj.ve,dim,length(obj.nodes)); % rearrange velocities for force update            

            for i = 1:length(Parts{obj.part}.elForm.N)
                
                [ftemp,obj.sig{i}] = Parts{obj.part}.material.fInternal(obj,i); % call stress function of material model (ftemp is a dim*#nodes x 1 VECTOR of FLOATS)
                
                fint = fint + ftemp*Parts{obj.part}.elForm.w(i); % (dim*#nodes x 1 VECTOR of FLOATS)
                
                %obj.damage(i) = Parts{obj.part}.damage.damageCalc(obj,i); % calculate damage at integration point
            
            end % calculate stress and internal forces
            
            f = fext - fint; % element nodal forces (dim*#nodes x 1 VECTOR of FLOATS)
            
            for i = 1:length(obj.DOFids)
                
                if obj.DOFids(i) == 0 % DOFid of 0 indicates that node cannot move in that direction. Skip it.
                    continue
                end                
                Fe(obj.DOFids(i)) = f(i); % (GDOF x 1 VECTOR of FLOATS)
                
            end % scatter nodal forces into global vector            
            
        end % calculate internal and external forces
        
        function [m,pe] = localNodeStress(obj)
            
            global Parts dim;
            
            numInts = length(Parts{obj.part}.elForm.dN); % number of integration points
            SM = Parts{obj.part}.elForm.StressMatrix; % pull out the stress matrix
            Ns = Parts{obj.part}.elForm.Ns; % pull out the stress matrix
            mm = zeros(size(SM{1})); % volumetric weighting matrix
            pe = zeros(length(SM{1}),1); % load vector           
            
            for i = 1:numInts
                
                sigma = obj.VoigtStress(i); % convert stress at this int point to voigt form
                J = obj.B.Jacobian(obj.xyz,Parts{obj.part}.elForm,i); % Current Jacobian determinent at this int point
                pe = pe + Ns{i}'*sigma; % calc load vector contribution of this int point
                mm = mm + SM{i}*J; % calc volumetric weighting matrix contribution of this int point                
                
            end
            
            m = sum(mm,2); % row summed "lumped" vector for decoupled, faster computing later
            
        end % calculate element contribution to nodal stress components
                    
    end
    
    methods(Static)
        
        function E = greenStrain(H)
            
            E = .5*(H + H' + H'*H); % warned against calculating E from F by textbook because roundoff errors can create serious problems             
            
        end  % calculate Green Strain Tensor
        
        function D = rateDef(L)
            
            D = .5*(L' + L);
            
        end % calculate Rate of Deformation tensor
        
        function NodeStress()
           
            global Elements Nodes dim GDOFall;
            
            P = zeros(GDOFall,1); % initialize global load vector
            Ms = zeros(GDOFall,1); % initialize global volumetric weighting vector
            numberSlots = (3*dim-3);
            
            for i = 1:length(Elements)
                
                [me,pe] = Elements{i}.localNodeStress(); % extract local element load and weighting vectors
                currentNodes = Elements{i}.nodes;
                
                for j = 1:length(currentNodes)
                    
                    for k = 1:numberSlots
                        
                        P(Nodes{currentNodes(j)}.DOFall(k)) = P(Nodes{currentNodes(j)}.DOFall(k)) + pe(j*numberSlots-(numberSlots-k)); % scatter to global vector
                        Ms(Nodes{currentNodes(j)}.DOFall(k)) = Ms(Nodes{currentNodes(j)}.DOFall(k)) + me(j*numberSlots-(numberSlots-k)); % scatter to global vector
                        
                    end
                    
                end
                
                sig_nodes = P./Ms; % solve for nodal stresses
             
                lastSlot = 0;
             
                for j = 1:length(Nodes)

                    Nodes{j}.sig = sig_nodes((lastSlot + 1):(lastSlot + numberSlots)); 
                    lastSlot = lastSlot + numberSlots;

                end
                
            end
                       
        end % L2 projection of stresses to nodes
        
    end       
        
end




















