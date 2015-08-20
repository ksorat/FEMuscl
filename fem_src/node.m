classdef node
% Create a node object. Contains coordinates and DOF IDs.

    properties        
                
        ID; % node ID (INT)
        XYZ = zeros(3,1); % Material Coordinates (3 x 1 ARRAY of FLOATS)
        xyz = zeros(3,1); % Spatial Coordinates (3 x 1 ARRAY of FLOATS)
        DOFids = zeros(3,1); % Degree of freedom catalog (3 x 1 ARRAY of FLOATS)
        DOFall; % DOF count for nodal projections
        SPCids = zeros(3,1); % Constraint catelog (1 for on, 0 for off) (3 x 1 ARRAY of INTs or some boolean)
        dim; % Dimension of the model (2D, 3D) (INT)
        owners = zeros(100,1); % List of elements that own this node (100 x 1 ARRAY of INTs) 100 allowing for an oddly large amount of elements sharing this node
        neighbors; % neighboring nodes and the elements that hold them
        damage = 0; % Damage value at node (FLOAT)
        sig; % stress tensor at node (see constructor)
        
    end
    
    methods
        
        function obj = node(ID,X,Y,varargin)
            % Inputs: X,Y,varargin - FLOATS
            % varargin is MATLAB's way of allowing a function to accept a
            % variable number of inputs to a function. This node could be
            % defined with just 2 inputs for a 2D model if we felt like it.
            
            global dim;
            
            if (nargin > 0)
                obj.ID = ID;
                obj.XYZ(1) = X; 
                obj.XYZ(2) = Y;
                obj.xyz(1) = X;
                obj.xyz(2) = Y;
                obj.dim = dim;
                
                if dim == 3                    
                    obj.XYZ(3) = varargin{1};
                    obj.xyz(3) = varargin{1};
                    obj.sig = zeros(6,1); % setup sig for 3D Voigt tensor (6 x 1 MATRIX of FLOATS)                    
                else
                    obj.sig = zeros(3,1); % setup sig for 2D Voigt tensor (3 x 1 MATRIX of FLOATS)
                end
                
            end
            
        end % Constructor        
        
        function DOF = getDOF(obj) % Get DOF ids
            
            DOF = obj.DOFids;
            
        end
        
        function obj = setxyz(obj,xyz) % Set spatial coordinates
            
            global dim;
            obj.xyz(1) = xyz(1);
            obj.xyz(2) = xyz(2);
            if dim == 3                
                obj.xyz(3) = xyz(3);
            end
            
        end
        
        function obj = setDOFs(obj) % Set DOF ids
            % Updates the global DOF counter (GDOF) and assigns the DOF IDs
            % to the node.
            global GDOF GDOFall dim;
            for i = 1:dim
                
                if ~obj.SPCids(i)
                    GDOF = GDOF + 1;
                    obj.DOFids(i) = GDOF;
                end
                
            end            
            
            obj.DOFall = zeros((3*dim-3),1);
            
            for i = 1:(3*dim-3)
                
                GDOFall = GDOFall + 1;
                obj.DOFall(i) = GDOFall;
                
            end
            
        end
        
        function obj = setSPCs(obj,dof) % Set SPC constraint on a specific DOF
            
            obj.SPCids(dof) = 1;
            
        end 
        
        function xyz = getSCoords(obj)
            
            xyz = obj.xyz;
            
        end % Get spatial coordinates

        function xyz = getMCoords(obj)
            
            xyz = obj.XYZ;
            
        end % Get material coordinates
        
        function obj = setOwners(obj,owner)
            
            last = find(obj.owners,1,'last'); % look for a single last non-zero entry in the owners vector (INT)
            
            if isempty(last) % if "last" comes out as an empty set
                idx = 1; % then this node is being assigned its first owner (INT)
            else
                idx = last + 1; % otherwise, set the target index to the next open owner slot (INT)
            end
            
            obj.owners(idx) = owner; % fill in the owner element ID
            
        end % set element owners
        
        function obj = setNeighbors(obj)
            
            global Elements Parts
            
            obj.owners = obj.owners(obj.owners>0); % reduce the vector down to the size of the # of owners
            obj.neighbors = zeros(length(obj.owners),4); % list of neighbor nodes and element that contains them (#owners x 4 ARRAY of INTs)
                                    
            for i = 1:length(obj.owners)
                
                obj.neighbors(i,1) = obj.owners(i); % mark the element ID that contains these neighbor nodes
                nodes = Elements{obj.owners(i)}.nodes; % nodes contained in this owner element
                place = find(nodes == obj.ID); % query which node in the owner element that this node is (INT)
                connections = Parts{Elements{obj.owners(i)}.part}.elForm.connections; % read the connection map from the given element formulation
                [~,numConn] = size(connections); % read the number of columns in the connection map, telling us the number of other nodes connected to this node in the element
                
                for j = 1:numConn
                   
                    obj.neighbors(i,j+1) = nodes(connections(place,j)); % mark down the neighbors contributed by this element
                    
                end                
                
            end            
            
        end % finds neighboring nodes; built this as part of my first way to split nodes; now obsolete
        
        function obj = split(obj)
            
            global Elements Nodes u v a Fext dim Lengths

            if dim == 2

                stress = [obj.sig(1) obj.sig(3);...
                          obj.sig(3) obj.sig(2)]; % 2D expanded stress tensor at node (2 x 2 MATRIX of FLOATs)

            elseif dim == 3

                stress = [obj.sig(1) obj.sig(6) obj.sig(5);...
                          obj.sig(6) obj.sig(2) obj.sig(4);...
                          obj.sig(5) obj.sig(4) obj.sig(3)]; % 3D expanded stress tensor at node (3 x 3 MATRIX of FLOATs)

            end

            [PrincVec,PrincStress] = eig(stress); % solve the eigenvalue problem to produce the principal stresses and their directions                
            PrincStress = diag(PrincStress); % convert to a vector
            [~,idx] = max(abs(PrincStress)); % find the max principal stress
            vector = PrincVec(:,idx); % pull out the max principal stress vector; note: MATLAB ensures this is a unit vector during its eigen analysis
            if dim == 2
                vector = [vector;0]; % add a z value onto the eigen vector to remain consistent with the xyz format of the rest of the code
            end

            P1 = obj.xyz + 0.5*min(Lengths)*vector; % create and displace first phantom point
            P2 = obj.xyz - 0.5*min(Lengths)*vector; % create and displace second phantom point                
            splitOrStay = zeros(size(obj.owners)); % a boolean array flagging an owner who will (1) split or (0) stay

            for i = 1:length(obj.owners)

                center = mean(Elements{obj.owners(i)}.xyz)'; % determine the average center value of the nodes that comprise this owner element
                d1 = norm(center - P1); % distance from the center to P1
                d2 = norm(center - P2); % distance from the center to P2

                if d2 <= d1                        
                    splitOrStay(i) = 1; % mark owner element for splitting                        
                end

            end
            
            if sum(splitOrStay) == 0 % if there aren't any elements that can split, bail out
                return
            end
            
            if sum(splitOrStay) == length(splitOrStay) % if every element wants to split, also bail out
                return
            end

            splitOrStay = logical(splitOrStay); % convert the vector into a boolean array (gotta do it this way in MATLAB)

            % Copy the original Node to paste onto the new clone Node
            tempNode = obj; % create the new node object by copying this current node
            tempNode = tempNode.setDOFs; % expand global DOFs and assign them to new node
            Nodes{end+1} = tempNode; % Add new node onto master node list
            newID = length(Nodes); % ID of the new copied node
            
%             if newID==256
%                 keyboard
%             end
            
            Nodes{newID}.ID = newID; % assign the correct ID to the new Node object

            % Initialize copy-vectors for the global attributes of the
            % new Node object
            tempU = zeros(dim,1); % temporary displacements vector
            tempV = zeros(dim,1); % temporary velocity vector
            tempA = zeros(dim,1); % temporary acceleration vector
            tempF = zeros(dim,1); % temporary external force vector

            % Pull out the original Node's information from the global
            % vectors
            for j = 1:dim

                tempU(j) = u(obj.DOFids(j)); % copy displacements
                tempV(j) = v(obj.DOFids(j)); % copy velocities
                tempA(j) = a(obj.DOFids(j)); % copy accelerations
                tempF(j) = Fext(obj.DOFids(j)); % copy external forces

            end

            % Tack the new Node's information back onto the global
            % vectors
            % sorry, I don't know how to dynamically expand the
            % size of the global history variable vectors in C.
            % MATLAB makes it easy:
            u = [u;tempU]; % concatenate the new displacements for new node with global displacement vector
            v = [v;tempV]; % concatenate the new velocities for new node with global velocity vector
            a = [a;tempA]; % concatenate the new accelerations for new node with global acceleration vector
            Fext = [Fext;tempF]; % concatenate the new forces for new node with global external force vector

            for i = 1:length(splitOrStay)

                if splitOrStay(i)

                    nodeNumber = find(Elements{obj.owners(i)}.nodes==obj.ID); % find the spot in the Element object's node list where we need to replace the original ID with new one
                    Elements{obj.owners(i)}.nodes(nodeNumber) = newID; % replace the ID
                    Elements{obj.owners(i)} = Elements{obj.owners(i)}.pullDOFs;
                    

                end

            end
            
            Nodes{newID}.owners = obj.owners(splitOrStay); % reassign owners in the new node
            obj.owners = obj.owners(~splitOrStay); % reassign owners in the original node
            
                         
        end % splits this node by creating a clone and assigning the clone to some of the owner elements, then updating all global history variables
        
    end
    
    methods(Static)
        
        function damageCheck()
            
            global Nodes Elements Parts M Damage GDOF
            
            splitFlag = false; % flag marking if anything has split
            element.NodeStress(); % run the static method to project stresses to the nodes
            
            for i = 1:length(Nodes)
                
                Nodes{i}.owners = Nodes{i}.owners(Nodes{i}.owners>0); % reduce the vector down to the size of the # of owners

                if length(Nodes{i}.owners) == 1 % if this node only belongs to one element, just skip it                    
                    continue                    
                end

                damageIDs = zeros(length(Nodes{i}.owners),1);
                
                for j = 1:length(Nodes{i}.owners)

                    damageIDs(j) = Parts{Elements{Nodes{i}.owners(j)}.part}.damage; % mark down the damage ID from this element
                    
                end
                
                uniqueDamage = unique(damageIDs); % pull out the unique damage models
                
                if length(uniqueDamage) > 1 % if there is more than one damage model associated with the node, just skip it                     
                    continue                    
                end
                
                try
                nodeDamage = Damage{uniqueDamage}.damageCalc(Nodes{i}.sig);
                catch err
                    keyboard
                end

                 if nodeDamage >= 1

                    splitFlag = true; % mark the flag as true
                    Nodes{i} = Nodes{i}.split(); % tell the node to split
                    %fprintf('\tNode %d split! (Better go catch it)\n', i);
                    global RedoS2G;
                    RedoS2G = true; %Signal a recalculation of surface geometry
                end
                
            end
            
                % Recalculate the Global Mass Matrix to redistribute the
                % mass of the system with the new nodes involved
                
                if splitFlag
                    
                    M = zeros(GDOF); % reinitialize global mass matrix to new size due to increase in DOFs

                    for i = 1:length(Elements)                    
                        Elements{i} = Elements{i}.massMatrix(Parts{Elements{i}.part}.elForm); % run mass matrix routine on all elements
                    end % rebuild the global mass matrix
                    M = diag(M); % turn M matrix into a M vector (a GDOF x 1 VECTOR of FLOATS)
                    
                end
            
        end
        
    end
end



























