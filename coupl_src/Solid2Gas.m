%Takes data from the solid solver and strips it down to data for the gas
%solver



%NOTE, having <=3 owners DOES NOT guarantee boundary node with splitting
%present

function lvlDef = Solid2Gas(Model)
global RedoS2G;

if isempty(RedoS2G)
    RedoS2G = true;
end

if (RedoS2G)
    %This work only done when topology changes
    
    %Calculate surface geometry
    fprintf('\t(Re)Calculating solid surface geometry\n');
    %Calculate Link Matrix, split nodes and recalculate until topology is
    %good
    L = LinkMatrix();
    
    %Pull orbits from link matrix
    lvlDef = findOrbits(L);
    RedoS2G = false;
    
else
    lvlDef = Model.Init.lvlDef; %Use old structure
end

%This work always must be done
lvlDef = reBoundary(lvlDef); %Refresh positions/velocities of boundary

function L = LinkMatrix()

global Nodes;
goodTopology = false;

while (~goodTopology)
    %Get link matrix
    L = calcLinkMatrix();
    
    %Check link matrix
    [n1 n2 ~] = find(L);
    [nUn Iu] = unique(n1); %Find unique values
    Chk = length(n1) - length(nUn);
    if (Chk == 0)
        goodTopology = true;
    else
        fprintf('\t\tBad topology, recalculating.\n');
        %Find problem node (start w/ first)
        nBad = n1; 
        nBad(Iu) = [];
        nBad
        NumNod1 = length(Nodes);
        for nb=1:length(nBad)
            
            fixTopology( nBad(nb) );
            if (nBad(nb) == 150)
                %keyboard
            end
        end
        %fixTopology( nBad(1) );
        NumNod2 = length(Nodes);
        if (NumNod1 == NumNod2)
            keyboard
        end
        
        if ( NumNod2 > 260 )
            %keyboard
        end
        %keyboard
    end
end

%Forces a split of node nB to hopefully fix topology
function fixTopology (nB)

global Nodes Elements M GDOF Parts;


Nodes{nB} = Nodes{nB}.split();

%Recalculate mass matrix/element connectivity
M = zeros(GDOF); % reinitialize global mass matrix to new size due to increase in DOFs

for i = 1:length(Elements)
    if isempty(Elements{i}) continue; end
    Elements{i} = Elements{i}.massMatrix(Parts{Elements{i}.part}.elForm); % run mass matrix routine on all elements
end % rebuild the global mass matrix
M = diag(M); % turn M matrix into a M vector (a GDOF x 1 VECTOR of FLOATS)


function L = calcLinkMatrix()

global Nodes Elements;

NumNod = length(Nodes);
NumElt = length(Elements);

%Start by finding boundary nodes/links (CCW)

L = spalloc(NumNod,NumNod,4*NumElt);
for n=1:NumElt
    eltn = Elements{n};
    if isempty(Elements{n}) continue; end
    enods = [ eltn.nodes eltn.nodes(1) ]; %Embiggen for all lines
    for l=1:4
        n1 = enods(l); n2 = enods(l+1);
        if ( L(n2,n1) == 1 )
            L(n2,n1) = 0; %Cancel link b/c double
        else
            L(n1,n2) = 1;
        end
        
    end
end

function lvlDef = findOrbits(L)

%Now we have sparse matrix with links (CCW)
%Each closed orbit corresponds to a different distinct object
%Traverse orbits, pull out of P
numObj = 0;
while ( nnz(L) > 0 )
    numObj = numObj+1;
    [n1 n2 ~] = find(L); %Arrows (CCW) n1->n2
    sub_id = n1(1); %Starting point of link
    nex_id = n2(1); %Ending point of link
    
    while ( nex_id ~= sub_id(1) ) %Traverse until closure
        sub_id = [sub_id nex_id];
        I = find( n1 == nex_id ); %Switch so nex_id is starting point
        if (length(I) > 1)
            disp('Topology failure!  This shouldn''t happen!');
            keyboard
        end
        nex_id = n2(I);
    end
    sub_id = [sub_id sub_id(1)]; %Close circuit
    %Pull entries from connectivity matrix
    Nsub = length(sub_id);
    for n=1:Nsub-1
        n1 = sub_id(n); n2 = sub_id(n+1);
        L(n1,n2) = 0;
    end
    lvlDef.obsDat{numObj}.nid = sub_id;
end
lvlDef.numObs = numObj;
    
function lvlDef = reBoundary(lvlDef)

global Nodes v;

numObj = lvlDef.numObs;

for no=1:numObj
    obs = lvlDef.obsDat{no};
    Nv = length(obs.nid);
    obs.xv = zeros(1,Nv);
    obs.yv = obs.xv;
    obs.vx = obs.xv;
    obs.vy = obs.xv;
    obs.vxid = obs.xv;
    obs.vyid = obs.xv;
        
    for n=1:Nv
        nod = Nodes{obs.nid(n)};
        obs.xv(n) = nod.xyz(1);
        obs.yv(n) = nod.xyz(2);
        vxid = nod.DOFids(1);
        vyid = nod.DOFids(2);
        
        if (vxid ~= 0)
            obs.vx(n) = v(vxid);
        else
            obs.vx(n) = 0;
        end
        if (vyid ~= 0)
            obs.vy(n) = v(vyid);
        else
            obs.vy(n) = 0.0;
        end
        %Save return mapping
        obs.vxid(n) = vxid;
        obs.vyid(n) = vyid;
    end
    lvlDef.obsDat{no} = obs; %Return
end



% function lvlDef = Solid2Gas(Nodes,Elements,v,Model)
% 
% global RedoS2G;
% 
% if isempty(RedoS2G)
%     RedoS2G = true;
% end
% 
% if (RedoS2G) %Calculate surface geometry
%     fprintf('\t(Re)Calculating solid surface geometry\n');
%     
%     NumNod = length(Nodes);
%     NumElt = length(Elements);
%     
%     %Start by finding boundary nodes/links (CCW)
%     
%     P = spalloc(NumNod,NumNod,4*NumElt);
%     for n=1:NumElt
%         eltn = Elements{n};
%         enods = [ eltn.nodes eltn.nodes(1) ]; %Embiggen for all lines
%         for l=1:4
%             n1 = enods(l); n2 = enods(l+1);
%             if ( P(n2,n1) == 1 )
%                 P(n2,n1) = 0; %Cancel link b/c double
%             else
%                 P(n1,n2) = 1;
%             end
%             
%         end
%     end
%     
%     
% 
%     RedoS2G = false;
%     
% else %Use surface geometry from lvlDef
%     lvlDef = Model.Init.lvlDef; %Use old structure
%     
% end


% %Note, all this previous works only needs to be done when the topology
% %changes.  Remaining needs to be done always and only relies on
% %lvlDef.obsDat{n}.nid's



% function lvlDef = Solid2Gas(Nodes,Elements,v,Model)
% 
% global RedoS2G
% 
% if isempty(RedoS2G)
%     RedoS2G = true;
% end
% 
% if (RedoS2G) %Calculate surface geometry
%     fprintf('\t(Re)Calculating solid surface geometry\n');
%     NumNod = length(Nodes);
%     NumElt = length(Elements);
%     
%     
%     %Start by finding boundary nodes, this may be overcomplicated
%     
%     P = spalloc(NumNod,NumNod,4*NumElt); %Overestimate
%     for n=1:NumElt
%         eltn = Elements{n};
%         enods = [eltn.nodes eltn.nodes(1)]; %Embiggen for all lines
%         for l=1:4
%             n1 = enods(l); n2 = enods(l+1);
%             nM = max(n1,n2); nm = min(n1,n2); %Wolog, nM>nm (can't be equal)
%             P(nm,nM) = P(nm,nM) + 1;
%         end
%     end
%     
%     %At this point, P(n,m) = 0,1,2 = no connection, 1 connection, dbl
%     %connection.  Boundary nodes correspond to 1 connection
%     cInd = (P == 2);
%     P(cInd) = 0;
%     
%     [i,j,Ps] = find(P); %Find nonzero entries
%     
%     nid = sort(unique([i;j]))'; %These are the boundary nodes!
%     %Get boundary elements
%     ebd = [];
%     for n=1:length(nid)
%         nod = Nodes{ nid(n) };
%         owned = nod.owners;
%         %Clean out zeros
%         owned(owned == 0) = [];
%         ebd = [ebd owned'];
%         
%     end
%     ebd = sort(unique(ebd));
%     
%     Nv = length(nid); %Number of boundary nodes
%     next = zeros(size(nid));
%     %Now loop through boundary elements looking for connectivity
%     
%     for e=1:length(ebd)
%         elt = Elements{ebd(e)};
%         enods = [elt.nodes elt.nodes(1)]; %Embiggen for closed curve
%         %Loop through nodes, find CW connectivity (assuming CW listing)
%         
%         for n=1:4
%             enn = enods(n); ennp = enods(n+1); %Find line segment
%             
%             if ( ismember(enn,nid) & ismember(ennp,nid) )
%                 
%                 %This line segment connects two boundary nodes
%                 I = find( nid == enn); %Index of node enn
%                 if ( next(I) == 0 )
%                     next(I) = ennp; %Set ennp to be after enn
%                 else
%                     %There was already somebody there, which do I choose?
%                     
%                     disp('Double found next');
%                 end
%                 
%             end
%         end
%         
%     end
%     
%     
%     numObj = 0;
%     tot_nid = nid;
%     tot_next = next;
%     %Now find orbits until the list is exhausted
%     while ( length(tot_nid) > 0 )
%         numObj = numObj+1;
%         sub_id = tot_nid(1); %nid(1) = [];
%         nex_id = tot_next(1); %next(1) = [];
%         while ( nex_id ~= sub_id(1) )
%             if ( ismember(sub_id,nex_id) )
%                 %You've hit a circle
%                 nex_id = sub_id(1);
%             else
%                 sub_id = [sub_id nex_id];
%                 I = find(tot_nid == nex_id);
%                 nex_id = tot_next(I);
%             end
%             
%         end
%         if ( length(sub_id) < 4 ) %Stupid orbit
%             numObj = numObj-1;
%         else
%             lvlDef.obsDat{numObj}.nid = sub_id;
%         end
%         
%         for n=1:length(sub_id) %Remove entries from nid
%             idn = sub_id(n);
%             I = find( tot_nid == idn);
%             tot_nid(I) = [];
%             tot_next(I) = [];
%         end
%         
%     end
%     lvlDef.numObs = numObj;
%     RedoS2G = false;
% else
%     lvlDef = Model.Init.lvlDef; %Use old structure
% end
% %Note, all this previous works only needs to be done when the topology
% %changes.  Remaining needs to be done always and only relies on
% %lvlDef.obsDat{n}.nid's
% 
% numObj = lvlDef.numObs;
% 
% for no=1:numObj
%     obs = lvlDef.obsDat{no};
%     Nv = length(obs.nid);
%     obs.xv = zeros(1,Nv);
%     obs.yv = obs.xv;
%     obs.vx = obs.xv;
%     obs.vy = obs.xv;
%     obs.vxid = obs.xv;
%     obs.vyid = obs.xv;
%         
%     for n=1:Nv
%         nod = Nodes{obs.nid(n)};
%         obs.xv(n) = nod.xyz(1);
%         obs.yv(n) = nod.xyz(2);
%         vxid = nod.DOFids(1);
%         vyid = nod.DOFids(2);
%         
%         if (vxid ~= 0)
%             obs.vx(n) = v(vxid);
%         else
%             obs.vx(n) = 0;
%         end
%         if (vyid ~= 0)
%             obs.vy(n) = v(vyid);
%         else
%             obs.vy(n) = 0.0;
%         end
%         %Save return mapping
%         obs.vxid(n) = vxid;
%         obs.vyid(n) = vyid;
%     end
%     lvlDef.obsDat{no} = obs; %Return
% end
% 
% 
% if (numObj > 1)
%     disp('Breaking!');
%     
%     %keyboard
% end

    

% 
% function lvlDef = Solid2Gas(Nodes,Elements,v)
% 
% NumNod = length(Nodes);
% 
% %Start by finding boundary nodes, and boundary elements
% nid = []; ebd = [];
% 
% for n=1:NumNod
%     nod = Nodes{n};
%     nod = nod.setNeighbors();
%     own = length(nod.owners);
%         
%     if (own <= 3) %You're a boundary node
%         nid = [nid n];
%         ebd = [ebd nod.owners'];
%         
%     end
% end
% 
% ebd = sort(unique(ebd));
% 
% Nv = length(nid); %Number of boundary nodes
% next = zeros(size(nid));
% %Now loop through boundary elements looking for connectivity
% 
% for e=1:length(ebd)
%     elt = Elements{ebd(e)};
%     enods = [elt.nodes elt.nodes(1)]; %Embiggen for closed curve
%     %Loop through nodes, find CW connectivity (assuming CW listing)
% 
%     for n=1:4
%         enn = enods(n); ennp = enods(n+1); %Find line segment
%         
%         if ( ismember(enn,nid) & ismember(ennp,nid) )
% 
%             %This line segment connects two boundary nodes
%             I = find( nid == enn); %Index of node enn
%             next(I) = ennp; %Set ennp to be after enn
% 
%         end
%     end
%     
% end
% 
% numObj = 0;
% tot_nid = nid;
% %Now find orbits until the list is exhausted
% while ( length(tot_nid) > 0 )
%     numObj = numObj+1;
%     sub_id = tot_nid(1); %nid(1) = [];
%     nex_id = next(1); %next(1) = [];
%     while ( nex_id ~= sub_id(1) )
%         sub_id = [sub_id nex_id];
%         I = find(nid == nex_id);
%         nex_id = next(I);
%         %Remove entries
%         
%     end
%     lvlDef.obsDat{numObj}.nid = sub_id;
%     
%     for n=1:length(sub_id) %Remove entries from nid
%         idn = sub_id(n);
%         I = find( tot_nid == idn);
%         tot_nid(I) = [];
%     end
%     Chk = (sub_id == 0);
%     if (sum(Chk) > 0)
%         keyboard
%     end
% end
% 
% %Note, all this previous works only needs to be done when the topology
% %changes.  Remaining needs to be done always and only relies on
% %lvlDef.obsDat{n}.nid's
% 
% for no=1:numObj
%     obs = lvlDef.obsDat{no};
%     Nv = length(obs.nid);
%     obs.xv = zeros(1,Nv);
%     obs.yv = obs.xv;
%     obs.vx = obs.xv;
%     obs.vy = obs.xv;
%     obs.vxid = obs.xv;
%     obs.vyid = obs.xv;
%         
%     for n=1:Nv
%         nod = Nodes{obs.nid(n)};
%         obs.xv(n) = nod.xyz(1);
%         obs.yv(n) = nod.xyz(2);
%         vxid = nod.DOFids(1);
%         vyid = nod.DOFids(2);
%         
%         if (vxid ~= 0)
%             obs.vx(n) = v(vxid);
%         else
%             obs.vx(n) = 0;
%         end
%         if (vyid ~= 0)
%             obs.vy(n) = v(vyid);
%         else
%             obs.vy(n) = 0.0;
%         end
%         %Save return mapping
%         obs.vxid(n) = vxid;
%         obs.vyid(n) = vyid;
%     end
%     lvlDef.obsDat{no} = obs; %Return
% end
% 
% lvlDef.numObs = numObj;
% 
% if (numObj > 1)
%     disp('Breaking!');
% end

% function lvlDef = Solid2Gas_Old(Nodes,Elements,v)
% 
% NumNod = length(Nodes);
% 
% xv = [];
% yv = [];
% nid = [];
% for n=1:NumNod
%     nod = Nodes{n};
%     nod = nod.setNeighbors();
%     own = length(nod.owners);
%     if (own <=3)
%         xn = nod.xyz(1);
%         yn = nod.xyz(2);
%         %Less than 3 owners, so exterior
%         xv = [xv xn];
%         yv = [yv yn];
%         nid = [nid n];
%     end
% end
% 
% %Now need to get oriented version (will go with ccw)
% Nbdry = length(xv);
% 
% [xvc yvc I] = Orient(xv,yv);
% [xvc yvc] = poly2cw(xvc,yvc);
% 
% % %NOTE: THIS ONLY WORKS FOR CONVEX SHAPES
% 
% % xcm = sum(xv)/Nbdry;
% % ycm = sum(yv)/Nbdry;
% % xvo = xv-xcm; yvo = yv-ycm;
% % 
% % theta = atan2(yvo,xvo);
% % 
% % [ths I] = sort(theta);
% % xvc = xv(I);
% % yvc = yv(I);
% 
% %Currently assuming 1 object
% lvlDef.numObs = 1;
% obsDat.xv = xvc;
% obsDat.yv = yvc;
% obsDat.vx = zeros(1,Nbdry);
% obsDat.vy = zeros(1,Nbdry);
% obsDat.vxid = zeros(1,Nbdry);
% obsDat.vyid = zeros(1,Nbdry);
% 
% obsDat.nid = nid(I);
% 
% for n=1:Nbdry
%     nnod = obsDat.nid(n);
%     nod = Nodes{nnod};
%     vxid = nod.DOFids(1);
%     vyid = nod.DOFids(2);
%     if (vxid ~= 0)  
%         obsDat.vx(n) = v(vxid);
%     else
%         obsDat.vx(n) = 0;
%     end
%     if (vyid ~= 0)
%         obsDat.vy(n) = v(vyid);
%     else
%         obsDat.vy(n) = 0.0;
%     end
%     %Save return mapping
%     obsDat.vxid(n) = vxid;
%     obsDat.vyid(n) = vyid;
% end
% 
% lvlDef.obsDat{1} = obsDat;

%Close polygon


% Ns = 1000;
% x = linspace(0,200,Ns); y = linspace(0,200,Ns);
% [yy xx] = meshgrid(y,x);
% In = inpolygon(xx,yy,xvc,yvc);
% plot(obsDat.xv,obsDat.yv,'bo'); hold on;
% quiver(obsDat.xv,obsDat.yv,obsDat.vx,obsDat.vy,'r'); 
% axis equal
% hold off;
% 
% keyboard

%Takes collection of points (xv,yv) and returns the same points but ordered
%by distance to previous point
function [xvc yvc I] = Orient(xv,yv)
Nv = length(xv);
xvc = zeros(1,Nv);
yvc = zeros(1,Nv);
Used = false(1,Nv);
I = zeros(1,Nv);

%First point is arbitrary
xvc(1) = xv(1);
yvc(1) = yv(1);
I(1) = 1;
Used(1) = true;

for n=2:Nv
    xp = xvc(n-1);
    yp = yvc(n-1);
    dv = sqrt( (xp-xv).^2 + (yp-yv).^2 );
    dv(Used) = Inf; %Don't repick one
    [dvmin i] = min(dv);
    xvc(n) = xv(i);
    yvc(n) = yv(i);
    Used(i) = true;
    I(n) = i;
end


