
%Enforces boundary conditions on hydrodynamic variables
%For now only does doubly-periodic

function [Gas Grid Model] = EnforceBCs(Model,Grid,Gas)

Dirs = {'ibx','obx','iby','oby'};
Nd = length(Dirs);
Dvec = zeros(2,Nd);
Dvec(:,1) = [-1,0]; Dvec(:,2) = [1,0];
Dvec(:,3) = [0,-1]; Dvec(:,4) = [0,1];

%Loop over directions (n/s/e/w)
for n=1:Nd
   ComS = sprintf('bc = Model.bcs.%s;',Dirs{n});
   eval(ComS);
   dir.vec = Dvec(:,n);
   dir.str = Dirs{n};
   switch lower(bc)
       case{'wall'}
           Gas = WallBC(Model,Grid,Gas,dir);
       case{'periodic'}
           Gas = PeriodicBC(Model,Grid,Gas,dir);
       case{'reflect'}
           Gas = ReflectBC(Model,Grid,Gas,dir);
       case{'outflow'}
           Gas = OutflowBC(Model,Grid,Gas,dir);
       case{'poutflow'}
           Gas = POutflowBC(Model,Grid,Gas,dir);
       case{'inflow'}
           Gas = InflowBC(Model,Grid,Gas,dir);
       case{'pinflow'}
           Gas = PInflowBC(Model,Grid,Gas,dir);
       case{'injet'}
           Gas = InJetBC(Model,Grid,Gas,dir);
       otherwise
           disp('Unknown boundary condition');
           pause
   end
       
end

Grid = InitLvl(Model,Grid);
Gas = LvlGhost(Model,Grid,Gas);

% 
% %Calculate values at internal ghost zones (level set method) using physical
% %image points
% if (Model.lvlSet.present)
%     if (Model.lvlSet.anymobile)
%         Model = moveObject(Model,Grid,Gas);
%         
%         %Recalculate geometry
%         Grid = InitLvl(Model,Grid);
%     end
%     Gas = LvlGhost(Model,Grid,Gas);
% end


function Gas = InJetBC(Model,Grid,Gas,dir)

%Note, we assume dir = ibx
[Nx Ny] = size(Gas.D);

ng = Grid.ng;

is = Grid.is; ie = Grid.ie;
isd = Grid.isd; ied = Grid.ied;

js = Grid.js; je = Grid.je;
jsd = Grid.jsd; jed = Grid.jed;

switch (lower(dir.str))
    case{'ibx'}
        
        pulseRad = Model.Init.rad;
        pulseC = Model.Init.cent;
        Disc = Model.Init.disc;
        Shape = approxId(Grid.yc,pulseC,pulseRad,Disc);
        
        d0 = Model.Init.rho0;
        P0 = Model.Init.P0;
        DelP = Model.Init.DelP;
        vin = Model.Init.vin;
        %Rescale shape function for pressure to vary between P0 and P0*DelP
        Pshp = Shape*( P0*DelP - P0 ) + P0;
        
        for n=1:ng
            Gas.D(isd+n-1,:) = d0;
            Gas.P(isd+n-1,:) = Pshp;
            Gas.Vx(isd+n-1,:) = vin*Shape;
            Gas.Vy(isd+n-1,:) = 0.0;
        end
    otherwise
        display('P-Inflow not implemented for this direction');
        pause
end

function Gas = PInflowBC(Model,Grid,Gas,dir)
%Note, we assume dir = ibx
[Nx Ny] = size(Gas.D);

ng = Grid.ng;

is = Grid.is; ie = Grid.ie;
isd = Grid.isd; ied = Grid.ied;

js = Grid.js; je = Grid.je;
jsd = Grid.jsd; jed = Grid.jed;


switch (lower(dir.str))
    case{'ibx'}
        for n=1:ng
            Gas.D(isd+n-1,:) = Model.Init.rho0;
            Gas.P(isd+n-1,:) = Model.Init.P0*Model.Init.DelP; %P0*DelP
            Gas.Vx(isd+n-1,:) = 0.0;
            Gas.Vy(isd+n-1,:) = 0.0;
        end
        
    otherwise
        display('P-Inflow not implemented for this direction');
        pause
end

function Gas = POutflowBC(Model,Grid,Gas,dir)
%Note, we assume dir = obx
[Nx Ny] = size(Gas.D);

ng = Grid.ng;

is = Grid.is; ie = Grid.ie;
isd = Grid.isd; ied = Grid.ied;

js = Grid.js; je = Grid.je;
jsd = Grid.jsd; jed = Grid.jed;

switch (lower(dir.str))
    case{'obx'}
        for n=1:ng
            Gas.D(ie+n,:) = Model.Init.rho0;
            Gas.Vx(ie+n,:) = 0.0;
            Gas.Vy(ie+n,:) = 0;
            Gas.P(ie+n,:) = Model.Init.P0;
        end
        
    otherwise
        display('P-Outflow not implemented for this direction');
        pause
end

function Gas = InflowBC(Model,Grid,Gas,dir)
%Note, we assume dir = ibx
[Nx Ny] = size(Gas.D);

ng = Grid.ng;

is = Grid.is; ie = Grid.ie;
isd = Grid.isd; ied = Grid.ied;

js = Grid.js; je = Grid.je;
jsd = Grid.jsd; jed = Grid.jed;
vin = Model.bcs.inspeed;

switch (lower(dir.str))
    case{'ibx'}
        for n=1:ng
            %Gas.D(isd+n-1,:) = Gas.D(is,:);
            %Gas.P(isd+n-1,:) = Gas.P(is,:);
            Gas.D(isd+n-1,:) = Model.Init.rho0;
            Gas.P(isd+n-1,:) = Model.Init.P0;
            Gas.Vx(isd+n-1,:) = vin;
            Gas.Vy(isd+n-1,:) = 0.0;
        end
        
    otherwise
        display('Inflow not implemented for this direction');
        pause
end

function Gas = OutflowBC(Model,Grid,Gas,dir)
%Note, we assume dir = obx
[Nx Ny] = size(Gas.D);

ng = Grid.ng;

is = Grid.is; ie = Grid.ie;
isd = Grid.isd; ied = Grid.ied;

js = Grid.js; je = Grid.je;
jsd = Grid.jsd; jed = Grid.jed;

switch (lower(dir.str))
    case{'ibx'}
        for n=1:ng
            Gas.D(is-n,:) = Gas.D(is,:);
            Gas.Vx(is-n,:) = max( Gas.Vx(is,:), 0); %Diode
            Gas.Vy(is-n,:) = 0.0;
            Gas.P(is-n,:) = Gas.P(is,:);
        end
    case{'obx'}
        for n=1:ng
            Gas.D(ie+n,:) = Gas.D(ie,:);
            Gas.Vx(ie+n,:) = max(Gas.Vx(ie,:),0); %Diode outflow
            Gas.Vy(ie+n,:) = 0;
            Gas.P(ie+n,:) = Gas.P(ie,:);
        end
    case{'iby'}
        for n=1:ng
            Gas.D(:,js-n) = Gas.D(:,js);
            Gas.Vx(:,js-n) = 0.0;
            Gas.Vy(:,js-n) = min( Gas.Vy(:,js), 0);
            Gas.P(:,js-n) = Gas.P(:,js);
        end
    case{'oby'}
        for n=1:ng
        	Gas.D(:,je+n) = Gas.D(:,je);
            Gas.Vx(:,je+n) = 0.0;
            Gas.Vy(:,je+n) = max( Gas.Vy(:,je), 0);
            Gas.P(:,je+n) = Gas.P(:,je);
        end
    otherwise
        display('Outflow not implemented for this direction');
        pause
end

function Gas = WallBC(Model,Grid,Gas,dir)
Vars = {'D','Vx','Vy','P'};

fl = [1,0,0,1]; %Sign, reflect density/pressure zero out velocity

%Loop over hydro variables
for v=1:length(Vars)
    var = Vars{v};
    ComS = sprintf('Gas.%s = ReflectVar(Gas.%s,Grid,dir,fl(v));',var,var);
    eval(ComS);
end

function Gas = ReflectBC(Model,Grid,Gas,dir)

Vars = {'D','Vx','Vy','P'};
fl = [1,1,1,1]; %Sign, put negative in component to be reflected

switch (lower(dir.str))
    case{'ibx','obx'}
        fl(2) = -1;
    case{'iby','oby'};
        fl(3) = -1;
end

%Loop over hydro variables
for v=1:length(Vars)
    var = Vars{v};
    ComS = sprintf('Gas.%s = ReflectVar(Gas.%s,Grid,dir,fl(v));',var,var);
    eval(ComS);
end

function Qbc = ReflectVar(Q,Grid,dir,sgn)

[Nx Ny] = size(Q);
Qbc = Q;
ng = Grid.ng;
is = Grid.is; ie = Grid.ie;
js = Grid.js; je = Grid.je;

switch lower(dir.str)
    case{'ibx'}
        for n=1:ng
            Qbc(is-n,:) = sgn*Q(is+n-1,:);
        end
    case{'obx'}
        for n=1:ng
            Qbc(ie+n,:) = sgn*Q(ie-n+1,:);
        end
    case{'iby'}
        for n=1:ng
            Qbc(:,js-n) = sgn*Q(:,js+n-1);
        end
    case{'oby'}
        for n=1:ng
            Qbc(:,je+n) = sgn*Q(:,je-n+1);
        end
end

function Gas = PeriodicBC(Model,Grid,Gas,dir)

Vars = {'D','Vx','Vy','P'};
%Loop over hydro variables
for v=1:length(Vars)
   var = Vars{v};
   ComS = sprintf('Gas.%s = PeriodicVar(Gas.%s,Grid,dir);',var,var);
   eval(ComS);
end


function Qbc = PeriodicVar(Q,Grid,dir)
[Nx Ny] = size(Q);
Qbc = Q;
ng = Grid.ng;
is = Grid.is; ie = Grid.ie;
js = Grid.js; je = Grid.je;

switch lower(dir.str)
    case{'ibx'}
        for n=1:ng
            Qbc(is-n,:) = Q(ie-n+1,:);
        end
    case{'obx'}
        for n=1:ng
            Qbc(ie+n,:) = Q(is+n-1,:);
        end
    case{'iby'}
        for n=1:ng
            Qbc(:,js-n) = Q(:,je-n+1);
        end
    case{'oby'}
        for n=1:ng
            Qbc(:,je+n) = Q(:,js+n-1);
        end
end


%Enforces doubly periodic on the variable Q
function Qbc = AllPeriodic(Q,Grid)
[Nx Ny] = size(Q);
Qbc = Q;

ng = Grid.ng;
is = Grid.is; ie = Grid.ie;
js = Grid.js; je = Grid.je;

for n=1:ng
    Qbc(is-n,:) = Q(ie-n+1,:); %Inner i
    Qbc(ie+n,:) = Q(is+n-1,:); %Outer i
    Qbc(:,js-n) = Q(:,je-n+1); %Inner j
    Qbc(:,je+n) = Q(:,js+n-1); %Outer 
end
 