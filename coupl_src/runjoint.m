
function [Grid Gas Nodes Elements]  = runjoint(Model)

clear global

%Initialize solid
global Parts Nodes Elements Materials ElForms Damage Lengths Speeds GDOF GDOFall u v Fext Fmax M dim a t;
global filename DoDamage;
DoDamage = Model.Init.DoDamage;
filename = Model.Init.filename;

tsScale = 0.9;
Solver = explicitSolverTF();
InitSolid();

lvlDef = Solid2Gas(Model);
%keyboard
Model.Init.lvlDef = lvlDef; 
Model.lvlSet.present = true;
Model.lvlSet.anymobile = true;

Model = FixMod(Model);

%Construct grid
%xs/xe = x-bounds, is/ie = physical indices, isd/ied = physical+ghost
Grid = BuildGrid(Model);

%Initialize gas
Gas = InitGas(Model,Grid);


%Enforce Gas BC's on initial setup
[Gas Grid Model] = EnforceBCs(Model,Grid,Gas);

%Calculate initial Gas timestep
Grid.dt = CalcDT(Model,Grid,Gas);
Tfin = Model.Tfin;


%Enter timeloop
while (Grid.t<Tfin)
    %Evolve from t->t+dt

    %Evolve gas subject to solid
    [Gas Model] = Integrate2D(Model,Grid,Gas);
    
    %Calculate force exerted on solid by gas
    Fext = Gas2Solid(Model.Init.lvlDef);
    
    %Evolve solid with external forces from gas
    Solver.run(tsScale,Grid.dt);
    

    %Bring in new solid positions/velocities
    lvlDef = Solid2Gas(Model);
    Model.Init.lvlDef = lvlDef; 
    
    %keyboard
    %Print diagnostics if necessary
    if (mod(Grid.Ts,Grid.tsDiag) == 0)    
        printDiag(Grid,Gas);
        if (Model.Pic.view)
            makeFig(Model,Grid,Gas);
        end
            
    end
    
    %Enforce BCs
    [Gas Grid Model] = EnforceBCs(Model,Grid,Gas);
    
    %Update time
    Grid.t = Grid.t+Grid.dt;
    Grid.Ts = Grid.Ts+1;
    
    
    %Calc new timestep
    Grid.dt = CalcDT(Model,Grid,Gas);
end


    
function InitSolid()


global Parts Nodes Elements Materials ElForms Damage Lengths Speeds GDOF GDOFall u v Fext Fmax M dim a t;
global filename

    dim = 2; % model dimensions (either 2 or 3) (INT)
    GDOF = 0; % total global degrees of freedom. Could work well as a static INT OR LONG (this can get to be a super big number in some models) for the node object
    GDOFall = 0; % total global degrees of freedom for the stress projection. Could work well as a static INT OR LONG (this can get to be a super big number in some models) for the node object


    % Setup Part, Material, and Element Profiles

    Parts = cell(1,1); % open up a Parts array
    Materials = cell(1,1); % open up a Materials array
    ElForms = cell(1,1); % open up an Element Formulation array
    Damage = {}; % open up an Element Formulation array
    Speeds = zeros(length(Materials),1); % initialize the sound speeds array
    ElForms{1} = Q4; % Creating the standard Q4 element
    Damage{1} = maxPrincipalDamage(.0008);
    %Damage{1} = maxPrincipalDamage(.008);
    Materials{1} = LinearElasticTL_planeStress(1,1.1e-6,.05,.45); % create material model (ID,density,modulus,poisson)
    %Materials{1} = LinearElasticTL_planeStress(1,1.1e-6,.15,.3); % create material model (ID,density,modulus,poisson)
    Parts{1} = Part(Materials{1},ElForms{1}); % create part object
    Parts{1}.damage = 1;
    for i = 1:length(Speeds)

        if length(Speeds) == 1

            Speeds = Materials{1}.c; % c is a float held within each Materials object
            break

        end

        Speeds(i) = Materials{i}.c;

    end % populate sound speed array for timestep calcs
    
        % Load Input Deck

    [~,~,~,~,~] = inputdeckreader(filename); % read LS-DYNA input deck (old, lots of uselessness inside of it, but it works for now)

    % Initialize Time Integration Variables

    M = diag(M); % turn M matrix into an M vector for 2x speedup of a=F/M calculation (M matrix is diagonal anyway due to the way we evaluated it) (a GDOF x 1 VECTOR of FLOATS)
    u = zeros(GDOF,1); % displacement (a GDOF x 1 VECTOR of FLOATS) 
    v = zeros(GDOF,1); % velocity (a GDOF x 1 VECTOR of FLOATS)
    Fext = zeros(GDOF,1); % external force (a GDOF x 1 VECTOR of FLOATS)
    a = zeros(GDOF,1); % acceleration (a GDOF x 1 VECTOR of FLOATS)


function printDiag(Grid,Gas)
%Print diagnostics at given cadence

fprintf('\tTime = %3.3f : Step = %4d, dt = %3.2e\n',Grid.t,Grid.Ts,Grid.dt);


function Grid = BuildGrid(Model)

Bds = Model.Bds;
Grid.xs = Bds(1); Grid.xe = Bds(2); Grid.ys = Bds(3); Grid.ye = Bds(4);

if isfield(Model,'Res')
    dx = Model.Res(1); dy = Model.Res(2);
    Model.Nvec = [ ceil( (Grid.xe-Grid.xs)/dx ) ceil( (Grid.ye-Grid.ys)/dy ) ];
end
Nx = Model.Nvec(1);
Ny = Model.Nvec(2);
dx = (Grid.xe-Grid.xs)/Nx;
dy = (Grid.ye-Grid.ys)/Ny;

Grid.Nxp = Nx; Grid.Nyp = Ny;
Grid.dx = dx; Grid.dy = dy;
Grid.ng = Model.ng;

%Construct dimensions
%xc = cell-centered
%xi = interface values
% is/ie: physical cells
% isd/ied: all cells


[Grid.xi Grid.xc Grid.is Grid.ie Grid.isd Grid.ied] = conDimNg(Grid.xs,Grid.xe,Grid.Nxp,Grid.ng);
Grid.Nx = length(Grid.xc);

[Grid.yi Grid.yc Grid.js Grid.je Grid.jsd Grid.jed] = conDimNg(Grid.ys,Grid.ye,Grid.Nyp,Grid.ng);
Grid.Ny = length(Grid.yc);

Grid.dt = 0;
Grid.t = 0;
Grid.Ts = 0;
Grid.C0 = 0.45;

Grid.tsDiag = Model.tsDiag;
fprintf('Building grid ...\n');
fprintf('\t(Dx,Dy) = (%f,%f)\n', dx,dy);
fprintf('\tAspect Ratio = %f\n', dx/dy);

function [xi xc is ie isd ied] = conDimNg(xs,xe,Nxp,ng)

Nx = Nxp + 2*ng;

xip = linspace(xs,xe,Nxp+1);
xcp = 0.5*(xip(1:end-1)+xip(2:end));
dx = xcp(2)-xcp(1);

isd = 1;
ied = Nx;
is = ng+1;
ie = ied-ng;

xc = zeros(1,Nx);
xi = zeros(1,Nx+1);

xc(is:ie) = xcp;

for n=1:ng
    xc(ie+n) = xc(ie+n-1) + dx;
    xc(is-n) = xc(is-n+1) - dx;
end
xi(1:end-1) = xc-dx/2;
xi(end) = xc(end)+dx/2;

function Model = FixMod(Model)

if ~isfield(Model,'ng')
    Model.ng = 3;
end

if ~isfield(Model,'solver')
    Model.solver = 'hllc';
end

if ~isfield(Model,'recon')
    Model.recon = 'ppm';
end

if ~isfield(Model.Init,'Gam')
    Model.Init.Gam = 5/3;
end

if ~isfield(Model,'Pic')
    Model.Pic.view = false;
    Model.Pic.dovid = false;
else
    if ~isfield(Model.Pic,'dovid')
        Model.Pic.dovid = false;
    end
end

%Handle video directory default
if (Model.Pic.dovid) & ~isfield(Model.Pic,'vid_dir')
        Model.Pic.vid_dir = 'Vids/scratch';
end

if ~isfield(Model,'tsDiag')
    Model.tsDiag = 10;
end
% 
% 
% if isfield(Model.Init,'lvlDef')
%     Model.lvlSet.present = true;
%     lvlDef = Model.Init.lvlDef;
%     anymobile = false;
%     allpoly = true;
%     
%     for n=1:lvlDef.numObs
%         obsDat = lvlDef.obsDat{n};
%         if ~isfield(obsDat,'mobile')
%             obsDat.mobile = false;
%         end
%         if (obsDat.mobile)
%             anymobile = true;
%         end
%         switch lower(obsDat.type)
%             case{'circle'}
%                 %Can add this later
%                 allpoly = false;
%             case{'poly'}
%                 obsDat = fixPoly(obsDat);
%         end
%         lvlDef.obsDat{n} = obsDat;
%     end
%     Model.lvlSet.anymobile = anymobile;
%     Model.lvlSet.allpoly = allpoly;
%     Model.Init.lvlDef = lvlDef;
% else
%     Model.lvlSet.anymobile = false;
%     Model.lvlSet.present = false;
%     Model.lvlSet.allpoly = false;
% end



%Set unset BCs
%All periodic if nothing set
if ~isfield(Model,'bcs')
    Model.bcs.ibx = 'periodic';
    Model.bcs.obx = 'periodic';
    Model.bcs.iby = 'periodic';
    Model.bcs.oby = 'periodic';
end

global SMALL_NUM;
global DEBUG;
global Pmin;
global Dmin;
global Nfig;


SMALL_NUM = 1.0e-8;
DEBUG = true;
Pmin = 1.0e-8;
Dmin = 1.0e-8;
Nfig = 0;


function obsDat = fixPoly(obsDat)

%Close polygon if necessary
if ((obsDat.xv(1) ~= obsDat.xv(end)) || (obsDat.yv(1) ~= obsDat.yv(end))) 
    obsDat.xv = [obsDat.xv ; obsDat.xv(1)];
    obsDat.yv = [obsDat.yv ; obsDat.yv(1)];
end
Nv = length(obsDat.xv);


%Create/finalize velocity data if necessary
if ( ~isfield(obsDat,'vx') | ~isfield(obsDat,'vy') )
    obsDat.vx = zeros(1,Nv);
    obsDat.vy = zeros(1,Nv);
    obsDat.Vcm_x = 0;
    obsDat.Vcm_y = 0;
end

if isfield(obsDat,'v0') %Initialize constant velocity
    obsDat.vx(:) = obsDat.v0(1);
    obsDat.vy(:) = obsDat.v0(2);
    obsDat.Vcm_x = obsDat.v0(1);
    obsDat.Vcm_y = obsDat.v0(2);
end
    
%Zero out initial acceleration
ax = zeros(1,Nv);
obsDat.ax = ax; obsDat.ay = ax;

%Create rigid body goodies
obsDat.omega = 0;
obsDat.Fx = zeros(1,Nv);
obsDat.Fy = zeros(1,Nv);
    