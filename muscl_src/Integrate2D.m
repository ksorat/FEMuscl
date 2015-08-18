function [Gas Model] = Integrate2D(Model,Grid,Gas)


%Integrate from Grid.t->Grid.t+Grid.dt

dt = Grid.dt;
halfdt = Grid.dt/2;

%Step 1 - Advance from t->t+0.5*dt

%Step 1a - Calculate L/R states at each interface using piecewise-constant
%reconstruction in primitive variables

%Start by calculating directed states using PCM
[nW eW sW wW] = ReconStates(Grid,Gas,'pcm');

%Now convert directed states to interface LR's (w/ x->y in Y-interfaces to
%trick the Riemann solver)

[Wxl Wxr Wyl Wyr] = DirstateToIntLR(Grid,Gas,nW,eW,sW,wW);

%Step 1b - Turn L/R into fluxes
%Construct Fx / Fy, fluxes in x and y
%Fx(q,i,j) = flux @ interface i-1/2,j of conserved quantity q
%Fy(q,i,j) = flux @ interface i,j-1/2 of consderved quantity q
%Note, |Fx| = [Nv, Nx+1,Ny]
%Also note, calculating Fx/Fy are independent and a good source of
%parallelism

[Fx Fxl] = RiemannFlux(Wxl,Wxr,Model,Model.solver,Grid);
[Fxy Fxyl] = RiemannFlux(Wyl,Wyr,Model,Model.solver,Grid); %Needs to be untwisted, 2<->3

Fy = UntwistFlux(Fxy);
Fyl = UntwistFlux(Fxyl);

%Step 1c - Evolve to half-timestep states
[D_hf Vx_hf Vy_hf P_hf] = FluxAdvance(Grid,0.5*Grid.dt,Gas.D,Gas.Vx,Gas.Vy,Gas.P,Fx,Fy,Fxl,Fyl,Model);

%Step 2 - Calculate dircted states using PLM w/ half-timestep states
%Step 2a - Package half-timestep states into Gas_hf
Gas_hf = Gas;
Gas_hf.D = D_hf;
Gas_hf.Vx = Vx_hf; Gas_hf.Vy = Vy_hf;
Gas_hf.P = P_hf;
%Step 2b - Calculate directed states
[nW eW sW wW] = ReconStates(Grid,Gas_hf,Model.recon);
%Step 2c - Convert directed states to interface LR's (w/ x<->y twist)
[Wxl Wxr Wyl Wyr] = DirstateToIntLR(Grid,Gas,nW,eW,sW,wW);

%Step 3 - Calculate fluxes
[Fx Fxl] = RiemannFlux(Wxl,Wxr,Model,Model.solver,Grid);
[Fxy Fxyl] = RiemannFlux(Wyl,Wyr,Model,Model.solver,Grid); %Needs to be untwisted, 2<->3

Fy = UntwistFlux(Fxy);
Fyl = UntwistFlux(Fxyl);

%Step 4 - Evolve from initial state (in Gas) -> final state over full
%timestep
%Add D_hf for time-centered density included in external forces

[D Vx Vy P] = FluxAdvance(Grid,Grid.dt,Gas.D,Gas.Vx,Gas.Vy,Gas.P,Fx,Fy,Fxl,Fyl,Model,D_hf);
Gas.D = D;
Gas.Vx = Vx; Gas.Vy = Vy;
Gas.P = P;

%If necessary, calculate forces exerted by fluid on object
if (Model.lvlSet.present)
    %Note, should maybe use the corrected fluxes when appropriate
    [outForce Model] = calcOutforce(Grid,Model,Fx,Fy);
end


%Uses flux to advance W_i -> W_o over a time dt
function [Do Vxo Vyo Po] = FluxAdvance(Grid,dt,Di,Vxi,Vyi,Pi,Fx,Fy,Fxl,Fyl,Model,varargin)

global SMALL_NUM;
fac = 1.2; %This is fairly ad hoc

if ( length(Fxl) > 0 )
    FluxCorrect = true;
else
    FluxCorrect = false;
end

[Nx,Ny] = size(Di);
%Note, |Fx| = [Nv,Nx+1,Ny]
%|Fy| = [Nv,Nx,Ny+1]

%Choose which density to use for external forces, if any

if (length(varargin) > 0)
    Dsrc = varargin{1};
else
    Dsrc = Di;
end

%Convert to conserved W_i -> U_i 
[Di Mxi Myi Ei] = Prim2Con(Di,Vxi,Vyi,Pi,Model);





if (Model.lvlSet.present)
    %Zero out interfaces connecting to "obj" cells (NOT Ghosts)
    Fx_obj = Grid.lvlSet.Fx_obj;
    Fy_obj = Grid.lvlSet.Fy_obj;
    Fx(:,Fx_obj) = 0.0;
    Fy(:,Fy_obj) = 0.0;
    if (FluxCorrect)
        Fxl(:,Fx_obj) = 0.0;
        Fyl(:,Fy_obj) = 0.0;
    end
end

%Advance in time, U_i -> U_o
dtox = dt/Grid.dx;
dtoy = dt/Grid.dy;

Do = Di(:,:) + squeeze(dtox*( Fx(1,1:Nx,:) - Fx(1,2:Nx+1,:) ) + ...
    + dtoy*( Fy(1,:,1:Ny) - Fy(1,:,2:Ny+1) ));

Mxo = Mxi(:,:) + squeeze(dtox*( Fx(2,1:Nx,:) - Fx(2,2:Nx+1,:) ) + ...
    + dtoy*( Fy(2,:,1:Ny) - Fy(2,:,2:Ny+1) ));

Myo = Myi(:,:) + squeeze(dtox*( Fx(3,1:Nx,:) - Fx(3,2:Nx+1,:) ) + ...
    + dtoy*( Fy(3,:,1:Ny) - Fy(3,:,2:Ny+1) ));

Eo = Ei(:,:) + squeeze(dtox*( Fx(4,1:Nx,:) - Fx(4,2:Nx+1,:) ) + ...
    + dtoy*( Fy(4,:,1:Ny) - Fy(4,:,2:Ny+1) ));

if ( FluxCorrect )
    %We've calculated high and low order fluxes, so correct if necessary
    %Need pressure
    %[Do Vxo Vyo Po] = Con2Prim(Do,Mxo,Myo,Eo,Model);
    [~, ~, ~, Po] = Con2Prim(Do,Mxo,Myo,Eo,Model);
    Ind = (Do < fac*SMALL_NUM) | (Po < fac*SMALL_NUM);
    Indn = find(Ind); NumC = sum(Ind(:));
    
    if (NumC > 0)
        fprintf('\tHLLC->HLLE @ %d cell(s).\n', NumC);
        [jj ii] = meshgrid(1:Ny,1:Nx);
        for nc=length(Indn)
            n = Indn(nc);
            i = ii(n); j = jj(n);
            %Redo update with HLLE fluxes
            
            Do(i,j) = Di(i,j) + dtox*( Fxl(1,i,j) - Fxl(1,i+1,j) ) + ...
                + dtoy*( Fyl(1,i,j) - Fyl(1,i,j+1) );
            Mxo(i,j) = Mxi(i,j) + dtox*( Fxl(2,i,j) - Fxl(2,i+1,j) ) + ...
                + dtoy*( Fyl(2,i,j) - Fyl(2,i,j+1) );
            Myo(i,j) = Myi(i,j) + dtox*( Fxl(3,i,j) - Fxl(3,i+1,j) ) + ...
                + dtoy*( Fyl(3,i,j) - Fyl(3,i,j+1) );
            Po(i,j) = Pi(i,j) + dtox*( Fxl(4,i,j) - Fxl(4,i+1,j) ) + ...
                + dtoy*( Fyl(4,i,j) - Fyl(4,i,j+1) );            
        end
    end
    
end
%Add forces if necessary
if isfield(Model,'force')
    %There is a force, incorporate it
    %Calculate forces @
    %fxc/fyc, both components at center
    %fw/fe, fx at west and east
    %fn/fs, fy at north and south
    
    switch lower(Model.force.type)
        case{'simpgrav'}
            %Simple gravity, f = (0,-g)
            fxc = zeros(Nx,Ny);
            fyc = -1*Model.force.g*ones(Nx,Ny);
            
            %Apply force
            %[Mxo Myo Eo] = ApplyForce(dt,Grid,Mxo,Myo,Eo,Dsrc,Fx,Fy,fxc,fyc,fn,fe,fs,fw);
            [Mxo Myo Eo] = ApplyForce(dt,Grid,Model,Do,Mxo,Myo,Eo,Dsrc,fxc,fyc);
        otherwise
            %Do the general form
    end
    
end

%Convert back to primitive and return, U_o -> W_o
[Do Vxo Vyo Po] = Con2Prim(Do,Mxo,Myo,Eo,Model);



%Simpler force application.  Take input states and D-src and only centered
%forces
function [Mxo Myo Eo] = ApplyForce(dt,Grid,Model,Di,Mxi,Myi,Ei,Dsrc,fxc,fyc)

[Di Vxi Vyi Pi] = Con2Prim(Di,Mxi,Myi,Ei,Model);
Mxo = Mxi + dt*fxc.*Dsrc;
Myo = Myi + dt*fyc.*Dsrc;

[Do Mxo Myo Eo] = Prim2Con(Di,Mxo./Di,Myo./Di,Pi,Model);

%Apply effect of forces to momenta/energy.  Take in input states and output
%updated states.  Use dsrc as density associated with the time-centered
%force
%Note, don't take dt from Grid struct b/c it may be only half step update
function [Mxo Myo Eo] = ApplyForce_Old(dt,Grid,Mxi,Myi,Ei,Dsrc,Fluxx,Fluxy,fxc,fyc,fn,fe,fs,fw)

%Apply forces to momenta
Mxo(:,:) = Mxi(:,:) + dt*fxc.*Dsrc;
Myo(:,:) = Myi(:,:) + dt*fyc.*Dsrc;

Nx = Grid.Nx;
Ny = Grid.Ny;

%Apply force to energy
mFw = squeeze(Fluxx(1,1:Nx,:)); %Mass flux from west
mFe = squeeze(Fluxx(1,2:Nx+1,:)); %Mass flux from east
mFs = squeeze(Fluxy(1,:,1:Ny)); %Mass flux from south
mFn = squeeze(Fluxy(1,:,2:Ny+1)); %Mass flux from north

%Eo(:,:) = Ei(:,:) + 0.5*dt*(mFw.*fw + mFe.*fe) + 0.5*dt*(mFs.*fs + mFn.*fn);
%Note, changing to simpler calculation but possibly not as accurate
Eo(:,:) = Ei(:,:) + dt*( Mxi(:,:).*fxc + Myi(:,:).*fyc );
    
%--------
%Converts directed states within each cell N/S/E/W to
%LR states at each interface
function [Wxl Wxr Wyl Wyr] = DirstateToIntLR(Grid,Gas,nW,eW,sW,wW)

%Create holders
[Nx,Ny] = size(Gas.D);

%X-interfaces = [Nx+1,Ny]
Wxl = zeros(4,Nx+1,Ny);
Wxr = Wxl;

%Y-Interfaces = [Nx,Ny+1]
Wyl = zeros(4,Nx,Ny+1);
Wyr = Wyl;

%Now map into X-Int LR
%Wxl(i,j) = L state @ i-1/2,j
Wxl(:,2:Nx+1,:) = eW(:,:,:);
Wxr(:,1:Nx,:)   = wW(:,:,:);

%Fill in gaps
Wxl(:,1,:)    = Wxr(:,1,:);
Wxr(:,Nx+1,:) = Wxl(:,Nx+1,:);

%Now map into Y-Int LR (untwisted)
%Wyl(i,j) = S state @ i,j-1/2
Wyl(:,:,2:Ny+1) = nW(:,:,:);
Wyr(:,:,1:Ny)   = sW(:,:,:);

%Fill in gaps
Wyl(:,:,1)    = Wyr(:,:,1);
Wyr(:,:,Ny+1) = Wyl(:,:,Ny+1);

%Now twist y, ie Vx<->Vy to trick Riemann solver
Tmp = Wyl(2,:,:);
Wyl(2,:,:) = Wyl(3,:,:); Wyl(3,:,:) = Tmp;

Tmp = Wyr(2,:,:);
Wyr(2,:,:) = Wyr(3,:,:); Wyr(3,:,:) = Tmp;


function[Vxo Vyo] = CorrectVel(Grid,Vxo,Vyo)
Nx = Grid.Nx; Ny = Grid.Ny;
sol = Grid.solid;
In = Grid.solid.In;

%Start by setting velocity of solid cells, eventually talk to solid model
Vxo(In) = 0.0;
Vyo(In) = 0.0;
bf3d = sol.BdryFlu;
bf1d = sol.bf1d; %1d vector of fluid boundary cells
%Create 2d i/j matrices
[ii,jj] = meshgrid(1:Nx,1:Ny);
ii = ii'; jj = jj';
dirx = [1 0 -1 0];
diry = [0 1 0 -1];

for n=1:length(bf1d);
    elt = bf1d(n);
    i = ii(elt);
    j = jj(elt);
    N = 0; vx_avg = 0; vy_avg = 0;
    
    for d=1:length(dirx)
        %Loop through cardinal neighbors and decide who's relevant
        dx = dirx(d);
        dy = diry(d);
        if ~sol.BdryFlu(i+dx,j+dy)
            %This neighbor is relevant, ie not another boundary cell
            N = N+1;
            vx_avg = vx_avg + Vxo(i+dx,j+dy);
            vy_avg = vy_avg + Vyo(i+dx,j+dy);
        end
        Vxo(i,j) = (1/N)*vx_avg;
        Vyo(i,j) = (1/N)*vy_avg;
    end
            
end
%keyboard

    
function Fy = UntwistFlux(Fxy)

if ( length(Fxy) > 0 )
    Fy = Fxy;
    Fy(2,:,:) = Fxy(3,:,:);
    Fy(3,:,:) = Fxy(2,:,:);
else
    Fy = Fxy;
end