
function Gas = InitGas(Model,Grid)

%Parse problem
switch lower(Model.Init.problem)
    case{'blast'}
        Gas = InitBlast(Model,Grid);
    case{'khi'}
        Gas = InitKHI(Model,Grid);
    case{'rti'}
        Gas = InitRTI(Model,Grid);
    case{'imp'}
        Gas = InitImplosion(Model,Grid);
    case{'blastobs'}
        Gas = InitBlastObs(Model,Grid);
    case{'flow'}
        Gas = InitFlow(Model,Grid);
    otherwise
        disp('Unknown problem type');
        pause
end


function Gas = InitBlast(Model,Grid)

Gas = BuildGas(Model,Grid);
[xx,yy] = meshgrid(Grid.xc,Grid.yc);
xx = xx'; yy = yy';

if isfield(Model.Init,'cent')
    x0 = Model.Init.cent(1);
    y0 = Model.Init.cent(2);
else
    x0 = 0;
    y0 = 0;
end

rad = sqrt( (xx-x0).^2 + (yy-y0).^2);

Gas.D(:,:) = Model.Init.rho0;
Gas.Vx(:,:) = 0.0;
Gas.Vy(:,:) = 0.0;
Gas.P(:,:) = Model.Init.P0;

Inside = (rad < Model.Init.rad);
Gas.P(Inside) = Model.Init.P0*Model.Init.DelP;

function Gas = InitRTI(Model,Grid)

Gas = BuildGas(Model,Grid);
[xx,yy] = meshgrid(Grid.xc,Grid.yc);
xx = xx'; yy = yy';

g = Model.force.g;

%Set density above zero equal to 2, and below equal to 1
Ind = (yy < 0);
Gas.D(Ind) = 1.0;
Gas.D(~Ind) = 2.0;

Gas.P = 2.5 - g*Gas.D.*yy;

%Perturb Vy
Lx = Grid.xe-Grid.xs;
Ly = Grid.ye - Grid.ys;
a = Model.Init.amp;

Gas.Vy = a* ( 1 + cos(2*pi*xx/Lx) ) .* ( 1 + cos(2*pi*yy/Ly) ) / 4;

function Gas = InitKHI(Model,Grid)

Gas = BuildGas(Model,Grid);
[xx,yy] = meshgrid(Grid.xc,Grid.yc);
xx = xx'; yy = yy';

Gas.D(:,:) = Model.Init.rho0;
Gas.P(:,:) = Model.Init.P0;
Gas.Vy(:,:) = 0.0;


yqts = linspace(Grid.yi(Grid.js), Grid.yi(Grid.je+1), 5);

Mid = ( yy <= yqts(4) ) & ( yy >= yqts(2) );
Out = ( yy > yqts(4) ) | (yy <= yqts(2) );

Gas.D(Mid) = Gas.D(Mid)*Model.Init.DelD; %Increase density in middle
Gas.Vx(Mid) = Model.Init.Vx0;
Gas.Vx(Out) = -1*Model.Init.Vx0;

[Nx Ny] = size(xx);
r = 2*( rand(Nx,Ny) - 0.5 ); %Random numbers between -1/1
rMid = r.*Mid; %Restrict fluctuations to middle region

%Random perturbations
%Gas.P = (1+Model.Init.amp*rMid).*Gas.P;
%Gas.Vy = Model.Init.amp*Model.Init.Vx0*rMid;

%Single mode
Numwv = 4;
DelX = Grid.xi(Grid.ie+1)-Grid.xi(Grid.is);
xc = Grid.xc;
vy = Model.Init.amp*Model.Init.Vx0*sin(Numwv*pi*xc/DelX);
Vyall = repmat(vy',[1,Grid.Ny]);
Gas.Vy = Mid.*Vyall;

function Gas = InitImplosion(Model,Grid)

Gas = BuildGas(Model,Grid);
[xx,yy] = meshgrid(Grid.xc,Grid.yc);
xx = xx'; yy = yy';

%Vx = Vy = 0
Ind = (xx + yy ) > 0.5;
Gas.D(Ind) = 1;
Gas.P(Ind) = 1;

Gas.D(~Ind) = 0.125;
Gas.P(~Ind) = 0.14;


function Gas = InitFlow(Model,Grid)
Gas = BuildGas(Model,Grid);
Gas.D(:,:) = Model.Init.rho0;
Gas.Vx(:,:) = 0.0;
Gas.Vy(:,:) = 0.0;
Gas.P(:,:) = Model.Init.P0;

function Gas = BuildGas(Model,Grid)

if ~isfield(Model.Init,'Gam')
    Gas.Gam = 5/3;
else
    Gas.Gam = Model.Init.Gam;
end

Nx = Grid.Nx; Ny = Grid.Ny;

%Create primitive variables
Gas.D = zeros(Nx,Ny);
Gas.Vx = zeros(Nx,Ny);
Gas.Vy = zeros(Nx,Ny);
Gas.P = zeros(Nx,Ny);

Gas.Nv = 4;

        