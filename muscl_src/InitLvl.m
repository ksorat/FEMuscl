function Grid = InitLvl(Model,Grid)
%Take data from Model.Init.lvlDef to generate obstructions via the level
%set method

%This has been changed from the previous method.  Now is based on polygons

%lvlDef contains
%numObs = number of obstructions
%obsDat{i} = data for each obstruction
%obsDat{i}.type = circle/poly
%circle type
%       .radius
%       .center = xc/yc
%poly type
%       .xv = x vertices
%       .yv = y vertices

%For all types
%obsDat{i}.mobile = true/false, is this object moving
%obsDat{i}.veldat = data about velocity of this object

%Data to be returned into Grid.lvlSet
%sd = signed distance at each point
%ghost = boolean, true at each point that is an interior ghost
%obj = boolean, true at each point that is an interior solid
%fluid = boolean, true at each point that is a fluid cell (unnecessary)
%ghost1d = 1d indices of ghost cells
%ghost_sd = 1d signed distances of ghost cells
%gi/gj = 1d i/j coordinates of ghost cells
%nx/ny = normal vector toward boundary for each interior ghost
%vx/vy = velocity vector for each interior ghost
%        THIS IS A LIE, we are assuming rigid motion for now
%ds = grid spacing to use, max{dx,dy}
%dip = probe length

global DEBUG;

lvlDef = Model.Init.lvlDef;
Nx = Grid.Nx; Ny = Grid.Ny;
numObs =lvlDef.numObs;

lvlSet.ds = max(Grid.dx,Grid.dy);
lvlSet.ds_min = min(Grid.dx,Grid.dy);


sd = inf(Nx,Ny);
Vx = zeros(Nx,Ny);
Vy = zeros(Nx,Ny);
Nvecx = zeros(Nx,Ny);
Nvecy = zeros(Nx,Ny);

%Create embiggen'ed arrays
[yy xx] = meshgrid(Grid.yc,Grid.xc);
[jj ii] = meshgrid(1:Ny,1:Nx);

for n=1:numObs
    obsDat = lvlDef.obsDat{n};
    %In this loop update sd, nx/ny, and vx/vy
    
    xv = obsDat.xv;
    yv = obsDat.yv;
    aLvl = lvlPoly(obsDat,Grid);
  
   

    
    %aLvl has been created
    Ind = (aLvl.sd < sd); %Where is the object the closest thing
    sd(Ind) = aLvl.sd(Ind);
    Nvecx(Ind) = aLvl.nx(Ind);
    Nvecy(Ind) = aLvl.ny(Ind);
    
    Vx(Ind) = aLvl.vx(Ind);
    Vy(Ind) = aLvl.vy(Ind);
    
end

%Calculate things derived from sd/V/N
fluid = (sd > 0);
obj = ( sd < -2*sqrt(2)*lvlSet.ds );
ghost = (~fluid) & (~obj);

%Create 1d arrays w/ extra info about the interior ghosts
ghost1d = find(ghost); lvlSet.ghost1d = ghost1d;
lvlSet.gi = ii(ghost1d);
lvlSet.gj = jj(ghost1d);
lvlSet.ng = length(ghost1d);
lvlSet.ghost_sd = sd(ghost1d);
lvlSet.gVx = Vx(ghost1d); lvlSet.gVy = Vy(ghost1d);
lvlSet.gNx = Nvecx(ghost1d); lvlSet.gNy = Nvecy(ghost1d);

lvlSet.dip = 1.75*lvlSet.ds;
lvlSet.sd = sd;
lvlSet.fluid = fluid;
lvlSet.ghost = ghost;
lvlSet.obj = obj;

lvlSet = calcGeom(Grid,lvlSet);

Grid.lvlSet = lvlSet;

if (DEBUG & ( (Grid.t+Grid.dt) <eps) & (numObs == 1) ) %Only do once
    plot(xv,yv,'r-o'); hold on;
    quiver(xx(ghost1d),yy(ghost1d),lvlSet.gNx,lvlSet.gNy,'k');
    plot(xx(ghost1d),yy(ghost1d),'bx');
    axis equal
    hold off;
    drawnow; pause
end

%   if ( length(xv) > 55 )
%        keyboard
%   end
    
function aLvl = lvlPoly(obsDat,Grid)

Nx = Grid.Nx; Ny = Grid.Ny;
xc = Grid.xc; yc = Grid.yc;
[yy xx] = meshgrid(yc, xc);
[jj ii] = meshgrid(1:Ny,1:Nx);
ds = max(Grid.dx,Grid.dy);


xv = obsDat.xv;
yv = obsDat.yv;

%Close polygon if necessary
if ((xv(1) ~= xv(end)) || (yv(1) ~= yv(end))) %Maybe change to eps comparison?
    xv = [xv' ; xv(1)];
    yv = [yv' ; yv(1)];
end
Nv = length(xv);

%Construct quick bounding circle
% xcm = sum(xv(1:Nv-1))/(Nv-1);
% ycm = sum(yv(1:Nv-1))/(Nv-1);
% Rcm = sqrt( (xv-xcm).^2 + (yv-ycm).^2 );
% Rcm = max(Rcm(:)); %Find maximum distance from center to a vertex
% Rcm = Rcm + (Grid.ng+1)*ds; %Add some breathing room
% rr = sqrt( (xx-xcm).^2 + (yy-ycm).^2 ); %Distance from poly-center
% InBd = rr<=Rcm;

%Instead of bounding circle, using smoothing operator on In
In = inpolygon(xx,yy,xv,yv);
InS = TwoSmooth(1*In,4,4); %Smoothing window
InBd = (InS > 0);

%Create necessary arrays
sd = zeros(Nx,Ny);
nx = sd; ny = sd;
xb = sd; yb = sd;
Vx = sd; Vy = sd;



%Only calculate information for stuff in "In"
%For each point in In, calculate
%Signed distance, closest point on polygon (may be yourself), normal toward
%closest point

Inbd1d = find(InBd);
x1d = xx(Inbd1d); y1d = yy(Inbd1d);

[sd1d xb1d yb1d Nx1d Ny1d] = sd_poly(x1d,y1d,xv,yv,Grid);

% if (sum(InBd(:)) == 0)
%     keyboard
% end

sdM = 2*max(sd1d);
sd(InBd) = sd1d;
%sum(~InBd(:))
sd(~InBd) = sdM;


nx(InBd) = Nx1d;
ny(InBd) = Ny1d;

%For each point, (xb,yb) is the closest point on the boundary
aLvl.xb = xb1d;
aLvl.yb = yb1d;

%Grab velocities of each vertex
Vxv = obsDat.vx;
Vyv = obsDat.vy;

%Find Vx/Vy of each boundary point using velocities of vertices
Vxi = BoundaryInterp(xb1d,yb1d,xv,yv,Vxv); 
Vyi = BoundaryInterp(xb1d,yb1d,xv,yv,Vyv);

Vx(InBd) = Vxi;
Vy(InBd) = Vyi;

%Store values and get outta here
aLvl.sd = sd;
aLvl.nx = nx;
aLvl.ny = ny;
aLvl.vx = Vx;
aLvl.vy = Vy;




function aLvl = lvlPolyBad(obsDat,Grid)

nDs = 1;
ds = max(Grid.dx,Grid.dy);


Nx = Grid.Nx; Ny = Grid.Ny;
xc = Grid.xc; yc = Grid.yc;
[yy xx] = meshgrid(yc, xc);
[jj ii] = meshgrid(1:Ny,1:Nx);

xv = obsDat.xv;
yv = obsDat.yv;
%Close polygon if necessary
if ((xv(1) ~= xv(end)) || (yv(1) ~= yv(end))) %Maybe change to eps comparison?
    xv = [xv ; xv(1)];
    yv = [yv ; yv(1)];
end

xvr = xv(1:end-1); yvr = yv(1:end-1); %Unclosed vertices
sd = inf(Nx,Ny);
nx = zeros(Nx,Ny);
ny = zeros(Nx,Ny);

[In On] = inpolygon(xx,yy,xv,yv); %Calculate who's in/on the polygon

InS = TwoSmooth(1*In,3,3);
Inish = (InS > 0);

In1d = find(Inish);
Nin = length(In1d);

for n=1:Nin
    %Find the K closest vertices
    
    nn = In1d(n);
    x0 = xx(nn); y0 = yy(nn);  
    i0 = ii(nn); j0 = jj(nn);
    dv = sqrt( (xvr-x0).^2 + (yvr-y0).^2 ); %Distances between this point and all vertices
    [dvmins,I] = sort(dv); %Sort distances into ascending
    %Do signed distance
    if ( In(i0,j0) )
        sd(i0,j0) = -1*dvmins(1);
    else
        sd(i0,j0) = dvmins(1);
    end
    
    %How many points are within nDs*ds of you
    Close = ( dvmins < nDs*ds);
    K = max( sum(Close), 1); %At least choose one
    K = min(K,5); %Let's not get crazy
    Vecs = zeros(2,K);
    Wk = zeros(1,K);
    P = [0 0];
    %Now calculate normal by averaging directions to K closest vertices
    if (On(i0,j0)) %Don't use yourself if you are yourself
        k1 = 2; K = 3;
    else
        k1 = 1;
    end
    for k=k1:K
        kk = I(k);
        Vecs(1,k) = xv(kk) - x0;
        Vecs(2,k) = yv(kk) - y0;
        Nvec = sqrt( Vecs(1,k)^2 + Vecs(2,k)^2 );
        
        Vecs(:,k) = Vecs(:,k)/Nvec;
        W(k) = abs(nDs*ds-dvmins(k)); %Weight based on Rmax-d
        P(1) = P(1) + Vecs(1,k)*W(k);
        P(2) = P(2) + Vecs(2,k)*W(k);
    end
    %P = sum(Vecs,2);
    nvec = sqrt( P(1).^2 + P(2).^2 );
    nx(i0,j0) = P(1)/nvec;
    ny(i0,j0) = P(2)/nvec;
end
aLvl.sd = sd;
aLvl.nx = nx;
aLvl.ny = ny;

function aLvl = lvlPolyOld(obsDat,Grid)

Nx = Grid.Nx; Ny = Grid.Ny;
xc = Grid.xc; yc = Grid.yc;
[yy xx] = meshgrid(yc, xc);
[jj ii] = meshgrid(1:Ny,1:Nx);

xv = obsDat.xv;
yv = obsDat.yv;
%Close polygon if necessary
if ((xv(1) ~= xv(end)) || (yv(1) ~= yv(end))) %Maybe change to eps comparison?
    xv = [xv ; xv(1)];
    yv = [yv ; yv(1)];
end

sd = zeros(Nx,Ny);
nx = zeros(Nx,Ny);
ny = zeros(Nx,Ny);

[In On] = inpolygon(xx,yy,xv,yv); %Calculate who's in/on the polygon

%Note, this is currently stupid, I should fix it
for i=1:Nx
    for j=1:Ny
        x = xx(i,j); y = yy(i,j);
        [dij xp yp] = sd_poly(x,y,xv,yv);
        sd(i,j) = dij;
        if (abs(dij) < 2*eps)
            %On the edge of polygon
            %First, find which edge you're on
            dvec2 = (xv-x).^2 + (yv-y).^2; %Distance between you and all vertices
            %Find the two smallest, you're on that line
            [dmins, I] = sort(dvec2);
            xv0 = xv(I(1)); yv0 = yv(I(1));
            xv1 = xv(I(2)); yv1 = yv(I(2));
            %ie, (x,y) is on the line connecting (xv0,yv0) -> (xv1,yv1)
            %Now construct normal
            %...
            disp('Exact equality of cell centers/polygon not supported');
            pause
        else
            %Inside or outside
            %Calculate vector pointing from us to closest point (xp,yp)
            Px = (xp-x); Py = (yp-y);
            nvec = sqrt(Px^2 + Py^2);
            nx(i,j) = Px/nvec; ny(i,j) = Py/nvec;
        end
    end
end

aLvl.sd = sd;
aLvl.nx = nx;
aLvl.ny = ny;

function aLvl = lvlCircle(Grid,x0,y0,rad)

Nx = Grid.Nx; Ny = Grid.Ny;
xc = Grid.xc; yc = Grid.yc;
[yy xx] = meshgrid(yc, xc);
[jj ii] = meshgrid(1:Ny,1:Nx);

rr = sqrt( (xx - x0).^2  + (yy - y0).^2 ) ;
In =  rr <= rad;

sd = (rr-rad);
Px = (xx-x0); Py = (yy-y0);
nvec = sqrt(Px.^2 + Py.^2);

aLvl.sd = sd;
aLvl.nx = Px./nvec;
aLvl.ny = Py./nvec;
aLvl.vx = zeros(Nx,Ny); aLvl.vy = aLvl.vx;

%Calculates relevant values to move object
%Displacement from center, velocity of boundary

%if obsDat.mobile = true, then
%obsDat.dx (1/2 elements)
%obsDat.tau (1/2 elements)
%obsDat.func (1/2 elements, 0 = sin, 1 = sin)

%x = x0 + dx*func(2*pi*t/tau_x)
%y = ...
%delx = x-x0
function [delx dely vx vy] = moveObj(obsDat,t)

dx = obsDat.dx(1);
if (length(obsDat.dx) > 1)
    dy = obsDat.dx(2);
else
    dy = dx;
end

taux = obsDat.tau(1);
if (length(obsDat.tau) > 1)
    tauy = obsDat.tau(2);
else
    tauy = taux;
end

func = obsDat.func;
fx = obsDat.func(1);
if (length(func) > 1)
    fy = func(2);
else
    fy = fx;
end

[delx vx] = moveDir(dx,taux,fx,t);
[dely vy] = moveDir(dy,tauy,fy,t);


function [del v] = moveDir(ds,taus,fs,t)

if (fs == 0)
    del = ds*sin(2*pi*t/taus);
    v = (2*pi/taus)*ds*cos(2*pi*t/taus);
elseif (fs == 1)
    del = ds*cos(2*pi*t/taus);
    v = -1*(2*pi/taus)*ds*sin(2*pi*t/taus);
else
    disp('Unknown functional form');
    pause;
end
    


    