%Sets internal ghost zones for the level set method

function Gas = LvlGhost(Model,Grid,Gas)

if (~Model.lvlSet.present)
    %How did you even get here?
    disp('You shouldn''t be in LvlGhost');
end

lvlSet = Grid.lvlSet;

%Loop through interior ghost zones
for n=1:lvlSet.ng
    %i/j indices of this ghost zone
    ig = lvlSet.gi(n);
    jg = lvlSet.gj(n);
    
    %x/y point of this ghost zone (centered)
    xg = Grid.xc(ig);
    yg = Grid.yc(jg);
    
    %Signed distance/normals of this ghost zone from boundary
    sdn = lvlSet.ghost_sd(n);
    nx = lvlSet.gNx(n);
    ny = lvlSet.gNy(n);
    
    %Distance from ghost zone center to image point
    L = abs(sdn) + lvlSet.dip;
    xip = xg + L*nx;
    yip = yg + L*ny;
    
    %Now we have image point, find 4 closest cells to image point
    Img = conImage(xip,yip,Grid);
    
    %Now we have info about the 4 closest cells, prep the bilinear
    %interpolation
    
    %Find the density and pressure at the image point
    %Assuming neuman BC w/ normal deriv=0 at boundary
    Dip = calcInterp(Img,Gas.D,Grid);
    Pip = calcInterp(Img,Gas.P,Grid);
    
    
    %Now calculate velocity to zero out slip at boundary
    Vxip = calcInterp(Img,Gas.Vx,Grid);
    Vyip = calcInterp(Img,Gas.Vy,Grid);
    scl = L/(L - abs(sdn) );
    
    Vxbi = lvlSet.gVx(n);
    Vybi = lvlSet.gVy(n);
    
    Vxgc = Vxip - scl*( Vxip - Vxbi );
    Vygc = Vyip - scl*( Vyip - Vybi );

    if (Img.g2g)
        Dip = Gas.D(ig,jg);
        Pip = Gas.P(ig,jg);
        Vxgc = Gas.Vx(ig,jg);
        Vygc = Gas.Vy(ig,jg);
    end

    Gas.D(ig,jg) = Dip;
    Gas.P(ig,jg) = Pip;
    
    Gas.Vx(ig,jg) = Vxgc;
    Gas.Vy(ig,jg) = Vygc;
    
end

function [Vx Vy] = calcVBi_trash(Model,Grid)

if (Model.lvlset.mobile)
    %Assuming sinusoidal structure
    ds = Model.Init.lvlDef.ds;
    tau = Model.Init.lvlDef.tau;
    t = Grid.t;
    
    scl = (2*pi/tau)*cos(2*pi*t/tau);
    Vx = ds(1)*scl;
    Vy = ds(2)*scl;
else
    Vx = 0.0; Vy = 0.0;
end

%Calculates value of Z @ image point given 4 data points
function Zip = calcInterp(Img,Z,Grid)

%Create helpers
i1 = min(Img.i); i2 = max(Img.i);
j1 = min(Img.j); j2 = max(Img.j);

x1 = Grid.xc(i1); x2 = Grid.xc(i2);
y1 = Grid.yc(j1); y2 = Grid.yc(j2);

z11 = Z(i1,j1);
z12 = Z(i1,j2);
z21 = Z(i2,j1);
z22 = Z(i2,j2);


%Calculate interpolated value at x/y
x = Img.xip; y = Img.yip;
scli = (x2-x1)*(y2-y1);
scl = 1/scli;

Zip = scl* ( z11*(x2-x)*(y2-y) + z21*(x-x1)*(y2-y) + z12*(x2-x)*(y-y1) + z22*(x-x1)*(y-y1) );


%Creates Img data structure, includes
%Img.i = 4 values of i for the 4 closest cells
%Img.j = same, but with j
%Img.x = same, but with x values of cell centers
%Img.y = same, but with y values of cell centers

function Img = conImage(xip,yip,Grid)

%Image point is in cell i/j
i = find(xip>=Grid.xi,1,'last');
j = find(yip>=Grid.yi,1,'last');


x = Grid.xc(i); y = Grid.yc(j);

if (xip >= x)
    %look to i+1
    i1 = i+1;
else
    i1 = i;
    i = i-1;
end

if (yip >= y)
    j1 = j+1;
else
    j1 = j;
    j = j-1;
end

Img.i = [i i i1 i1];
Img.j = [j j1 j j1];

Img.x = Grid.xc(Img.i);
Img.y = Grid.yc(Img.j);

Img.xip = xip; Img.yip = yip;

%Check for ghost 2 ghost contact
g2g = false;
for n=1:4
    if Grid.lvlSet.ghost(Img.i(n),Img.j(n))
        g2g = true;
    end
end

Img.g2g = g2g;
