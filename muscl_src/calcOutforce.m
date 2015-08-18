%Takes in x/y interface fluxes
%Returns outForce which contains

function [outForce Model] = calcOutforce(Grid,Model,Fl_x,Fl_y)

Nx = Grid.Nx; Ny = Grid.Ny;
dx = Grid.dx; dy = Grid.dy;
hdx = 0.5*dx; hdy = 0.5*dy;
xc = Grid.xc; yc = Grid.yc;
%Create embiggen'ed arrays
[yy xx] = meshgrid(yc,xc);
[jj ii] = meshgrid(1:Ny,1:Nx);

Mx_x = dx*squeeze( Fl_x(2,:,:) ); %Mx flux in x direction
Mx_y = dy*squeeze( Fl_y(2,:,:) ); %Mx flux in y direction
My_x = dx*squeeze( Fl_x(3,:,:) ); %My flux in x direction
My_y = dy*squeeze( Fl_y(3,:,:) ); %My flux in y direction

lvl = Grid.lvlSet;
lvlDef = Model.Init.lvlDef;

%Go back and think about scaling here

Any = lvl.West | lvl.East | lvl.South | lvl.North;
Ntot = sum( lvl.West(:) + lvl.East(:) + lvl.South(:) + lvl.North(:) );
any1d = find(Any);
N = length(any1d);

xp = zeros(1,Ntot);
yp = xp; Fx = xp; Fy = xp;

k = 1;

%Note, Ntot != N, Ntot is total number of points (may have 2 interfaces in
%one cell)
%Construct a series of points and the Fx/Fy forces at each
for n=1:N %Loop over cells that are on boundary
    nc = any1d(n); %Cell number
    i = ii(nc); j = jj(nc);
    if lvl.West(i,j)
        xp(k) = xc(i) + hdx;
        yp(k) = yc(j);
        Fx(k) = -Mx_x(i+1,j);
        Fy(k) = -My_x(i+1,j);
        k = k+1;
    end
    
    if lvl.East(i,j)
        xp(k) = xc(i) - hdx;
        yp(k) = yc(j);
        Fx(k) = Mx_x(i,j);
        Fy(k) = My_x(i,j);
        k = k+1;
    end
    
    if lvl.South(i,j)
        xp(k) = xc(i);
        yp(k) = yc(j)+hdx;
        Fx(k) = -Mx_y(i,j+1);
        Fy(k) = -My_y(i,j+1);
        k = k+1;
    end
    
    if lvl.North(i,j)
        xp(k) = xc(i);
        yp(k) = yc(j)-hdx;
        Fx(k) = Mx_y(i,j);
        Fy(k) = My_y(i,j);
        k = k+1;
    end    
    
end


% plot(xp,yp,'ro'); hold on;
% for n=1:lvlDef.numObs
%     obsDat = lvlDef.obsDat{n};
%     plot(obsDat.xv,obsDat.yv,'bx');
% end
% hold off;

%Now we have to map from interface points that house forces to the vertices
%of the polygon
%This is currently crude and should likely be improved later

%Start by creating one big list of vertices and their corresponding
k = 1;
xv = []; yv = []; obj = [];
for n=1:lvlDef.numObs
    Nv = length(lvlDef.obsDat{n}.xv);
    xv = [xv lvlDef.obsDat{n}.xv];
    yv = [yv lvlDef.obsDat{n}.yv];
    obj(k: k+Nv-1) = n;
    k = k+Nv;
end

Nvtot = length(xv);
fxv = zeros(1,Nvtot);
fyv = zeros(1,Nvtot);

%Now loop through interface points, for each interface point find which
%vertex you're closest to and apply your force to that
% 
for n=1:Ntot
   x = xp(n); y = yp(n);
   dv = sqrt( (xv-x).^2 + (yv-y).^2 );
   [dvmin,I] = min(dv);
   %You're closest to the Ith vertex, give it your force
   fxv(I) = fxv(I) + Fx(n);
   fyv(I) = fyv(I) + Fy(n);
   %For giggles, draw vector from you to your vertex
%    hold on; plot([x xv(I)],[y yv(I)],'g'); hold off
end

%Correct for closed loop
fxv(1) = fxv(1) + fxv(end);
fxv(end) = fxv(1);
fyv(1) = fyv(1) + fyv(end);
fyv(end) = fyv(1);

%Finally put the forces back in

for n=1:lvlDef.numObs
    Ind = find(obj == n);
    lvlDef.obsDat{n}.Fx = -1*fxv(Ind);
    lvlDef.obsDat{n}.Fy = -1*fyv(Ind);
end


Fxtot = sum(fxv); Fytot = sum(fyv);
%fprintf('\tTotal Fx = %d / Total Fy = %d\n', Fxtot,Fytot);

Model.Init.lvlDef = lvlDef;
outForce.Fxtot = Fxtot;
outForce.Fytot = Fytot;

