
function Grid = InitSolid(Model,Grid)
%Take data from Model.Init.obsDef to generate obstructions via legoland
%approach

%Constructs details of obstruction region given data-type obsDef
%obsDef contains
%obsDef.numObs = number of obstructions
%obsDef.obsType(1:numObs) = type of each obstruction
%obsDef.obsParam(1:numObs,1:4) = details for each, requires 4 numbers
%block-type, the 4 numbers are [xmin xmax ymin ymax] 
%wedge-type, the 4 numbers are [x0 y0 theta r] 
%circle-type, the 4 numbers are [xc yc rad xx?] 


%Generate mesh
[xx,yy] = meshgrid(Grid.xc,Grid.yc);
xx = xx'; yy = yy';

obsDef = Model.Init.obsDef;
numObs = obsDef.numObs; %The number of obstructions
obsType = obsDef.obsType;
obsParam = obsDef.obsParam;

for i=1:numObs
    switch lower(obsType{i})
        case{'circle'}
            %Do circle type
            xc = obsParam(i,1);
            yc = obsParam(i,2);
            radc = obsParam(i,3);
            rad = sqrt( (xx-xc).^2 + (yy-yc).^2);
            Mark{i} = (rad <= radc);
            
            
        case{'block'}
            %Do block type
            xb = [obsParam(i,1) obsParam(i,2)];
            yb = [obsParam(i,3) obsParam(i,4)];
            
            delxb = xb(2)-xb(1);
            delyb = yb(2)-yb(1);
%             if (delxb < 1.5*NumG) || (delyb < 1.5*NumG)
%                 disp('Under-resolved obstruction, halted');
%                 keyboard
%                 %pause
%             end
            Mark{i} = (xx > xb(1)) & (xx < xb(2)) & (yy > yb(1)) & (yy < yb(2));
        case{'trix'}
            %Create x-aligned isosceles triangle
            %4 points are x1/y1 x2/y2, west vertix and north-east points
            x1 = obsParam(i,1);
            y1 = obsParam(i,2);
            x2 = obsParam(i,3);
            y2 = obsParam(i,4);
            x3 = x2;
            y3 = y1-(y2-y1);
            
            %Now use geometry
            atop = (y2 - y3)*(xx-x3) + (x3 - x2)*(yy - y3);
            bot = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);
            a = atop/bot;
            btop = (y3-y1)*(xx-x3) + (x1-x3)*(yy-y3);
            b = btop/bot;
            c = 1 - a - b;
            
            Mark{i} = (a >= 0) & (a <= 1) & (b >= 0) & (b <= 1) & (c >=0) & (c <= 1);
            
            
        case{'wedge'}
            %Do wedge type obstruction
            %Define points of triangle (x1,y1),(x2,y2) ...
            x1 = obsParam(i,1);
            y1 = obsParam(i,2);
            theta = obsParam(i,3); %assuming degrees!
            rad = obsParam(i,4);
            
            %Define other points of triangle
            x2 = x1 + rad*cosd(theta);
            y2 = y1 + rad*sind(theta);
            x3 = x2;
            y3 = y1 - rad*sind(theta);

            
            %Now use geometry
            atop = (y2 - y3)*(xx-x3) + (x3 - x2)*(yy - y3);
            bot = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);
            a = atop/bot;
            btop = (y3-y1)*(xx-x3) + (x1-x3)*(yy-y3);
            b = btop/bot;
            c = 1 - a - b;
            
            Mark{i} = (a >= 0) & (a <= 1) & (b >= 0) & (b <= 1) & (c >=0) & (c <= 1);
    end
end

Markfin = Mark{1};
for i=1:numObs
    Markfin = Mark{i} | Markfin;
end

%obReg.reg = Markfin;
%obReg.reg1d = find(Markfin);


solid.In = Markfin;
%solid.In1d = find(Markfin); %Are you necessary?


smoothIn = SmoothCross(1*solid.In);
bdry = (smoothIn < 1) & (smoothIn > 0);
solid.BdrySol = bdry & solid.In; %Boundary of solid
solid.BdryFlu = bdry & (~solid.In); %Boundary of fluid
solid.bf1d = find(solid.BdryFlu);
%Now calculate necessary quantities from solid geometry
%Fx_sol - x interfaces that connect fluid/solid
%Fy_sol - y ...
%BdryFlu = Bdry_xm U Bdry_xp U Bdry_ym U Bdry_yp
%For each fluid cell on the boundary identify its neighbor configuration,
%ie is the xp (plus-x) neighbor a solid, etc
%Note, fluid cell can be on corner (both an x/y neighbor) but can't have
%two x neighbors


[solid.Fx_sol solid.Fy_sol solid.Bdry_xm solid.Bdry_xp solid.Bdry_ym solid.Bdry_yp] = FindGeom(Grid,solid);

Grid.solid = solid;



function [Fx_sol Fy_sol Bdry_xm Bdry_xp Bdry_ym Bdry_yp Fx_pc Fy_pc] = FindGeom(Grid,sol)
Nx = Grid.Nx;
Ny = Grid.Ny;
Nxi = Grid.Nx+1;
Nyi = Grid.Ny+1;

Fx_sol = false(Nxi,Ny); %Number of x-interfaces
Fy_sol = false(Nx,Nyi); %Number of y-interfaces

Bdry_xp = false(Nx,Ny);
Bdry_xm = false(Nx,Ny);
Bdry_yp = false(Nx,Ny);
Bdry_ym = false(Nx,Ny);


%Loop through boundary

bdry1d  = sol.bf1d;
%Create 2d i/j matrices
[ii,jj] = meshgrid(1:Nx,1:Ny);
ii = ii'; jj = jj';

In = sol.In;

for i=1:Nx
    for j=1:Ny
        if In(i,j)
            Fx_sol(i,j) = true;
            Fx_sol(i+1,j) = true;
            Fy_sol(i,j) = true;
            Fy_sol(i,j+1) = true;
        end
    end
end

% %Loop through fluid boundary cells and construct nec. info
% for n=1:length(bdry1d)
%     elt = bdry1d(n);
%     i = ii(elt);
%     j = jj(elt);
%     N = 0;
%     if In(i+1,j)
%         Bdry_xp(i,j) = true;
%         Fx_sol(i+1,j) = true;
%         N = N+1;
%     end
%     if In(i-1,j)
%         Bdry_xm(i,j) = true;
%         Fx_sol(i,j) = true;
%         N = N+1;
%     end
%     if In(i,j+1)
%         Bdry_yp(i,j) = true;
%         Fy_sol(i,j+1) = true;
%         N = N+1;
%     end
%     if In(i,j-1)
%         Bdry_ym(i,j) = true;
%         Fy_sol(i,j-1) = true;
%         N = N+1;
%     end
% 
% end

%Make sure no cell is on both a plus/minus boundary
Chk = Bdry_xm+Bdry_xp+Bdry_ym+Bdry_yp;
Chk = Chk > 2;
if (sum(Chk(:)) > 0)
    disp('Poorly resolved boundary');
    keyboard
end

% %Check coverage
% Chk = Bdry_xm | Bdry_xp | Bdry_ym | Bdry_yp;
% Chk = Chk - sol.BdryFlu;


function Zs = SmoothCross(Z)
[Nx Ny] = size(Z);
Zs = Z;
Zs(2:Nx-1,2:Ny-1) = 0.25*( Z(1:Nx-2,2:Ny-1) + Z(3:Nx,2:Ny-1) + Z(2:Nx-1,1:Ny-2) + Z(2:Nx-1,3:Ny) );

%--------------------
%Below this line is old method
% 
% %Use data from Model to initialize the ib nodes
% %Format of Model.ib
% %numObs = number of obstructions
% %Model.ib.obs{i} = data from the ith obstruction
% %Model.ib.obs{i}.type = 'circle' / polygon
% %   If circle, subfields cent, rad, N (number of ib nodes)
% 
% 
% %Creates as a subfield of Grid, Grid.ib
% %Grid.ib.numObs
% %Grid.ib.obs{i}
% %Grid.ib.obs{i}.x / .y / .Vx / .Vy (positions and velocities of each node
% %Grid.ib.obs{i}.N = nodes in this obstruction
% %Grid.ib.numNodes = Total number of nodes
% %Grid.ib.Bdry = Boolean array, true at cell ij if ij contains node or edge
% %connecting nodes
% %Grid.ib.Inside = Boolean aray, true if cell is enclosed in boundary
% 
% numObs = Model.ib.numObs; %Number of obstructions
% 
% Grid.ib.numObs = numObs;
% Grid.ib.numNodes = 0;
% 
% Nx = Grid.Nx; Ny = Grid.Ny;
% [xx yy] = meshgrid(Grid.xc,Grid.yc); xx = xx'; yy = yy';
% %xxi = xx(IntX,IntY); yyi = yy(IntX,IntY);
% 
% %Create boolean arrays
% %Bdry - Intersected by curve Gamma, ie points that do what solid tells them
% %Inside - Interior of curve Gamma, ie not really fluid cells
% %Nbrhd - Region within one cell of boundary, fluid cells that can inform
% %behavior of Bdry points along with solid.  
% 
% Grid.ib.Bdry = false(Nx,Ny); 
% Grid.ib.Inside = false(Nx,Ny);
% Grid.ib.Nbrhd = false(Nx,Ny);
% 
% for i=1:numObs
%     obsDat = Model.ib.obs{i};
%    %Create i-th obstruction
%    switch lower(obsDat.type)
%        case{'circle'}
%            N = obsDat.N; x0 = obsDat.cent(1); y0 = obsDat.cent(2); rad = obsDat.rad;
%            tau = linspace(0,2*pi,N+1);
%            tau = linspace(0,2*pi,N);
%            tau = tau(1:N);
%            Grid.ib.obs{i}.x = rad*cos(tau) + x0;
%            Grid.ib.obs{i}.y = rad*sin(tau) + y0;
%            %Assuming stationary boundary for now
%            Grid.ib.obs{i}.vx = zeros(1,N);
%            Grid.ib.obs{i}.vy = zeros(1,N);
%            Grid.ib.obs{i}.N = N;
%            Grid.ib.numNodes = Grid.ib.numNodes + N;
%        case{'poly'}
%            gridObs = PolyObs(Grid,obsDat);
%            Grid.ib.obs{i} = gridObs;
%            Grid.ib.numNodes = Grid.ib.numNodes + gridObs.N; 
%            
% 
%    end
%    
%    [LocBdry LocInside LocNbr] = obsGeom(Grid,Grid.ib.obs{i});
%    
%    Grid.ib.Bdry = Grid.ib.Bdry | LocBdry;
%    Grid.ib.Inside = Grid.ib.Inside | LocInside;
%    Grid.ib.Nbrhd = Grid.ib.Nbrhd | LocNbr;
%    
%            
% end
% 
% 
% function [Bdry Inside Nbr] = obsGeom(Grid,Obs)
% 
% Nsm = 2; %Radius of discrete smoothing window
% Nx = Grid.Nx; Ny = Grid.Ny;
% [xx yy] = meshgrid(Grid.xc,Grid.yc); xx = xx'; yy = yy';
% 
% Inside = inpolygon(xx,yy,Obs.x,Obs.y);
% Bdry = obsBoundary(Grid,Obs); %Note, not a very good calculation here
% 
% %Now figure out neighborhood of interest
% %Also not a very good calculation, but take union of Inside and bdry, then
% %smooth and find intermediate values
% 
% %wNbr = TwoSmooth(1*(Inside | Bdry),Nsm,Nsm); %Smooth boolean
% wNbr = TwoSmooth(1*(Bdry),Nsm,Nsm); %Smooth boolean
% Nbr = (wNbr > 0) & (wNbr < 1);
% 
% function Bdry = obsBoundary(Grid,Obs)
% %Note, this is likely more accurate than other calculation but slower
% Nx = Grid.Nx; Ny = Grid.Ny;
% dx = Grid.dx; dy = Grid.dy;
% 
% Bdry = false(Nx,Ny);
% x = Obs.x;
% y = Obs.y;
% [is ie js je] = Localize(Grid,x,y);
% 
% for i=is:ie
%     for j=js:je
%         xc = Grid.xc(i); yc = Grid.yc(j);
%         cellx = [xc - 0.5*dx, xc-0.5*dx, xc+0.5*dx, xc+0.5*dx];
%         celly = [yc - 0.5*dx, yc+0.5*dx, yc+0.5*dx, yc-0.5*dx];
%         [xcut ycut] = polyxpoly(cellx,celly,x,y);
%         if length(xcut) > 0
%             Bdry(i,j) = true;
%         end
%     end
% end
%         
% 
% 
% 
% %Given a single boundary curve, find min/max i/j on the grid which contain
% %it
% function [is ie js je] = Localize(Grid,x,y)
% Nx = Grid.Nx; Ny = Grid.Ny;
% xi = Grid.xi; yi = Grid.yi;
% 
% xmin = min(x);
% xmax = max(x);
% ymin = min(y);
% ymax = max(y);
% 
% is = find( xmin >= xi(1:Nx),1,'last' ) - 1;
% ie = find( xmax <= xi(2:Nx+1),1,'first' ) + 1;
% js = find( ymin >= yi(1:Ny),1,'last' ) - 1;
% je = find( ymax <= yi(2:Ny+1),1,'first' ) + 1;
% 
% 
% function bdry = obsBoundaryOld(Grid,Obs)
% Nx = Grid.Nx; Ny = Grid.Ny;
% x = Obs.x;
% y = Obs.y;
% 
% bdry = false(Nx,Ny);
% 
% %Loop through nodes
% for n=1:Obs.N
%     xn = Obs.x(n);
%     yn = Obs.y(n);
%     [in jn] = FindCell(Grid,xn,yn);
%     bdry(in,jn) = true;
% end
% 
% %Loop through edges
% %Write this
% 
% 
% function [i j] = FindCell(Grid,xn,yn)
% 
% xi = Grid.xi; yi = Grid.yi;
% Nx = Grid.Nx; Ny = Grid.Ny;
% 
% xlog = (xn >= xi(1:Nx-1) ) & (xn < xi(2:Nx) );
% i = find(xlog,1,'first');
% 
% ylog = (yn >= yi(1:Ny-1) ) & (yn < yi(2:Ny) );
% j = find(ylog,1,'first');
% 
% function gridObs = PolyObs(Grid,obsDat)
% Nper = obsDat.Np; %Number of points per segment
% Nseg = obsDat.Ns; %Number of segments
% xs = obsDat.xs;
% ys = obsDat.ys;
% %Close the circle
% xs = [xs obsDat.xs(1)];
% ys = [ys obsDat.ys(1)];
% 
% gridObs.x = [];
% gridObs.y = [];
% for n=1:Nseg
%     x0 = xs(n);
%     y0 = ys(n);
%     x1 = xs(n+1);
%     y1 = ys(n+1);
%     delx = x1-x0; dely = y1-y0;
%     if (abs(delx) < eps)
%         %Handle vertical line
%         xhf = 0.5*(x1+x0);
%         x = xhf*ones(1,Nper);
%         y = linspace(y0,y1,Nper+1); y = y(1:Nper);
%     else
%         m = dely/delx;
%         x = linspace(x0,x1,Nper+1); x = x(1:Nper);
%         y = m*(x-x0)+y0;
%     end
%     gridObs.x = [gridObs.x x];
%     gridObs.y = [gridObs.y y];
% end
% 
% 
% gridObs.N = length(gridObs.x);
% gridObs.vx = zeros(size(gridObs.x));
% gridObs.vy = gridObs.vx;
% 


