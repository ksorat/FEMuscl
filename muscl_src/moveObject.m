
%Moves the objects defined in Model.lvlDef
function Model = moveObject(Model,Grid,Gas)

%Useful grabs
lvlDef = Model.Init.lvlDef;
Nx = Grid.Nx; Ny = Grid.Ny;
numObs =lvlDef.numObs;

lvlSet.ds = max(Grid.dx,Grid.dy);
lvlSet.ds_min = min(Grid.dx,Grid.dy);
dt = Grid.dt;

for n=1:numObs
    obsDat = lvlDef.obsDat{n};
    M = 1; %Total mass of polygon
    I = 0.1; %Moment of inertia
    
    if isfield(obsDat,'Fx') %Negative sign from fluid->solid
        %size(obsDat.vx),size(obsDat.Fx)
        Nv = length(obsDat.xv);
        xvr = obsDat.xv(1:Nv-1); yvr = obsDat.yv(1:Nv-1);
        vxr = obsDat.vx(1:Nv-1); vyr = obsDat.vy(1:Nv-1);
        Fxr = obsDat.Fx(1:Nv-1);
        Fyr = obsDat.Fy(1:Nv-1);
        
        xcm = sum(xvr)/(Nv-1);
        ycm = sum(yvr)/(Nv-1);
        
        Fx_cm = sum(Fxr); Fy_cm = sum(Fyr);
        
        Tz = 0;
        for v=1:Nv-1
            xp = xvr(n)-xcm;
            yp = yvr(v)-ycm;
            T = cross( [xp yp 0], [Fxr(v) Fyr(v) 0] );
            Tz = Tz + T(3);
        end %Now have total torque
        
        ax_cm = -Fx_cm/M; ay_cm = -Fy_cm/M;
        alpha = -Tz/I;
        
        obsDat.Vcm_x = obsDat.Vcm_x + dt*ax_cm;
        obsDat.Vcm_y = obsDat.Vcm_y + dt*ay_cm;
        
        obsDat.omega = obsDat.omega + dt*alpha;
        for v=1:Nv-1
            xp = xvr(v) - xcm;
            yp = yvr(v) - ycm;
            
            Vrot = -1*cross( [xp yp 0] , [0 0 obsDat.omega] );
            
            vxr(v) = obsDat.Vcm_x + Vrot(1);
            vyr(v) = obsDat.Vcm_y + Vrot(2);
            
            theta = dt*obsDat.omega;
            xr = cos(theta)*xp - sin(theta)*yp;
            yr = sin(theta)*xp + cos(theta)*yp;
            
            xvr(v) = xcm + dt*obsDat.Vcm_x + xr;
            yvr(v) = ycm + dt*obsDat.Vcm_y + yr;
        end
        obsDat.vx(1:Nv-1) = vxr; obsDat.vx(end) = obsDat.vx(1);
        obsDat.vy(1:Nv-1) = vyr; obsDat.vy(end) = obsDat.vy(1);
        
        obsDat.xv(1:Nv-1) = xvr; obsDat.xv(end) = obsDat.xv(1);
        obsDat.yv(1:Nv-1) = yvr; obsDat.yv(end) = obsDat.yv(1);
            
    else
        obsDat.xv = obsDat.xv + dt*obsDat.vx;
        obsDat.yv = obsDat.yv + dt*obsDat.vy;
        
    end
    
    lvlDef.obsDat{n} = obsDat;
end
Model.Init.lvlDef = lvlDef;