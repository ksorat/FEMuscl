
%Calculates interface geometries from lvlSet data
%Returns lvlSet.(Dir)
%Marks fluid cells that are adjacent to ghosts

%West Fluid->Ghost, Flu(i,j) & Gh(i+1,j)
%East Ghost->Fluid, Flu(i,j) & Gh(i-1,j)
%North Ghost->Fluid Flu(i,j) & Gh(i,j-1)

%Also calculate "fake" interfaces
%A fake interface is anything connected to an object cell
%These interfaces will zero out flux
%Fx_obj(i,j) refers to the interface connecting cell (i,j) -> (i-1,j)

function lvlSet = calcGeom(Grid,lvlSet)
Gh = lvlSet.ghost;
Flu = lvlSet.fluid;

gh1d = lvlSet.ghost1d;

Nx = Grid.Nx; Ny = Grid.Ny;
Nxi = Nx+1; Nyi = Ny+1;

Fx_obj = false(Nxi,Ny);
Fy_obj = false(Nx,Nyi);
[jj ii] = meshgrid(1:Ny,1:Nx);

%Could easily make these sparse
West = false(Nx,Ny); East = West; North = West; South = West;

Ngh = length(gh1d);

%Loop through ghost cells
for n=1:Ngh
    i = lvlSet.gi(n);
    j = lvlSet.gj(n);
    %Check neighbors
    if Flu(i+1,j)
        East(i+1,j) = true;
    end
    if Flu(i-1,j)
        West(i-1,j) = true;
    end
    if Flu(i,j+1)
        North(i,j+1) = true;
    end
    if Flu(i,j-1)
        South(i,j-1) = true;
    end
        
end

%Loop through object cells
Obj = lvlSet.obj;
obj1d = find(Obj);
Nob = length(obj1d);
obji = ii(obj1d); objj = jj(obj1d);
for n=1:Nob
    i = obji(n); j = objj(n);
    Fx_obj(i,j) = true;
    Fx_obj(i+1,j) = true;
    Fy_obj(i,j) = true;
    Fy_obj(i,j+1) = true;
end

lvlSet.East  = East;
lvlSet.West  = West;
lvlSet.North = North;
lvlSet.South = South;

lvlSet.Fx_obj = Fx_obj;
lvlSet.Fy_obj = Fy_obj;



