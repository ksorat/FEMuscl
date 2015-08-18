%Calculates timestep on grid
function dt = CalcDT(Model,Grid,Gas)

%Find min grid size
dxmin = min(Grid.dx,Grid.dy);

Vabs = sqrt(Gas.Vx.^2 + Gas.Vy.^2);

Cs = Prim2Cs(Gas.D,Gas.Vx,Gas.Vy,Gas.P,Model);

if (Model.lvlSet.anymobile)
    Vxm = max(abs(Grid.lvlSet.gVx(:)));
    Vym = max(abs(Grid.lvlSet.gVy(:)));
    Vabs = Vabs + max(Vxm,Vym);
end

dt2d = dxmin./(Cs + Vabs);

dt = Grid.C0*min(dt2d(:));

if (dt < 1.0e-8)
    disp('Unreasonably small timestep');
    keyboard
end


