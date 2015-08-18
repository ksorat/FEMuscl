
function Cs = Prim2Cs(D,Vx,Vy,P,Model)

Gam = Model.Init.Gam;

Cs = realsqrt( Gam*P./D );

