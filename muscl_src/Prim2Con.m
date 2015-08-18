function [D Mx My E] = Prim2Con(D,Vx,Vy,P,Model)
Gam = Model.Init.Gam;

global SMALL_NUM;

%K: Removed check here
%D = max(D,SMALL_NUM);
%P = max(P,SMALL_NUM);

Mx = D.*Vx;
My = D.*Vy;
K = 0.5*D.*( Vx.^2 + Vy.^2);
e = P/(Gam-1);
E = e + K;
