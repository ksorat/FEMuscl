function [Flux Fluxl] = RiemannFlux(Wl,Wr,Model,solver,Grid)
%Calculates Riemann flux using 'solver'
% |Wl| = |Wr| = [Nvar,Nxi,Nyi]
% Assume L/R based on interfaces and not cells, ie Wl(i,j) is the
% L state of the (i-1/2,j) interface

%Note, we assume 2D and that we are solving along the x-direction
%Y-direction is handled by interchanging x<->y

switch lower(solver)
    case{'hlle'}
        Flux = RiemannFlux_HLLE(Wl,Wr,Model);
        Fluxl = [];
    case{'hllc'}
        Flux = RiemannFlux_HLLC(Wl,Wr,Model);
        Fluxl = [];
    case{'hll'}
        %Calculate both HLLE/HLLC fluxes
        [Fluxl Flux] = RiemannFlux_HLL(Wl,Wr,Model);
        
    otherwise
        disp('Unknown Riemann solver type');
        pause
end

%Calculates both HLL(E+C) fluxes at once to avoid redundant calculation
function [Fle Flc] = RiemannFlux_HLL(Wl,Wr,Model)
Fle = zeros(size(Wl));
Flc = Fle;

%Calculate Roe averages
[droe vxroe vyroe hroe] = RoeAverages(Wl,Wr,Model);

evals = RoeEigenvals(droe,vxroe,vyroe,hroe,Model);

%Calculate min/max wavespeeds of L/R states
cfl = P2Cs(Wl,Model);
cfr = P2Cs(Wr,Model);

%Take min/max of Roe's eigenvalues and L/R states
ar = max( evals(4,:,:) , Wr(2,:,:) + cfr );
al = min( evals(1,:,:) , Wl(2,:,:) - cfl );

bp = max(ar,0.0);
bm = min(al,0.0);

%Calculate L/R fluxes
%Note, fluxes apply to *CONSERVED* not primitive vars
[Fl Fr] = CalcLR_Fluxes(Wl,Wr,bm,bp,Model);

%Compute the HLLE flux at each interface
scl = 0.5*(bp + bm)./(bp - bm);

for n=1:4
    Fle(n,:,:) = 0.5*( Fl(n,:,:) + Fr(n,:,:) ) + ( Fl(n,:,:) - Fr(n,:,:) ).*scl;
end

%Now switch gears and do HLLC calculation
%Calculate speed of contact wave and pressure
tl = Wl(4,:,:) + Wl(1,:,:).*Wl(2,:,:).* ( Wl(2,:,:) - al );
tr = Wr(4,:,:) + Wr(1,:,:).*Wr(2,:,:).* ( Wr(2,:,:) - ar );

dl = Wl(1,:,:).*Wl(2,:,:) - Wl(1,:,:).*al;
dr = -1*(Wr(1,:,:).*Wr(2,:,:) - Wr(1,:,:).*ar);

am = (tl - tr) ./ (dl + dr); %Contact wave speed
cp = (dl.*tr + dr.*tl) ./ (dl + dr); %Contact pressure at moving surface
cp = max(cp,0.0); %Enforce positivity on pressure

%Calculate scaling terms for L/R fluxes
%Distinguish right-moving contact wave from left
Ind = (am >= 0.0);
SclR = zeros(size(am)); SclL = SclR; SclM = SclR;

%Right-moving contact wave
SclL(Ind) = am(Ind)./(am(Ind)-bm(Ind));
SclR(Ind) = 0.0;
SclM(Ind) = -bm(Ind)./(am(Ind)-bm(Ind));

%Left-moving contact wave
SclL(~Ind) = 0.0;
SclR(~Ind) = -am(~Ind)./(bp(~Ind)-am(~Ind));
SclM(~Ind) = bp(~Ind)./(bp(~Ind)-am(~Ind));

%Compute HLLC flux
for n=1:4
    Flc(n,:,:) = SclL.* Fl(n,:,:) + SclR.*Fr(n,:,:);
end
%Add correction for contact
Flc(2,:,:) = Flc(2,:,:) + SclM.*cp; %Correct Mx w/ contact pressure
Flc(4,:,:) = Flc(4,:,:) + SclM.*cp.*am; %Correct energy


function Flux = RiemannFlux_HLLE(Wl,Wr,Model)

Flux = zeros(size(Wl));

%Calculate Roe averages
[droe vxroe vyroe hroe] = RoeAverages(Wl,Wr,Model);

evals = RoeEigenvals(droe,vxroe,vyroe,hroe,Model);

%Calculate min/max wavespeeds of L/R states
cfl = P2Cs(Wl,Model);
cfr = P2Cs(Wr,Model);

%Take min/max of Roe's eigenvalues and L/R states
ar = max( evals(4,:,:) , Wr(2,:,:) + cfr );
al = min( evals(1,:,:) , Wl(2,:,:) - cfl );

bp = max(ar,0.0);
bm = min(al,0.0);

%Calculate L/R fluxes
%Note, fluxes apply to *CONSERVED* not primitive vars
[Fl Fr] = CalcLR_Fluxes(Wl,Wr,bm,bp,Model);

%Compute the HLLE flux at each interface
scl = 0.5*(bp + bm)./(bp - bm);

for n=1:4
    Flux(n,:,:) = 0.5*( Fl(n,:,:) + Fr(n,:,:) ) + ( Fl(n,:,:) - Fr(n,:,:) ).*scl;
end

%Calculates fluxes using the HLLC 3-wave solver
function Flux = RiemannFlux_HLLC(Wl,Wr,Model,varargin)
%Start by doing the same as HLLE


Flux = zeros(size(Wl));

%Calculate Roe averages
[droe vxroe vyroe hroe] = RoeAverages(Wl,Wr,Model);

evals = RoeEigenvals(droe,vxroe,vyroe,hroe,Model);

%Calculate min/max wavespeeds of L/R states
cfl = P2Cs(Wl,Model);
cfr = P2Cs(Wr,Model);

%Take min/max of Roe's eigenvalues and L/R states
ar = max( evals(4,:,:) , Wr(2,:,:) + cfr );
al = min( evals(1,:,:) , Wl(2,:,:) - cfl );

bp = max(ar,0.0);
bm = min(al,0.0);

%Calculate L/R fluxes
%Note, fluxes apply to *CONSERVED* not primitive vars
[Fl Fr] = CalcLR_Fluxes(Wl,Wr,bm,bp,Model);

%Calculate speed of contact wave and pressure
tl = Wl(4,:,:) + Wl(1,:,:).*Wl(2,:,:).* ( Wl(2,:,:) - al );
tr = Wr(4,:,:) + Wr(1,:,:).*Wr(2,:,:).* ( Wr(2,:,:) - ar );

dl = Wl(1,:,:).*Wl(2,:,:) - Wl(1,:,:).*al;
dr = -1*(Wr(1,:,:).*Wr(2,:,:) - Wr(1,:,:).*ar);

am = (tl - tr) ./ (dl + dr); %Contact wave speed
cp = (dl.*tr + dr.*tl) ./ (dl + dr); %Contact pressure at moving surface
cp = max(cp,0.0); %Enforce positivity on pressure

%Calculate scaling terms for L/R fluxes
%Distinguish right-moving contact wave from left
Ind = (am >= 0.0);
SclR = zeros(size(am)); SclL = SclR; SclM = SclR;

%Right-moving contact wave
SclL(Ind) = am(Ind)./(am(Ind)-bm(Ind));
SclR(Ind) = 0.0;
SclM(Ind) = -bm(Ind)./(am(Ind)-bm(Ind));

%Left-moving contact wave
SclL(~Ind) = 0.0;
SclR(~Ind) = -am(~Ind)./(bp(~Ind)-am(~Ind));
SclM(~Ind) = bp(~Ind)./(bp(~Ind)-am(~Ind));

%Compute HLLC flux
for n=1:4
    Flux(n,:,:) = SclL.* Fl(n,:,:) + SclR.*Fr(n,:,:);
end

%Add correction for contact
%Don't correct if you're connected to a solid cell
% if (length(varargin) > 0)
%     Ns = 2;
%     Fi_sol = varargin{1};
%     %Spread Fi_sol to neighboring cells
%     Fi_sm = TwoSmooth(1*Fi_sol,Ns,Ns);
%     Fi_sol = ~(Fi_sm>0);
%     SclM(1,:,:) = squeeze(SclM).*Fi_sol; %Stupid matlab
% end

Flux(2,:,:) = Flux(2,:,:) + SclM.*cp; %Correct Mx w/ contact pressure
Flux(4,:,:) = Flux(4,:,:) + SclM.*cp.*am; %Correct energy

%Flux calculation complete!

%------
%Calculate fluxes along the lines bm/bp according to
%F_L - S_L U_L
%F_R - S_R U_R
function [Fl Fr] = CalcLR_Fluxes(Wl,Wr,bm,bp,Model)

%Create variables
Fl = zeros(size(Wl));
Fr = Fl;

%Mass flux (Mx - bm*d)
Fl(1,:,:) = Wl(1,:,:).*Wl(2,:,:) - bm.*Wl(1,:,:);
Fr(1,:,:) = Wr(1,:,:).*Wr(2,:,:) - bp.*Wr(1,:,:);

%Mx flux ( Mx*(Vx-bm) )
Fl(2,:,:) = Wl(1,:,:).*Wl(2,:,:).*( Wl(2,:,:) - bm );
Fr(2,:,:) = Wr(1,:,:).*Wr(2,:,:).*( Wr(2,:,:) - bp );

%My flux
Fl(3,:,:) = Wl(1,:,:).*Wl(3,:,:).*( Wl(2,:,:) - bm );
Fr(3,:,:) = Wr(1,:,:).*Wr(3,:,:).*( Wr(2,:,:) - bp );

%Add pressure correction in x
Fl(2,:,:) = Fl(2,:,:) + Wl(4,:,:);
Fr(2,:,:) = Fr(2,:,:) + Wr(4,:,:);

%Energy flux [ E*(Vx-bm) + P*Vx ]

Ul = stateP2C(Wl,Model); El = Ul(4,:,:);
Ur = stateP2C(Wr,Model); Er = Ur(4,:,:);


Fl(4,:,:) = El.*( Wl(2,:,:) - bm ) + Wl(4,:,:).*Wl(2,:,:);
Fr(4,:,:) = Er.*( Wr(2,:,:) - bp ) + Wr(4,:,:).*Wr(2,:,:);

function [droe vxroe vyroe hroe] = RoeAverages(Wl,Wr,Model)
Gam = Model.Init.Gam;

%droe = sqrt(d_L)*sqrt(d_R)
%vxroe = sqrt(d_L)*vx_L + sqrt(d_R)*vx_R )/(sqrt(d_L) + sqrt(d_R))


droe = realsqrt(Wl(1,:,:)).*realsqrt(Wr(1,:,:));
invden = 1./( realsqrt(Wl(1,:,:)) + realsqrt(Wr(1,:,:)) );
vxroe = (realsqrt(Wl(1,:,:)).*Wl(2,:,:) + realsqrt(Wr(1,:,:)).*Wr(2,:,:) ).*invden;
vyroe = (realsqrt(Wl(1,:,:)).*Wl(3,:,:) + realsqrt(Wr(1,:,:)).*Wr(3,:,:) ).*invden;

%Calculate L/R enthalpy
hL = CalcEnthalpy(Wl,Model);
hR = CalcEnthalpy(Wr,Model);

hroe = ( realsqrt(Wl(1,:,:)).*hL + realsqrt(Wr(1,:,:)).*hR ) .* invden;

function evals = RoeEigenvals(d,vx,vy,h,Model)

Gam = Model.Init.Gam;
[trash Nx Ny] = size(d);
evals = zeros(4,Nx,Ny);
TINY = 1.0e-12;
vsq = vx.^2 + vy.^2;
%Note, this relies on the assumption of an adiabatic EOS, ie P = (Gam-1)e
asq = (Gam-1)*max(h - 0.5*vsq, TINY);
a = realsqrt(asq);
%Is there an issue w/ singleton dimensions here?
evals(1,:,:) = vx - a;
evals(2,:,:) = vx;
evals(3,:,:) = vx;
evals(4,:,:) = vx + a;


function H = CalcEnthalpy(W,Model)

D = W(1,:,:);
P = W(4,:,:);

U = stateP2C(W,Model);
E = U(4,:,:);

H = (E+P)./D;

%Prim state 2 sound speed
function cs = P2Cs(W,Model)

cs = Prim2Cs( W(1,:,:),W(2,:,:),W(3,:,:),W(4,:,:), Model );

function U = stateP2C(W,Model)
U = zeros(size(W));
[D Mx My E] = Prim2Con( W(1,:,:),W(2,:,:),W(3,:,:),W(4,:,:), Model);
U(1,:,:) = D;
U(2,:,:) = Mx;
U(3,:,:) = My;
U(4,:,:) = E;

function W = stateC2P(U,Model)

W = zeros(size(U));
[D Vx Vy P] = Con2Prim( U(1,:,:),U(2,:,:),U(3,:,:),U(4,:,:), Model);
W(1,:,:) = D;
W(2,:,:) = Vx;
W(3,:,:) = Vy;
W(4,:,:) = P;

