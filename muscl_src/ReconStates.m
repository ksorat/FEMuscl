
%--------
%Calculates directed states within a cell, ie
%Density @ North/South/East/West interface
%according to a given method/monotonization routine
function [nW eW sW wW] = ReconStates(Grid,Gas,method)

%Construct a holder regardless of method
%Each holder is of size [Nv,Nx,Ny]
[Nx,Ny] = size(Gas.D);
nW = zeros(4,Nx,Ny);


switch lower(method)
    case {'pcm'}
        nW(1,:,:) = Gas.D;
        nW(2,:,:) = Gas.Vx;
        nW(3,:,:) = Gas.Vy;
        nW(4,:,:) = Gas.P;
        eW = nW;
        sW = nW;
        wW = nW;
    case {'plm'}
        %Create other holders
        eW = nW; sW = nW; wW = nW;
        [wW eW] = DirstatePLM(Grid,Gas,'x');
        [sW nW] = DirstatePLM(Grid,Gas,'y');
    case {'ppm'}
        %Create other holders
        eW = nW; sW = nW; wW = nW;
        [wW eW] = DirstatePPM(Grid,Gas,'x');
        [sW nW] = DirstatePPM(Grid,Gas,'y');
        
        
    otherwise
        disp('Unknown reconstruction method');
end

nW = Checkstate(nW);
eW = Checkstate(eW);
sW = Checkstate(sW);
wW = Checkstate(wW);

%----------
%Calculates minus/plus states (in dirvec direction) for all variables using
%the PPM method
function [mW pW] = DirstatePPMvec(Grid,Gas,dir)

%Setup bounds
[Nx,Ny] = size(Gas.D);
allgas = zeros(Gas.Nv,Nx,Ny);
allgas(1,:,:) = Gas.D;
allgas(2,:,:) = Gas.Vx;
allgas(3,:,:) = Gas.Vy;
allgas(4,:,:) = Gas.P;
mW = allgas; pW = allgas;

IntX = 2:Nx-1;
IntY = 2:Ny-1;
switch lower(dir)
    case{'x'}
        Xp = 3:Nx;
        Yp = IntY;
        Xm = 1:Nx-2;
        Ym = IntY;
    case{'y'}
        Xp = IntX;
        Yp = 3:Ny;
        Xm = IntX;
        Ym = 1:Ny-2;
end

[DelWl DelWr DelWc] = GetDels(allgas,IntX,IntY,Xm,Xp,Ym,Yp);
DelWmon = Monotonize(DelWl,DelWr,DelWc);

mW(:,IntX,IntY) = 0.5*( allgas(:,IntX,IntY) + allgas(:,Xm,Ym) ) - (1/6)* ( DelWmon(:,IntX,IntY) + DelWmon(:,Xm,Ym) );
pW(:,IntX,IntY) = 0.5*( allgas(:,Xp,Yp) + allgas(:,IntX,IntY) ) - (1/6)* ( DelWmon(:,Xp,Yp) + DelWmon(:,IntX,IntY) );


%Apply further monotonicity constraints
Ind = ( (pW(:,:,:) - allgas(:,:,:)).*(mW(:,:,:) - allgas(:,:,:)) <= 0);
mW(Ind) = allgas(Ind);
pW(Ind) = allgas(Ind);

Lam = 6*( pW(:,:,:) - mW(:,:,:) ) .* ( allgas(:,:,:) - 0.5*(pW(:,:,:)+mW(:,:,:)));
Beta = (pW(:,:,:) - mW(:,:,:)).^2;
   
Ind = (Lam > Beta);
mW(Ind) = 3*allgas(Ind) - 2*pW(Ind);
   
Ind = (Lam < -1*Beta);
pW(Ind) = 3*allgas(Ind) - 2*mW(Ind);


%Monotonize the deltas
function DelWmon = Monotonize(DelWl,DelWr,DelWc)

tmpMin = min(2*abs(DelWl),2*abs(DelWr));
fullMin = min(tmpMin,abs(DelWc));
DelWmon = sign(DelWc).*fullMin;

%Calculate L/R/C deltas
function [DelWl DelWr DelWc] = GetDels(allgas,IntX,IntY,Xm,Xp,Ym,Yp)
[Nv Nx Ny] = size(allgas);

DelWl = zeros(Nv,Nx,Ny);
DelWr = zeros(Nv,Nx,Ny);
DelWc = zeros(Nv,Nx,Ny);

DelWl(:,IntX,IntY) = allgas(:,IntX,IntY) - allgas(:,  Xm,  Ym);
DelWr(:,IntX,IntY) = allgas(:,  Xp,  Yp) - allgas(:,IntX,IntY);
DelWc(:,IntX,IntY) = 0.5*( allgas(:,Xp,Yp) - allgas(:,Xm,Ym) );

function [mW pW] = DirstatePPM(Grid,Gas,dir)

%Setup bounds 
[Nx,Ny] = size(Gas.D);
allgas = zeros(Gas.Nv,Nx,Ny);
allgas(1,:,:) = Gas.D;
allgas(2,:,:) = Gas.Vx;
allgas(3,:,:) = Gas.Vy;
allgas(4,:,:) = Gas.P;
mW = allgas; pW = allgas;

IntX = 2:Nx-1;
IntY = 2:Ny-1;
switch lower(dir)
    case{'x'}
        Xp = 3:Nx;
        Yp = IntY;
        Xm = 1:Nx-2;
        Ym = IntY;
    case{'y'}
        Xp = IntX;
        Yp = 3:Ny;
        Xm = IntX;
        Ym = 1:Ny-2;
end

%Note, this is somewhat poorly coded as it replicates parts of PLM

DelWl = zeros(Nx,Ny);
DelWr = zeros(Nx,Ny);
DelWc = zeros(Nx,Ny);
DelWmon = zeros(1,Nx,Ny); %Stupid matlab

for n=1:Gas.Nv
   %Calculate directed states along direction 'dir' for each variable
   DelWl(IntX,IntY) = allgas(n,IntX,IntY) - allgas(n,Xm,Ym);
   DelWr(IntX,IntY) = allgas(n,Xp,Yp) - allgas(n,IntX,IntY);
   DelWc(IntX,IntY) = 0.5*( allgas(n,Xp,Yp) - allgas(n,Xm,Ym) );
   
   %Calculate minima
   tmpmin = min(2*abs(DelWl),2*abs(DelWr));
   fullmin = min(tmpmin,abs(DelWc));
   DelWmon(1,:,:) = sign(DelWc).*fullmin; %Monotonized differences
   
   mW(n,IntX,IntY) = 0.5*( allgas(n,IntX,IntY) + allgas(n,Xm,Ym) ) - (1/6)* ( DelWmon(1,IntX,IntY) + DelWmon(1,Xm,Ym) );
   pW(n,IntX,IntY) = 0.5*( allgas(n,Xp,Yp) + allgas(n,IntX,IntY) ) - (1/6)* ( DelWmon(1,Xp,Yp) + DelWmon(1,IntX,IntY) );
   
   %Apply further monotonicity constraints
   Ind = ( (pW(n,:,:) - allgas(n,:,:)).*(mW(n,:,:) - allgas(n,:,:)) <= 0);
   mW(n,Ind) = allgas(n,Ind);
   pW(n,Ind) = allgas(n,Ind);
   
   Lam = 6*( pW(n,:,:) - mW(n,:,:) ) .* ( allgas(n,:,:) - 0.5*(pW(n,:,:)+mW(n,:,:)));
   Beta = (pW(n,:,:) - mW(n,:,:)).^2;
   
   Ind = (Lam > Beta);
   mW(n,Ind) = 3*allgas(n,Ind) - 2*pW(n,Ind);
   
   Ind = (Lam < -1*Beta);
   pW(n,Ind) = 3*allgas(n,Ind) - 2*mW(n,Ind);
   
end


%----------
%Calculates minus/plus states (in dirvec direction) for all variables using
%the PLM method
function [mW pW DelWmon] = DirstatePLM(Grid,Gas,dir)

[Nx,Ny] = size(Gas.D);
allgas = zeros(Gas.Nv,Nx,Ny);
allgas(1,:,:) = Gas.D;
allgas(2,:,:) = Gas.Vx;
allgas(3,:,:) = Gas.Vy;
allgas(4,:,:) = Gas.P;

mW = allgas; pW = allgas;

IntX = 2:Nx-1;
IntY = 2:Ny-1;
switch lower(dir)
    case{'x'}
        Xp = 3:Nx;
        Yp = IntY;
        Xm = 1:Nx-2;
        Ym = IntY;
    case{'y'}
        Xp = IntX;
        Yp = 3:Ny;
        Xm = IntX;
        Ym = 1:Ny-2;
end

for n=1:Gas.Nv
   %Calculate directed states along direction 'dir' for each variable
   DelWl = allgas(n,IntX,IntY) - allgas(n,Xm,Ym);
   DelWr = allgas(n,Xp,Yp) - allgas(n,IntX,IntY);
   DelWc = 0.5*( allgas(n,Xp,Yp) - allgas(n,Xm,Ym) );
   
   %Calculate minima
   tmpmin = min(2*abs(DelWl),2*abs(DelWr));
   fullmin = min(tmpmin,abs(DelWc));
   DelWmon = sign(DelWc).*fullmin; %Monotonized differences
   
   
   mW(n,IntX,IntY) = allgas(n,IntX,IntY) - 0.5*DelWmon;
   pW(n,IntX,IntY) = allgas(n,IntX,IntY) + 0.5*DelWmon;
end
        
function W = Checkstate(W)

global SMALL_NUM;

W(1,:,:) = max( W(1,:,:), SMALL_NUM );
W(4,:,:) = max( W(4,:,:), SMALL_NUM );


