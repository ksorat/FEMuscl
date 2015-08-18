function [xv yv] = naca(inP, N)
%Returns a polygon-ization of a naca airfoil w/ N points
%Uses data from inP
%inP.x0, inP.y0
%inP.T 
%inP.c

disp('NACA currently broken, go away'); pause;
x0 = inP.x0; y0 = inP.y0;
T = inP.T; c = inP.c;

if isfield(inP,'alpha')
    alpha = -1*inP.alpha; %Rotation in degrees
else
    alpha = 0;
end

if isfield(inP,'m')
    m = inP.m;
else
    m = 0;
end

if isfield(inP,'p')
    p = inP.p;
else
    p = c;
end

Lam = linspace(0,c,N);

a0= 0.2969;
a1=-0.1260;
a2=-0.3516;
a3= 0.2843;
a4=-0.1036;

Lscl = Lam/c;
Yt = 5*T*c*( a0*sqrt(Lscl) + a1*Lscl + a2*Lscl.^2 + a3*Lscl.^3 + a4*Lscl.^4);

Ind = (Lscl >= 0) & (Lscl <= p);
Ind1 = ~Ind;

Yc = zeros(size(Yt));

Yc(Ind) = m*(Lam(Ind)/(p*p)).*( 2*p - Lscl(Ind));
Yc(Ind1) = m*(c-Lam(Ind1)).*(1 + Lscl(Ind1) - 2*p)/ ( (1-p)^2 );

atheta = zeros(size(Yt));
atheta(Ind) = 2*m*(p - Lscl(Ind))/(p*p);
atheta(Ind1) = 2*m*(p-Lscl(Ind1))/( (1-p)^2 );

theta = atan(atheta);
xU = Lam-Yt.*sin(theta);
yU = Yc + Yt.*cos(theta);

xL = Lam+Yt.*sin(theta);
yL = Yc - Yt.*cos(theta);

%LamEps = 0.01;

% Lam1 = linspace(0,LamEps,N/3);
% Lam2 = linspace(LamEps,2-LamEps,N/3);
% Lam3 = linspace(2-LamEps,2,N/3);
% Lam = [Lam1(1:end-1) Lam2(1:end-1) Lam3];



xv = zeros(1,N); yv = zeros(1,N);



%Old stuff below here
% 
% for i=1:length(Lam)
%     lambda = Lam(i);
%     if (lambda <= 1)
%         xv(i) = lambda*c + x0;
%         pos = 1; 
%     else
%         %Lambda >1 & < 2
%         xv(i) = x0 + c*(2-lambda);
%         pos = -1; 
%     end
%     xscl = (xv(i) - x0)/c;
%     yv(i) = y0 + pos*(5*c*T)*(a0*sqrt(xscl) + a1*xscl + a2*(xscl^2) + a3*(xscl^3) + a4*(xscl^4));
% end
% 
% %Now rotate
% C = cosd(alpha); S = sind(alpha);
% 
% %Move so that origin is at center (c/2,0) of airfoil
% xvo = xv -  x0 - 0.5*c;
% yvo = yv-y0;
% 
% %Rotate
% xvp = C*xvo - S*yvo;
% yvp = S*xvo + C*yvo;
% 
% %Translate back
% xv = xvp+0.5*c + x0; yv = yvp+y0;
% 
