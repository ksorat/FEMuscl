%Takes data from the solid solver and strips it down to data for the gas
%solver

function lvlDef = Solid2Gas(Nodes,Elements,v)

NumNod = length(Nodes);

xv = [];
yv = [];
nid = [];
for n=1:NumNod
    nod = Nodes{n};
    nod = nod.setNeighbors();
    own = length(nod.owners);
    if (own <=3)
        xn = nod.xyz(1);
        yn = nod.xyz(2);
        %Less than 3 owners, so exterior
        xv = [xv xn];
        yv = [yv yn];
        nid = [nid n];
    end
end

%Now need to get oriented version (will go with ccw)
Nbdry = length(xv);

[xvc yvc I] = Orient(xv,yv);
[xvc yvc] = poly2cw(xvc,yvc);

% %NOTE: THIS ONLY WORKS FOR CONVEX SHAPES

% xcm = sum(xv)/Nbdry;
% ycm = sum(yv)/Nbdry;
% xvo = xv-xcm; yvo = yv-ycm;
% 
% theta = atan2(yvo,xvo);
% 
% [ths I] = sort(theta);
% xvc = xv(I);
% yvc = yv(I);

%Currently assuming 1 object
lvlDef.numObs = 1;
obsDat.xv = xvc;
obsDat.yv = yvc;
obsDat.vx = zeros(1,Nbdry);
obsDat.vy = zeros(1,Nbdry);
obsDat.vxid = zeros(1,Nbdry);
obsDat.vyid = zeros(1,Nbdry);

obsDat.nid = nid(I);

for n=1:Nbdry
    nnod = obsDat.nid(n);
    nod = Nodes{nnod};
    vxid = nod.DOFids(1);
    vyid = nod.DOFids(2);
    if (vxid ~= 0)  
        obsDat.vx(n) = v(vxid);
    else
        obsDat.vx(n) = 0;
    end
    if (vyid ~= 0)
        obsDat.vy(n) = v(vyid);
    else
        obsDat.vy(n) = 0.0;
    end
    %Save return mapping
    obsDat.vxid(n) = vxid;
    obsDat.vyid(n) = vyid;
end

lvlDef.obsDat{1} = obsDat;

%Close polygon


% Ns = 1000;
% x = linspace(0,200,Ns); y = linspace(0,200,Ns);
% [yy xx] = meshgrid(y,x);
% In = inpolygon(xx,yy,xvc,yvc);
% plot(obsDat.xv,obsDat.yv,'bo'); hold on;
% quiver(obsDat.xv,obsDat.yv,obsDat.vx,obsDat.vy,'r'); 
% axis equal
% hold off;
% 
% keyboard

%Takes collection of points (xv,yv) and returns the same points but ordered
%by distance to previous point
function [xvc yvc I] = Orient(xv,yv)
Nv = length(xv);
xvc = zeros(1,Nv);
yvc = zeros(1,Nv);
Used = false(1,Nv);
I = zeros(1,Nv);

%First point is arbitrary
xvc(1) = xv(1);
yvc(1) = yv(1);
I(1) = 1;
Used(1) = true;

for n=2:Nv
    xp = xvc(n-1);
    yp = yvc(n-1);
    dv = sqrt( (xp-xv).^2 + (yp-yv).^2 );
    dv(Used) = Inf; %Don't repick one
    [dvmin i] = min(dv);
    xvc(n) = xv(i);
    yvc(n) = yv(i);
    Used(i) = true;
    I(n) = i;
end


