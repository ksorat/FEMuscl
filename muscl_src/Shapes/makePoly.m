%Takes a set of vertices (xv/yv) in reduced form (only as many points as
%vertices) and densifies it (ie, adds interpolated points between vertices)

function [xvd yvd] = makePoly(xv,yv,ds)

%Close polygon if necessary
if ((xv(1) ~= xv(end)) || (yv(1) ~= yv(end))) %Maybe change to eps comparison?
    xv = [xv xv(1)];
    yv = [yv yv(1)];
end

Nv = length(xv);

K = 1;
for n=1:Nv-1
    x0 = xv(n); y0 = yv(n);
    x1 = xv(n+1); y1 = yv(n+1);
    D = sqrt( (x1-x0).^2 + (y1-y0).^2 );
    Num = round(D/ds);
    
    xd = linspace(x0,x1,Num); yd = linspace(y0,y1,Num);
    xvd(K:K+Num-1) = xd; yvd(K:K+Num-1) = yd;
    K = K+Num;
end