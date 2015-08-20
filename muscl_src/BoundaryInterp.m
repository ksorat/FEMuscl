%Interpolates values of Q (defined on xv/yv) onto xb/yb (Qi)

function Qi = BoundaryInterp(xb,yb,xv,yv,Q)

Nv = length(xv);
Np = length(xb);
Qi = zeros(1,Np);

%Remove redundant vertex at end
xv = xv(1:Nv-1);
yv = yv(1:Nv-1);

for n=1:Np %Loop over points (should vectorize)
    %The interpolation here is crude, could do better
    %Ie, PPM
    xp = xb(n);
    yp = yb(n);
    dv = sqrt( (xv-xp).^2 + (yv-yp).^2 ); %Distance to each vertex
    [dvmins,I] = sort(dv); 
    i1 = I(1); i2 = I(2); %Our point (xp,yp) is on line joining x1,y1->x2,y2
    x1 = xv(i1); y1 = yv(i1);
    x2 = xv(i2); y2 = yv(i2);
    Q1 = Q(i1); Q2 = Q(i2);
    
    L = sqrt( (x2-x1)^2 + (y2-y1)^2 );
    L1 = sqrt( (xp-x1)^2 + (yp-y1)^2 );
    m = (Q2-Q1)/L; 
    
    if (L < eps)
        Qi(n) = 0.5*(Q1+Q2);
    else
        Qi(n) = Q1 + m*L1;
    end

end