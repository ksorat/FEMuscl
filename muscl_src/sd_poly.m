%Calculates distance of input points (x,y) from polygon (xv,yv)

function [sd xb yb Nx Ny] = sd_poly(x,y,xv,yv,Grid)

%Close polygon if necessary
if ((xv(1) ~= xv(end)) || (yv(1) ~= yv(end))) %Maybe change to eps comparison?
    xv = [xv ; xv(1)];
    yv = [yv ; yv(1)];
end

Nv = length(xv); 
Np = length(x);

%Create data structures
sd = zeros(1,Np);
xb = sd; yb = sd; Nx = sd; Ny = sd;

[In On] = inpolygon(x,y,xv,yv);

%Find parameters for segments
%Ax + By + C = 0
A = -diff(yv); B = diff(xv);
C = yv(2:Nv).*xv(1:Nv-1) - xv(2:Nv).*yv(1:Nv-1);
AB = 1./( A.^2 + B.^2); %For projections

%Test for vert/horiz segments
id_v = find( diff(xv) == 0 );
id_h = find( diff(yv) == 0 );


for n=1:Np %Loop over relevant points, maybe switch to vectorized op?
    xn = x(n); yn = y(n);
    %Project this point onto each line of the polygon
    vv = (A*xn + B*yn + C);
    xp = xn - (A.*AB).*vv;
    yp = yn - (B.*AB).*vv;
    %Fix for horiz/vertical
    xp(id_v) = xv(id_v);
    yp(id_h) = yv(id_h);
    
    %Find which projections actually lie on line segment
    idx = (((xp>=xv(1:Nv-1)) & (xp<=xv(2:Nv))) | ((xp>=xv(2:Nv)) & (xp<=xv(1:Nv-1))));
    idy = (((yp>=yv(1:Nv-1)) & (yp<=yv(2:Nv))) | ((yp>=yv(2:Nv)) & (yp<=yv(1:Nv-1))));
    Id = idx & idy;
    
    %Find distance between this point and each vertex
    dv = sqrt((xv(1:Nv-1)-xn).^2 + (yv(1:Nv-1)-yn).^2);
    
    if ( ~any(Id) ) %None of the projections hit the polygon segments
        [d,I] = min(dv);
        x_poly = xv(I);
        y_poly = yv(I);
    else
        %You hit a rib, now use that projection
        dproj = sqrt( (xp(Id) - xn).^2 + (yp(Id) - yn).^2);
        [min_dv,I1] = min(dv);
        [min_dp,I2] = min(dproj);
        [d,I] = min( [min_dv min_dp] );
        if (I == 1)
            x_poly = xv(I1);
            y_poly = yv(I1);
        else %I == 2
            idxs = find(Id);
            x_poly = xp(idxs(I2));
            y_poly = yp(idxs(I2));
        end
    end
        
    sgn = 1;
    if In(n)
        sgn = -1;
    end
    sd(n) = sgn*d; %Get signed distance
    xb(n) = x_poly;
    yb(n) = y_poly;
    
    if On(n)
        %Add this, find two closest vertices and calculate normal to line
        %connecting them (need to figure out orientation)
        %If you *are* a vertex, average normals from 2/3 nearest
        disp('Point is on poly vertex, not yet supported');
        keyboard
    else
        %Calculate normal vector pointing from (xn,yn) -> (xb,yb)
        Px = (x_poly-xn); Py = (y_poly-yn);
        nvec = sqrt(Px^2 + Py^2);
        Nx(n) = Px/nvec; Ny(n) = Py/nvec;
    end

end




