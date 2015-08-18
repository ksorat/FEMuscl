
function [xv yv] = makeCircle(InP,N)

r0 = InP.R;
x0 = InP.C(1);
y0 = InP.C(2);

tau = linspace(0,2*pi,N)';

xv = x0 + r0*cos(tau);
yv = y0 + r0*sin(tau);