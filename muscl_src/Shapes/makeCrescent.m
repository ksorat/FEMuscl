
%Creates a crescent shape and returns the xv/yv of the discretization
%Uses DelR*DelTh/N as typical length scale
function [xv yv] = makeCrescent(InP, N)

r0 = InP.R(1); r1 = InP.R(2);
th0 = InP.Theta(1); th1 = InP.Theta(2); %Assuming degrees

Th = linspace(th0,th1,N);

%Start by constructing inward arc, counter clock
xv_in = r0*cosd(Th);
yv_in = r0*sind(Th);

%Calculate effective discretization (ds)
DelTh_Rad = (th1-th0)*(pi)/180; %angle of arc in radiations
ds = r0*DelTh_Rad/N;
%Move to outer arc
ArcN = (r1-r0)/ds;
ceil(ArcN);

rArc = linspace(r0,r1,ArcN);

xv_arc = rArc*cosd(th1);
yv_arc = rArc*sind(th1);

OutN = ceil(r1*DelTh_Rad/ds);
Thp = linspace(th1,th0,OutN);
xv_out = r1*cosd(Thp);
yv_out = r1*sind(Thp);

rArcP = linspace(r1,r0,ArcN);
xv_arcp = rArcP*cosd(th0);
yv_arcp = rArcP*sind(th0);

% plot(xv_in,yv_in,'ro'); hold on;
% plot(xv_arc,yv_arc,'b+'); 
% plot(xv_out,yv_out,'ro');
% plot(xv_arcp,yv_arcp,'b+'); 
% 
% axis equal; hold off

xv = [ xv_in xv_arc(2:end) xv_out(2:end) xv_arcp(2:end) ];
yv = [ yv_in yv_arc(2:end) yv_out(2:end) yv_arcp(2:end) ];