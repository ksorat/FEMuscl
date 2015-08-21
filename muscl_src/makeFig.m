function makeFig(Model,Grid,Gas)

global Nfig;
global Elements;
global Nodes;


Gam = Model.Init.Gam;
%Set bounds
if (Model.Pic.pg)
    %Print ghosts
    is = Grid.isd;
    ie = Grid.ied;
    js = Grid.jsd;
    je = Grid.jed;
else
    %Don't print ghosts
    is = Grid.is;
    ie = Grid.ie;
    js = Grid.js;
    je = Grid.je;
end

x = Grid.xc(is:ie);
y = Grid.yc(js:je);

Pos = false; %Is this a positive definite quantity?

switch lower(Model.Pic.val)
    case{'d'}
        Z = Gas.D;
        varS = 'Density';
        Pos = true;
    case{'p'}
        Z = Gas.P;
        varS = 'Pressure';
        Pos = true;
    case{'logp'}
        Z = log10(Gas.P);
        varS = 'Log Pressure';
    case{'vx'}
        Z = Gas.Vx;
        varS = 'X-Velocity';
    case{'vy'}
        Z = Gas.Vy;
        varS = 'Y-Velocity';
    case{'spd'}
        Z = sqrt( Gas.Vx.^2 + Gas.Vy.^2);
        varS = 'Speed';
        Pos = true;
    case{'k'}
        Z = 0.5*Gas.D.*( Gas.Vx.^2 + Gas.Vy.^2);
        varS = 'Kinetic Energy';
        Pos = true;
    case{'ma','mach'}
        Spd = sqrt( Gas.Vx.^2 + Gas.Vy.^2);
        Cs = Prim2Cs(Gas.D,Gas.Vx,Gas.Vy,Gas.P,Model);
        Z = Spd./Cs;
        varS = 'Mach Number';
        Pos = true;
    case{'sd'}
        Z = Grid.lvlSet.sd;
        varS = 'Signed Distance';
    case{'in','inside'}
        Z = 1*( Grid.lvlSet.sd <= 0);
        varS = 'Inside Object';
    otherwise
        disp('Unknown diagnostic');
        pause;
end

Zp = Z(is:ie,js:je);
kcolor(x,y,Zp'); %axis equal;
%set(gcf,'units','normalized','outerposition',[0 0 0.5 0.75]);
Xdom = Grid.xi(end)-Grid.xi(1);
Ydom = Grid.yi(end)-Grid.yi(1);
aspRat = Xdom/Ydom;
FixIn = 18;
set(gcf,'units','inches','outerposition', [ 0 0 FixIn FixIn/aspRat]); 

%axis equal;
if isfield(Model.Pic,'cax')
    %caxis(Model.Pic.cax);
    cAx = Model.Pic.cax;
else
    nStd = 3;
    zM = mean(Zp(:));
    zStd = std(Zp(:));
    cAx = [ (zM-nStd*zStd) (zM+nStd*zStd) ];
end
if (Pos)
    cAx(1) = 0.0;
end


%Draw solid info
hold on
for i=1:length(Elements)
    if isempty(Elements{i}) continue; end
    xel = Elements{i}.xyz(:,1);
    yel = Elements{i}.xyz(:,2);
    fill(xel,yel,'w');
end
hold off
% if (Model.lvlSet.present)
%     sd = Grid.lvlSet.sd(is:ie,js:je);
%     hold on;
%     
%     if (Model.Pic.pg)
%         
%         for n=1:Model.Init.lvlDef.numObs
%             obsDat = Model.Init.lvlDef.obsDat{n};
%             
%             plot(obsDat.xv,obsDat.yv,'wo');
%         end
%         
%         
%     else
%         %Fill region
%         for n=1:Model.Init.lvlDef.numObs
%             obsDat = Model.Init.lvlDef.obsDat{n};
%             fill(obsDat.xv,obsDat.yv,'w');
%             %plot(obsDat.xv,obsDat.yv,'wo');
%             Nv = length(obsDat.xv);
%             xvr = obsDat.xv(1:Nv-1); yvr = obsDat.yv(1:Nv-1);
%             
%             x1 = sum(xvr)/(Nv-1);
%             y1 = sum(yvr)/(Nv-1);
%             
%             x0 = obsDat.xv(1);
%             y0 = obsDat.yv(1);
%             
%             %plot([x0 x1],[y0 y1],'k');
%         end
%         
%     end
%     
%     hold off;
% end


xlabel('X'); ylabel('Y');
titS = sprintf('%s @ t=%3.3f', varS, Grid.t);
title(titS);

axis equal
colorbar();
caxis(cAx);
drawnow;
if (Model.Pic.dovid)
    Figfile = sprintf('%s/Vid.%04d.png', Model.Pic.vid_dir,Nfig);
    export_fig(Figfile);
end

Nfig = Nfig+1;

function FillPolys(x,y,sd)

C = contourc(x,y,sd',[0 0]);

N = length(C);
ic = 1;
while (ic < N)
    icp = ic+1;
    num = C(2,ic); %How many elements are in this
    xobs = C(1,icp:icp+num-1);
    yobs = C(2,icp:icp+num-1);
    ic = icp+num;
    fill(xobs,yobs,'w');
end