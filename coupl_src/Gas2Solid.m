%Takes in level set data (from gas solver) and plugs this into Fext
%using forces calculated from gas

function Fext = Gas2Solid(lvlDef)

global GDOF;
Fext = zeros(GDOF,1);

Nobj = lvlDef.numObs;

for n=1:Nobj
    obs = lvlDef.obsDat{n};
    Nv = length(obs.xv); %Number of vertices
    for k=1:Nv-1 %Because last point is redundant
        vxid = obs.vxid(k);
        vyid = obs.vyid(k);
        if (vxid ~= 0)
            Fext(vxid) = obs.Fx(k);
        end
        if (vyid ~= 0)
            Fext(vyid) = obs.Fy(k);
        end
    end
end
