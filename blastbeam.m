%Test of blast (gas) + beam (solid) coupling

clear; close all;
Model.recon = 'plm'; Model.solver = 'hll';

Low = true;
%Generic output
Pic.view = true;
Pic.val = 'd';
Pic.pg = false;
Pic.dovid = false;
Pic.cax = [0.0 5.0e-6];

%generic initialization
Init.rho0 = 1e-6;
Init.P0 = 1e-6;
Init.DelP = 1000;
Model.Tfin = 20;
if (Low)
    Model.Bds = [-40 101 -10 175];
    Model.Nvec = round( [768 1024]/8);
    Init.filename = 'beam_low_res.k';
    Init.rad = 10;
    Init.cent = [50 60];
else
   Model.Bds = [50 300 -50 200];
    Model.Nvec = round( [512 1024]/16);
   Init.filename = 'beam.k';
end
Init.problem = 'blast';


Model.bcs.ibx = 'outflow'; Model.bcs.obx = 'outflow';
Model.bcs.iby = 'outflow'; Model.bcs.oby = 'outflow';


Model.Init = Init; Model.Pic = Pic;
[Grid Gas Nodes Elements] = runjoint(Model);
