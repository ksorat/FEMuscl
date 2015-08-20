
%Test of jet (gas) + beam (solid) coupling

clear; close all;
Model.recon = 'ppm'; Model.solver = 'hll';
Model.tsDiag = 10;
Init.DoDamage = true;

Low = true;
%Generic output
Pic.view = true;
Pic.val = 'd';
Pic.pg = false;
Pic.dovid = false;
Pic.cax = [0 5.0e-6];

%generic initialization
Init.rho0 = 1e-6;
Init.P0 = 1e-6;
Init.DelP = 1;
Model.Tfin = 100;
if (Low)
    Model.Bds = [-39 59 -10 175];
    Model.Bds = [-30 59 -10 175];
    Model.Nvec = round( [512 1024]/16);
    Init.filename = 'beam_low_res.k';
else
   Model.Bds = [50 300 -50 200];
   Model.Bds = [130 200 -50 200];
    Model.Nvec = round( [512 1024]/16);
   Init.filename = 'beam.k';
end
Init.problem = 'flow';


Init.Min = 50;
Init.cent = 100; Init.rad = 20; Init.disc = false;
Model.bcs.ibx = 'injet'; Model.bcs.obx = 'outflow';
Model.bcs.iby = 'outflow'; Model.bcs.oby = 'outflow';
Cs = sqrt( (5/3)*Init.P0/Init.rho0 ); Init.vin = Init.Min*Cs;

Model.Init = Init; Model.Pic = Pic;
[Grid Gas Nodes Elements] = runjoint(Model);

