%% FEM Engine
% Main FEM program calling all objects and subroutines

%% Preliminaries

close all
clear all
clear classes
clc

global Parts Nodes Elements Materials ElForms Damage Lengths Speeds GDOF GDOFall u v Fext Fmax M dim a t;

restart = false; % true or false for restarting a run
endtime = 0.5; % solver endtime (FLOAT)
tsScale = 0.9; % scale for internally calculated timestep (90% or lower of the calculated CFL condition is a good safety factor) (FLOAT)
filename = 'beam.k'; % input deck file name (STRING)

if restart
    
    load restart.mat
    Solver.run(tsScale,endtime); % run the solver
    
else

    %% Initialize Global Variables
    % It probably isn't necessary to create so many "global" variables. I did
    % it for lazy convenience. Parts, Nodes, Elements, Materials, and Elforms are
    % arrays of objects, so those are definitely global. Lengths is a global
    % vector I created that each element object populates with its smallest leg
    % length (used to calculated the timestep to satisfy the Courant
    % condition). Speeds is a global vector that holds the sound speeds of each
    % different material object created. In the end the minimum "Length" and
    % maximum "Speed" are used to calculate the solver timestep.    

    dim = 2; % model dimensions (either 2 or 3) (INT)
    GDOF = 0; % total global degrees of freedom. Could work well as a static INT OR LONG (this can get to be a super big number in some models) for the node object
    GDOFall = 0; % total global degrees of freedom for the stress projection. Could work well as a static INT OR LONG (this can get to be a super big number in some models) for the node object


    %% Setup Part, Material, and Element Profiles

    Parts = cell(1,1); % open up a Parts array
    Materials = cell(1,1); % open up a Materials array
    ElForms = cell(1,1); % open up an Element Formulation array
    Damage = {}; % open up an Element Formulation array
    Speeds = zeros(length(Materials),1); % initialize the sound speeds array
    ElForms{1} = Q4; % Creating the standard Q4 element
    Damage{1} = maxPrincipalDamage(.0008);
    Materials{1} = LinearElasticTL_planeStress(1,1.1e-6,.05,.45); % create material model (ID,density,modulus,poisson)
    
    Parts{1} = Part(Materials{1},ElForms{1}); % create part object
    Parts{1}.damage = 1;
    for i = 1:length(Speeds)

        if length(Speeds) == 1

            Speeds = Materials{1}.c; % c is a float held within each Materials object
            break

        end

        Speeds(i) = Materials{i}.c;

    end % populate sound speed array for timestep calcs

    %% Load Input Deck

    [~,~,~,~,~] = inputdeckreader(filename); % read LS-DYNA input deck (old, lots of uselessness inside of it, but it works for now)

    %% Initialize Time Integration Variables

    M = diag(M); % turn M matrix into an M vector for 2x speedup of a=F/M calculation (M matrix is diagonal anyway due to the way we evaluated it) (a GDOF x 1 VECTOR of FLOATS)
    u = zeros(GDOF,1); % displacement (a GDOF x 1 VECTOR of FLOATS) 
    v = zeros(GDOF,1); % velocity (a GDOF x 1 VECTOR of FLOATS)
    Fext = zeros(GDOF,1); % external force (a GDOF x 1 VECTOR of FLOATS)
    a = zeros(GDOF,1); % acceleration (a GDOF x 1 VECTOR of FLOATS)
    Fmax = .01; % maximum force at the end of the explicit simulation (linearly interpolated from 0 to endtime for now) -- this was to make my beam model work. Need a fancier force input than this to pass LBM data
    Solver = explicitSolverTF(); % set implicit or explicit solver

    %% Run Time Loop

    Solver.run(tsScale,endtime); % run the solver
    %plotter2D % plot the end deformed condition
end

%% Stress Plot
element.NodeStress();
figure
for i = 1:length(Elements)
    XX = [Elements{i}.xyz(1,1) Elements{i}.xyz(2,1) Elements{i}.xyz(3,1) Elements{i}.xyz(4,1) Elements{i}.xyz(1,1)];
    YY = [Elements{i}.xyz(1,2) Elements{i}.xyz(2,2) Elements{i}.xyz(3,2) Elements{i}.xyz(4,2) Elements{i}.xyz(1,2)];
%    cc = [Elements{i}.sig{1}(2,2) Elements{i}.sig{2}(2,2) Elements{i}.sig{3}(2,2) Elements{i}.sig{4}(2,2) Elements{i}.sig{1}(2,2)];
%    cc = [sig_nodes(Elements{i}.nodes(1)*3-1) sig_nodes(Elements{i}.nodes(2)*3-1) sig_nodes(Elements{i}.nodes(3)*3-1) sig_nodes(Elements{i}.nodes(4)*3-1) sig_nodes(Elements{i}.nodes(1)*3-1)];
    cc = [Nodes{Elements{i}.nodes(1)}.sig(2) Nodes{Elements{i}.nodes(2)}.sig(2) Nodes{Elements{i}.nodes(3)}.sig(2) Nodes{Elements{i}.nodes(4)}.sig(2) Nodes{Elements{i}.nodes(1)}.sig(2)];
    patch(XX,YY,cc)
end
axis equal
colorbar