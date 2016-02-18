clear variables;
close all;

% define this problem as static
p.mode = 'static';




%% set material parameters
% parameters for material 1 (St for steel)
p.ESt = 200e6; % E-Modul / Young's modulus / [kN/m^2]
p.nuSt = 0.3; % Querdehnzahl / poissons ratio / [-]
p.rhoSt = 7.85; % Dichte / density / [g/cm^3]

% parameters for material 1 (Rub for rubber)
p.ERub = 0.2e6; % E-Modul / Young's modulus / [kN/m^2]
p.nuRub = 0.4; % Querdehnzahl / poissons ratio / [-]
p.rhoRub = 0.92; % Dichte / density / [g/cm^3]




%% feti solver parameters

% select type of coarse grid
% options: 0: no coarse grid (not possible / wont work for static problems)
% options: 1: use RBM (Rigid Body Modes)
% options: 2: use GenEO, see mode selection
p.Coarse = 1;

% select which geneo modes to use for GenEO coarse grid
% e.g. [1 4 7]: use the first, fourth and seventh GenEO mode to build the coarse space
p.GeneoModes = [1 2 3];

% select scaling type
% options: 'eye'    (no scaling)
% options: 'multi'  (multiplicity)
% options: 'K'      (k-scaling)
% options: 'SchurK' (like k-scaling, using the schur-complement of K instead of K)
% options: 'D'      (only for dynamic, like k-scaling, using D instead of K)
% options: 'SchurD' (only for dynamic, like k-scaling, using the schur-complement of D instead of K)
p.Scaling = 'K';

% if true/1: do matrix computations substructure by substructure
% does not mean parallelization in time
% strongly recommended to be set true/1
p.Parallel = 1;

% select preconditioner
% options: 'K'       (lumped)
% options: 'SchurK'  (dirichlet)
% options: 'D'       (like lumped, using D instead of K)
% options: 'SchurD'  (like dirichlet, using the schur-complement of D instead)
p.Preconditioner = 'K';

% select the matrix Q
% options: 'eye'  (set Q to identity)
% options: 'pre'  (use the preconditioner as Q)
p.SwitchQ = 'pre'; % eye or pre

% choose tolerance (stopping criterion for CG iterations: res < tol)
p.tol = 1e-8;

% select the solver
% the so-named function will be used as solver
% that solver function shall be located in
% FETI/staticsolver or FETI/dynamicsolver
p.Solver = 'FETIstaticsolverExp2';

% do not call the solver
% (e.g. if you just want eigenspectrum plots of the operator)
p.NoSolve = 0;




%% post processing
% do postprocessing (plotting) after solving
p.PostProc = 1;

% give this experiment some name
% the directory in which output files are written will be called like this
p.Experiment = mfilename;

% write plots to svg file
p.WriteSVG = 0;

% do not write to any files at all
p.NoWrite = 1;

% plot correlations... ?
p.PlotCorrel = 0;

% plot eigenspectrum of operator
p.PlotEig = 0;

% plot residual
p.PlotRes = 1;

% plot deformed structure
p.PlotDeform = 1;

% select which deformations to plot (as array)
% options: 1 Actual State
% options: 2 Actual State Without RBM
% options: 3 Actual State Exploded
% options: 4 Undeformed State
% options: 5 Undeformed State Exploded
% options: 6 Actual State Without RBM Exploded
% example: [1 3 4] plots actual state, actual state exploded and undeformed
% state
p.StaticPlots = [1];

% always close the last deform plot (with rbm, without rbm etc.) before
% plotting the next e.g. if they should just be saved to files
p.PlotDeformClosePrev = 0;

% which iterations to plot
p.PlotIterations = [];

% whether to plot the last (final/converged) iteration
p.PlotLastIteration = 1;

% plot geneo modes
p.PlotGeneo = 0;

% select substructures for which geneo modes should be plotted
p.PlotGeneoSubs = 3;

% number of geneo modes to plot
p.PlotGeneoMax = 16;

% display eigenfrequencies and damping parameters
p.DisplayEigen = 0;
p.DisplayDamping = 0;




%% geometry of the structure
p.Nely = 4; % number of elements in y-direction
p.Nelx = 7; % number of elements in x-direction
p.elThick = 0.1; % thickness of structure in meters
p.elHeight = 1; % height and length of one FE element in meters
p.Nsy = 3; % number of substructures in y direction
p.Nsx = 5; % number of substructures in x direction

% every substructure is rectangular, build out of Nely times Nelx FE elements
% for every row of FE elements (counting starts from the bottom, going in y direction)
% you can choose the material:
% e.g.
% p.SteelRowNrsOddNsx = [2 4];
% means that the 2nd and fourth row will be steel, rest rubber
%
% there are two variables of this kind
% SteelRowNrsOddNsx sets the material for substructures which have an odd
% number counting in x direction
% SteelRowNrsEvenNsx sets the material for substructures which have an even
% number counting in x direction
%
p.SteelRowNrsOddNsx = [2 4];
p.SteelRowNrsEvenNsx = [2 4];

% if ChangeForEvenNsy is set true/1, the material setting for Odd/Even is
% switched the other way round for every row of substructures
p.ChangeForEvenNsy = 0;




%% loading
p.Loadcase = 1;
p.bendforce = 5000;
p.axforce = 10000;
p.clamping = 4;




%% other
p.Figh = struct( ...
    'structure', 1, ...
    'mac', 2, ...
    'eig', 3, ...
    'its', 4, ...
    'geneomodes', 5 ...
    );




%% set up cases
% cases have to be defined for every simulation run
% this means the variable Params is a structure array, in which every
% element contains a different set of parameters p
%
% usually you want to define some general settings above, and then create
% the different cases by just copying the general case p like
% Params(<CaseNr>).p = p;
% and then adapting the values for this case e.g. changing the changing the
% coarse grid setting
% Params(<CaseNr>).p.Coarse = 2;
%
% create Case Nr. 1:
% Case = 1;
% Params(Case).p = p;
% % change tolerance for Case Nr. 1:
% Params(Case).p.tol = 1e-8;
% % you HAVE to set the Case and some description:
% Params(Case).p.Case = Case;
% Params(Case).p.description = 'some case';
% % the variable CaseRuns sets how often that same case is run e.g. to get
% % some average value
% Params(Case).p.CaseRuns = 1;
% 
% Case = Case + 1;
% Params(Case).p = p;
% % change tolerance for Case Nr. 1:
% Params(Case).p.tol = 1e-6;
% % you HAVE to set the Case and some description:
% Params(Case).p.Case = Case;
% Params(Case).p.description = 'some other case';
% % the variable CaseRuns sets how often that same case is run e.g. to get
% % some average value
% Params(Case).p.CaseRuns = 1;

CaseNr = 0;

CaseNr = CaseNr + 1;
% case 1:
Params(CaseNr).p = p;
Params(CaseNr).p.Case = CaseNr;
Params(CaseNr).p.Odd = 1;
Params(CaseNr).p.description = '';
Params(CaseNr).p.Nely = 4; % number of elements in y-direction of a substructure
Params(CaseNr).p.Nelx = 7; % number of elements in x-direction of a substructure
Params(CaseNr).p.Nsy = 2; % number of substructures in y-direction
Params(CaseNr).p.Nsx = 4; % number of substructures in x-direction
Params(CaseNr).p.CaseRuns = 1;
Params(CaseNr).p.NoPattern = 1;



%% execute calculations
% display information for all cases to be run
for Case = 1:length(Params)
    display(['Case length=' num2str(Params(1).p.elHeight*Params(1).p.Nelx*Params(1).p.Nsx) ' m;   height=' num2str(Params(1).p.elHeight*Params(1).p.Nely*Params(1).p.Nsy) ' m']);
end




% start calculation
%parpool(min(length(Params)));
for Case = 1:length(Params)
    FETI(Params(Case).p);
end

%close all;