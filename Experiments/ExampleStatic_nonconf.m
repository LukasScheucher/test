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
%p.Solver = 'FULLsolver';
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
p.StaticPlots = 1;

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

% note,whether the mesh is nonconforming (0=no, 1=yes)
p.mesh_method='NTS-LM'; % Implemented Methods:
                        % - Srd-LM (Conforming meshes)
                        % - NTS-LM
                        % - Mortar
                        
p.plot='disp';      % Choose the plot, you want to see: '' - default, 'disp' - displacements, 'strain' - strains, 'stress' - stresses
p.strain_dir=1;     % Direction of strain: 1: epsilon_11, 2: epsilon_22, 3: epsilon_12, 4: von Mises stress (use only for stresses)
p.addNTSLMs=0;      % Add additional LMs for NTS (experimental)
p.plot_int=0;       % Interface to plot stresses/strains
p.novertDBC=0;      % no vertical DBC (!Choose only for Clamping=4 and load in x-direction!)
p.axconststress=0;  % For patch test: half the axial load at outer nodes due to assembly to get constant displacement field
p.plain=1;          % 1: plain stress, 2: plain strain

p.max_iteration=5; % stop solver, if iteration counter = p.max_iteration*Nlm
%% geometry of the structure
% Parameters for nonconforming meshes (Note: Choose a suitable method with p.mesh_method)
%p.Height = 2; % cantilever height in meters
%p.Length = 2; % cantilever length in meters
%p.sizes = [1, 1, 1, 1;     % B-Matrix Notizenbeispiel_NTS2x2.xlsx
%            1, 1, 1, 1];
%p.elcount = [1, 2, 2, 3];
%p.sizes = [3,1,2;
%            1,0.5,0.5];
%p.elcount = [2,4,2];
p.sizes = [6/5, 1, 4/5, 3;
               1,1,1,0.5]; % size of subsstructures in meters from first substructure to last; first line = length, second line = height
p.elcount = [5, 2, 5, 2]; % element count in y-direction for each substructure from first to last
p.Height = 1.5; % cantilever height in meters
p.Length = 3; % cantilever length in meters
%p.sizes = [4/3, 1, 2/3, 2, 1; % B-Matrix_assembliert.xlsx
%            1, 1, 1, 0.5, 0.5];
%p.elcount = [3, 1, 3, 1, 2];
%p.sizes = [1, 2, 2;
%            1, 0.5, 0.5];
%p.elcount = [3, 1, 2];

p.elThick = 0.1;
p.StaticIterations = 1; % Do not solve
p.geom_tol = 1e-9; % Global tolerance for control of floating-point operations (e.g. geometric positioning vectors)
    
% Parameters for conforming meshes (Srd-LM)
p.Nely = 2; % number of elements in y-direction
p.Nelx = 6; % number of elements in x-direction
p.elThick = 0.1; % thickness of structure in meters
p.elHeight = 0.5; % height and length of one FE element in meters
p.Nsy = 2; % number of substructures in y direction
p.Nsx = 1; % number of substructures in x direction

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

p.SteelRowNrsOddNsx = [];
p.SteelRowNrsEvenNsx = [];

% if ChangeForEvenNsy is set true/1, the material setting for Odd/Even is
% switched the other way round for every row of substructures
p.ChangeForEvenNsy = 0;




%% loading
p.Loadcase = 1;
%p.Loadcase = 1;
p.bendforce = 100;
p.axforce = 500;
p.axforcefield_max = 1500;
p.axforcefield_offset = -2000;
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
p.Nely0=p.Nely;
p.Nelx0=p.Nelx;
p.elHeight0=p.elHeight;
for z=1:1
    CaseNr = CaseNr + 1;
    % case 1:
    Params(CaseNr).p = p;
    Params(CaseNr).p.Case = CaseNr;
    Params(CaseNr).p.Odd = 1;
    Params(CaseNr).p.description = '';
    Params(CaseNr).p.mesh_method=p.mesh_method;%1;
    if strcmp(Params(CaseNr).p.mesh_method,'Srd-LM')
        Params(CaseNr).p.Nely = double(p.Nely0 + int32((z - 1)*p.Nely0/2)); % number of elements in y-direction of a substructure
        Params(CaseNr).p.Nelx = double(int32(Params(CaseNr).p.Nely*p.Nelx0/p.Nely0)); % number of elements in x-direction of a substructure
        disp('p.Nelx: ')
        disp(Params(CaseNr).p.Nelx)
        disp('p.Nely: ')
        disp(Params(CaseNr).p.Nely)
        Params(CaseNr).p.elHeight = p.elHeight0*Params(CaseNr).p.Nely/p.Nely0; % number of elements in x-direction of a substructure
        Params(CaseNr).p.Nsy = p.Nsy; % number of substructures in y-direction
        Params(CaseNr).p.Nsx = p.Nsx; % number of substructures in x-direction
    else
        Params(CaseNr).p.Height = p.Height;
        Params(CaseNr).p.Length = p.Length;
        Params(CaseNr).p.sizes = p.sizes;
        Params(CaseNr).p.elcount = p.elcount;
    end
    Params(CaseNr).p.CaseRuns = z;
    Params(CaseNr).p.NoPattern = 0;
end


%% execute calculations
% display information for all cases to be run
for Case = 1:length(Params)
    if strcmp(Params(CaseNr).p.mesh_method,'Srd-LM')
        display(['Case length=' num2str(Params(Case).p.elHeight*Params(Case).p.Nelx*Params(Case).p.Nsx) ' m;   height=' num2str(Params(Case).p.elHeight*Params(Case).p.Nely*Params(Case).p.Nsy) ' m']);
    else
        display(['Case length=' num2str(Params(Case).p.Length) ' m;   height=' num2str(Params(Case).p.Height) ' m']);
    end
end



% start calculation
% parpool(min(length(Params)));
convergence=zeros(2,length(Params));
disp('Convergence')
disp(convergence)
n=1;
for Case = 1:length(Params)
    [p] = FETI(Params(Case).p);
    disp(['p.nonconforming: ' num2str(p.nonconforming)])
    disp(['Tracking: ' num2str(p.tracking,8)])
    disp(['globalassembly: ' num2str(p.globalassembly)])
    if p.globalassembly==1
        convergence(1,n)=size(p.L_man,2);
    else
        convergence(1,n)=size(p.B2,2)-rank(p.B2);
    end
    convergence(2,n)=p.tracking;
    disp('Convergence')
    disp(convergence)
    disp(['Nlm: ' num2str(p.Nlm)])
    n=n+1;
end

figure(10)
hold on
plot(convergence(1,1:end-1),convergence(2,1:end-1))
plot(convergence(1,end),convergence(2,end),'*r')
hold off
%close all;