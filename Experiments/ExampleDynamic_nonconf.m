clear variables;
% close all;

p.mode = 'dynamic';

%% feti parameters
% steel
p.ESt = 200e6; % in kN/m^2
p.nuSt = 0.3;
p.rhoSt = 7.85; % in g/cm^3

% rubber
p.ERub = 0.2e6; % in kN/m^2
p.nuRub = 0.5;
p.rhoRub = 0.92; % in g/cm^3



%% feti
p.Coarse = 0;  % 0: no coarse grid 1: use RBM; 2: use GenEO, see mode selection
p.GeneoModes = [1 2 3];  % number of Geneo Modes % [1 2 3 4 5 6]
p.Scaling = 'D';  % eye, multi, K, SchurK, D, SchurD
p.PostProc = 1;
p.Parallel = 1;  % do matrix computations substructure by substructure
p.Integrator = 3;  % 1: manual; 2: newmark; 3: gen-alpha
p.Preconditioner = 'D'; % K, SchurK, D, SchurD
%p.Solver = 'FETIdynamicsolver'; % feti or Full _Te_Het
p.Solver = 'FULLsolver'; % direct solver

p.WriteVideo = 0;
p.WriteSVG = 0;
p.PlotDeform = 1;
p.PlotCorrel = 0;
p.PlotEig = 0;
p.NoWrite = 1;
p.PlotRes = 0; 
p.DisplayTimeSteps = 0;

p.DisplayEigen = 0;
p.DisplayDamping = 0;

% plot Iterations
p.opt.plotMAC = 0;


p.tol = 1e-5;

p.timestepNo = 1;

p.Experiment = mfilename;
p.NoSolve = 0;





p.Step = 2e-4; 
p.End = 500*p.Step; %100*p.Step;
p.t = 0:p.Step:p.End;
p.Plotstep = 5*p.Step;% 2*p.Step;


%% feti
% note,whether the mesh is nonconforming (0=no, 1=yes)
p.mesh_method='Mortar'; % Implemented Methods:
                        % - Srd-LM (Conforming meshes)
                        % - NTS-LM
                        % - Mortar

p.max_iteration=5; % stop solver, if iteration counter = p.max_iteration*Nlm

p.trackpoint=[1]; % 

%% geometry of the structure
% Parameters for nonconforming meshes (Note: Choose a suitable method with p.mesh_method)
%p.Height = 2; % cantilever height in meters
%p.Length = 2; % cantilever length in meters
p.globalassembly=1;
%p.sizes = [1, 1, 1, 1;
%            1, 1, 1, 1];
%p.elcount = [1, 2, 2, 3];
%p.sizes = [3,1,2;
%            1,0.5,0.5];
%p.elcount = [2,4,2];
%p.sizes = [6/5, 1, 4/5, 3;
%               1,1,1,0.5]; % size of subsstructures in meters from first substructure to last; first line = length, second line = height
%p.elcount = [5, 2, 5, 2]; % element count in y-direction for each substructure from first to last
p.Height = 1.5; % cantilever height in meters
p.Length = 3; % cantilever length in meters
%p.sizes = [4/3, 1, 2/3, 2, 1;
%            1, 1, 1, 0.5, 0.5];
%p.elcount = [3, 1, 3, 1, 2];
p.sizes = [3, 3; 
    0.75, 0.75];
p.elcount = [2, 2];
%p.sizes = [1, 2, 2;
%            1, 0.5, 0.5];
%p.elcount = [3, 1, 2];

p.geom_tol = 1e-9; % Global tolerance for control of floating-point operations (e.g. geometric positioning vectors)

p.elThick = 0.1;
p.StaticIterations = 1; % Do not solve
p.Nsy = 3;
p.Nsx = 3;
% Parameters for conforming meshes (Srd-LM)
p.Nely = 1; % number of elements in y-direction %4
p.Nelx = 2; % number of elements in x-direction %7
p.Nsy = 1; % number of substructures in y-direction %2 
p.Nsx = 2; % number of substructures in x-direction %5
p.elThick = 0.1; % in meters
p.elHeight = 0.5; % in meters#


p.SteelRowNrsOddNsx = []; %[2 4] % the SteelrowNr'th element row will be steel, rest rubber
p.SteelRowNrsEvenNsx = []; %[2 4]
p.ChangeForEvenNsy = 0;

%% NEW
% Irreguläre Heterogenität
p.NoPattern = 0; % an/aus -> falls die Option an ist wird das Problem nicht nach dem 
                 % odd-even Schmea aufgebaut. Die Substrukturen 
                 % in p. sind dann die Einzigen mit dem  Muster aus p.SteelRowNrsOddNsx
                 % alle restlichen Substrukturen haben p.SteelRowNrsEvenNsx
                 % als Vorlage -> gedacht um einzelne "schwierige" Substrukturen in
                 % einem ansonsten homogenen Problem zu erzeugen
p.Odd = []; % 5,6
% Einspannung ändern
p.clamping = 4; % 1-oben, 2-rechts, 3-unten, 4-links

% Belastung ändern
p.Loadcase = [1]; % 1-rechts_oben y, 2-rechts_unten y, 3-Links_unten y, 4-links_oben y,
                  % 5-rechts_oben x, 6-rechts_unten x, 7-Links_unten x, 8-links_oben x,
                  % 9-mittig_oben1 y, 10-mittig_oben2 y
                  % Mehrfachwahl möglich
                
                 
% plot GenEo Modes
p.PlotGeneo = 0;
p.PlotGeneoSubs = [p.Odd];
p.PlotGeneoMax = 12; %size(p.GeneoModes,2)* p.Nsy * p.Nsx;

%% loading

BendforceAmp = -2.0; % in kN (10 ~ 1 ton) %-2.0;
BendforceShock = -1/p.Step;
BendforceStatic = -100.0; % in kN (10 ~ 1 ton)
%p.Loadcase = 6;
period = 2.0;
if 0
    p.Loading = sin(2*pi./period*p.t).*BendforceAmp;
elseif 0
    p.Loading(1) = BendforceShock;
    p.Loading(2:length(p.t)) = 0.0;
elseif 1
    frac = 8;
    p.Loading(1:ceil(length(p.t)/frac)) = p.t(1:ceil(length(p.t)/frac))./p.t(ceil(length(p.t)/frac))*BendforceStatic;
    p.Loading(ceil(length(p.t)/frac):ceil(length(p.t)*2/frac)) = BendforceStatic;
    p.Loading(ceil(length(p.t)*2/frac):length(p.t)) = 0;
elseif 0
    p.Loading(1:ceil(length(p.t)/2)) = p.t(1:ceil(length(p.t)/2))./p.t(ceil(length(p.t)/2))*BendforceStatic;
    p.Loading(ceil(length(p.t)/2):ceil(length(p.t))) = BendforceStatic;
elseif 0
    frac = 8;
    p.Loading(1:ceil(length(p.t)/frac)) = BendforceStatic;
    p.Loading(ceil(length(p.t)/frac):ceil(length(p.t)*2/frac)) = BendforceStatic;
    p.Loading(ceil(length(p.t)*2/frac):length(p.t)) = 0;
elseif 0
    t_loadend = 25;
    p.Loading(1:t_loadend) = BendforceStatic;
    p.Loading(t_loadend:length(p.t)) = 0;
end



%% numerical parameters
temp = p.Integrator;
switch temp
    case 1
        p.GAalpham = 0.3; % GAalpham < GAalphaf < 0.5
        p.GAalphaf = 0.4;
        p.GAgamma = 0.5 + p.GAalphaf - p.GAalpham;
        p.GAbeta = 0.25 + 0.5*(p.GAalphaf - p.GAalpham);
    case 2
        p = struct( ...
            'GAalphaf', 0, ...
            'GAalpham', 0, ...
            'GAgamma', 0.5, ...
            'GAbeta', 0.25 ... % GAbeta > 0.25
            );
    case 3
        % after Chung and Hulbert optimal for 2nd order ODEs:
        rhoInf = 0.9;
        p.GAalphaf = rhoInf/(rhoInf+1); % GAalpham < GAalphaf < 0.5
        p.GAalpham = (2*rhoInf - 1)/(rhoInf + 1);
        p.GAgamma = 0.5 + p.GAalphaf - p.GAalpham;
        p.GAbeta = 0.25*(p.GAgamma + 0.5)^2;
end

%% other
p.Figh = struct( ...
    'structure', 1, ...
    'mac', 2, ...
    'eig', 3, ...
    'its', 4, ...
    'geneomodes', 5 ...
    );




CaseNr = 0;
CaseNr = CaseNr + 1;

% case 1: 
p.description = 'stehend_TeX'; %'heterogeneous stripes, classical solver with geneo [1 2 3] grid';

p.Case = CaseNr;
p.Subcase = 1;
p.CaseRuns = 1; %1;

if strcmp(p.mesh_method,'Srd-LM')
    length = p.elHeight*p.Nelx*p.Nsx;
    height = p.elHeight*p.Nely*p.Nsy;
else
    length = p.Length;
    height = p.Height;
end
FETI(p);



