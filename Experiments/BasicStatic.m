clear variables;
close all;

p.mode = 'static';

%% feti parameters
% steel
p.ESt = 200e6;
p.nuSt = 0.3;
p.rhoSt = 7850.0;

% rubber parameters
p.ERub = 0.2e6;
p.nuRub = 0.5;
p.rhoRub = 920.0;

p.SteelRowNrsOddNsx = [1 2 3 4]; % the SteelrowNr'th element row will be steel, rest rubber
p.SteelRowNrsEvenNsx = [1 2 3 4];
p.ChangeForEvenNsy = 0;


%% feties
p.Coarse = 2;  % 0: no coarse grid 1: use RBM; 2: use GenEO, see mode selection
p.GeneoModes = 4;  % number of Geneo Modes
p.GeneoPlotSubs = 3;
p.Scaling = 'K';  % eye, multi, K, SchurK, D, SchurD
p.PostProc = 1;
p.Parallel = 1;  % do matrix computations substructure by substructure
p.Integrator = 3;  % 1: manual; 2: newmark; 3: gen-alpha
p.Preconditioner = 'K'; % K, SchurK, D, SchurD
p.Solver = 'feti'; % feti or Full
p.SwitchQ = 'eye';

p.WriteVideo = 0;
p.WriteSVG = 0;
p.DeformPlot = 0;
p.CorrelPlot = 0;
p.EigPlot = 0;
p.NoWrite = 1;

p.Loadcase = 1;
p.timestepNo = 1;

p.Experiment = mfilename;
p.EigOnly = 0;
p.PlotIterations = [];
p.PlotLastIteration = 1;
p.PlotGeneo = 0;


p.t = 1;


%% feti
p.Nely = 4; % number of elements in y-direction
p.Nelx = 4; % number of elements in x-direction
p.elThick = 0.002; % in meters
subHeight = 0.090;
p.elHeight = subHeight/p.Nely;
p.Nsy = 1;
p.Nsx = 3;


%% loading
p.Loadcase = 1;
p.bendforce = -100;


%% other
p.Figh = struct( ...
    'structure', 1, ...
    'mac', 2, ...
    'eig', 3, ...
    'its', 4, ...
    'geneomodes', 5 ...
    );

Params(1).p = p;




%% execute calculations
for CaseNr = 1:length(Params)
    display(['Case length=' num2str(Params(1).p.elHeight*Params(1).p.Nelx*Params(1).p.Nsx) ' m;   height=' num2str(Params(1).p.elHeight*Params(1).p.Nely*Params(1).p.Nsy) ' m']);
end

%parpool(min(length(Params)));
for CaseNr = 1:length(Params)
    FETIstatic(Params(CaseNr).p);
end

%close all;