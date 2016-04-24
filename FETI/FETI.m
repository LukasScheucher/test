function [p] = FETI(p)
display(char(13));
display('-------------------------------------------------------------------');
display(['Case ' num2str(p.Case) ': ' p.description]);

% length = Params(CaseNr).p.elHeight*Params(CaseNr).p.Nelx*Params(CaseNr).p.Nsx;
% height = Params(CaseNr).p.elHeight*Params(CaseNr).p.Nely*Params(CaseNr).p.Nsy;
% display(['length=' num2str(length) ' m;   height=' num2str(height) ' m']);
% toller neuer kommentar

[p] = FETIsetup(p);

%% setup B matrices (interface connections)
[p] = FETIgeometry(p);

disp('Building geometry finished')

[p] = FETIload(p);

[p] = FETIassembly(p);

disp('Assembly finished')

if strcmp(p.Solver,'FULLsolver')~=1
    [p] = FETIscaling(p);
    [p] = FETIpreconditioner(p);
    [p] = FETIcoarse(p);
end


%% pre-solve plots
if p.PlotCorrel
    FETIplotMACcoarse(p)
end
if p.PlotEig
    FETIplotOperatorEig(p);
end
if isfield(p,'CondiNr') && p.CondiNr
if strcmp(p.mode,'dynamic')
    display(['Condition number: ', num2str(cond(p.FItransPre*p.Ptrans'*p.FItrans*p.Ptrans))]);
elseif strcmp(p.mode,'static')
    display(['Condition number: ', num2str(cond(p.FIPre*p.P'*p.FI*p.P))]);
end
end
    
%% solve
if (p.NoSolve)
    return;
end
disp('Solver started')
for i = 1:p.CaseRuns
    p.Subcase = i;
    p = eval([p.Solver '(p)']);
end
% p.Odd_cS = zeros(10);
% for i= 1:10
%     for j = 1:10
%         p.Odd = i;
%         display(['p.Odd: ' num2str(p.Odd)]);
%         p.cS = j;
%         display(['p.cS: ' num2str(p.cS)]);
%         p = eval([p.Solver '(p)']);
%     end
% end

%% post-solve plots
if p.PlotRes %strcmp(p.Solver,'FULLsolver')==0
    FETIplotRes(p);
end

if p.PostProc
    [p]=FETIpost(p);
end

if isfield(p,'PlotOperatorEig') && p.PlotOperatorEig
    FETIplotOperatorEig(p);
end


end