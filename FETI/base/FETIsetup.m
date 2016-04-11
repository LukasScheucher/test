function [p] = FETIsetup(p)
%% some basic setup

Checklist = FETIchecklist;
for i=1:size(Checklist,1)
    if isfield(p,Checklist{i,1})==0
        if isletter(Checklist{i,2})==0
            p.(Checklist{i,1})=str2double(Checklist{i,2});
        else
            p.(Checklist{i,1})=(Checklist{i,2});
        end
        disp(['Due to missing variable in Exmaple-file p.' Checklist{i,1} ' has been set to standard value: '])
        disp(p.(Checklist{i,1}))
    end
end

if p.cal_stress==1
    p.cal_strains=1;
end

if strcmp(p.mesh_method,'Srd-LM')
    p.nonconforming=0;
else
    p.nonconforming=1;
end

if strcmp(p.Solver,'FULLsolver') || strcmp(p.mode,'dynamic')
    p.globalassembly=1;
else
    p.globalassembly=0;
end

if strcmp(p.mode,'static')
    p.t = 1;
end

if exist(p.Experiment,'dir') ~= 7
    mkdir(p.Experiment);
end


if p.Coarse == 0
    p.CoarseGrid = 'NoCoarse';
elseif p.Coarse == 1
    p.CoarseGrid = 'RBM';
elseif p.Coarse == 2
    p.CoarseGrid = 'GenEO';
end

end