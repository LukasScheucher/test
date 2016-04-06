function Checklist=FETIchecklist
% Checklist whether all new parameters exist
% If not, set a defined standard value for conforming meshes

Checklist={'mesh_method','Srd-LM';
    'nonconforming','0';
    'max_iteration','1';
    'geom_tol','1e-9';
    'globalassembly','0';
    'cal_strains','0';
    };
disp(size(Checklist))

end