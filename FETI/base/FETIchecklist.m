function Checklist=FETIchecklist
% Checklist whether all new parameters exist
% If not, set a defined standard value for conforming meshes

Checklist={'mesh_method','Srd-LM';
    'nonconforming','0';
    'max_iteration','1';
    'geom_tol','1e-9';
    'globalassembly','0';
    'cal_strains','0';
    'cal_stress','0';
    'strain_dir','1';
    'addNTSLMs','0';
    'plot_int','0';
    'plain','1';
    'novertDBC','0';
    'axconststress','0';
    'plot','';
    'NoArrows','0';
    'forceoffset','0';
    };
disp(size(Checklist))

end