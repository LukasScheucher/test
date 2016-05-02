function y = FETIpN(x,transpose,p)
    if p.Coarse == 3 || p.Coarse == 4 || isempty(p.GI)
        y = x;
        return;
    end
    
    if ~isfield(p,'NCoarseOp')
        p.NCoarseOp = p.GI'*p.Q*p.GI;
    end
    
    if transpose
        y = x - p.GI*(p.NCoarseOp\(p.GI'*(p.Q*x)));
    else
        y = x - p.Q*(p.GI*(p.NCoarseOp\(p.GI'*x)));
    end
end