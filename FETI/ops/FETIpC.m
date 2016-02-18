function y = FETIpC(x,transpose,p)
    if p.Coarse < 2
        y = x;
        return;
    end
    
    if ~isfield(p,'CCoarseOp')
            p.CCoarseOp = p.GI0'*p.FI*p.GI0;
    end
    
    if transpose
        y = x - p.FI*(p.GI0*(p.CCoarseOp\(p.GI0'*x)));
    else
        y = x - p.GI0*(p.CCoarseOp\(p.GI0'*(p.FI*x)));
    end
end