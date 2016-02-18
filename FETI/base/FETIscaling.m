function [p] = FETIscaling(p)

switch p.Scaling
    case 'eye'
        p.WBi = eye(p.Nlm)*p.Bb; % p.WBi: p.WB matrix for identity scaling
        p.WB = p.WBi;
        for s = 1:p.Ns
            p.WBs{s} = p.Bbs{s};
        end
    case 'multi'
        %p.WBm = p.Winv*p.Bb; % p.WBm: p.WB matrix for multiplicity scaling
        BABBA = pinv(p.Bb*p.Bb')*p.Bb; % BABBA should be p.WBm
        p.WB = BABBA;
        for s = 1:p.Ns
            p.WBs{s} = pinv(p.Bb*p.Bb')*p.Bbs{s};
        end
    case 'K'
        %p.WBk = betaScale*p.Bdiag; % p.WBk: p.WB matrix for k-scaling
        A = inv(diag(diag(p.Kbb)));
        BABBA = pinv(p.Bb*A*p.Bb')*p.Bb*A; % BABBA should be p.WBk
        p.WB = BABBA;
        for s = 1:p.Ns
            A = inv(diag(diag(p.Kbb)));
            As = inv(diag(diag(p.Ksbb{s})));
            p.WBs{s} = pinv(p.Bb*A*p.Bb')*p.Bbs{s}*As;
        end
    case 'SchurK'
        A = inv(diag(diag(p.SK)));
        BABBA = pinv(p.Bb*A*p.Bb')*p.Bb*A;
        p.WB = BABBA;
    case 'D'
        if strcmp(p.mode,'dynamic')
            A = inv(diag(diag(p.Dbb))); % D-scaling
            BABBA = pinv(p.Bb*A*p.Bb')*p.Bb*A; % BABBA should be p.WBk
            p.WB = BABBA;
            
            temp = pinv(p.Bb*A*p.Bb');
            for s = 1:p.Ns
                As = inv(diag(diag(p.Dsbb{s})));
                p.WBs{s} = temp*p.Bbs{s}*As;
            end
        else
            display('D-Scaling only available in dynamic mode!');
            return;
        end
    case 'SchurD'
        if strcmp(p.mode,'dynamic')
            A = inv(diag(diag(p.SD))); % SD-scaling
            BABBA = pinv(p.Bb*A*p.Bb')*p.Bb*A;
            p.WB = BABBA;
        else
            display('Schur-D-Scaling only available in dynamic mode!');
            return;
        end
end

end