function [p] = FETIpreconditioner(p)
%% select preconditioner

switch p.Preconditioner
    case 'K'
        
        p.FIPreUnscaled = zeros(p.Nlm);
        p.FIPre = zeros(p.Nlm);
        for s = 1:p.Ns
            p.FIPreUnscaled = p.FIPreUnscaled + p.Bbs{s}*p.Ksbb{s}*p.Bbs{s}';
            p.FIPreS{s} = p.WBs{s}*p.Ksbb{s}*p.WBs{s}';
            p.FIPre = p.FIPre + p.FIPreS{s};
        end

        %p.FIPreUnscaledz = p.Bb*p.Kbb*p.Bb';
        %p.FIPrez = p.WB*p.Kbb*p.WB';
        
        %norm(p.FIPre - p.FIPrez)
        
    case 'SchurK'
        p.FIPreUnscaled = p.Bb*p.SK*p.Bb';
        p.FIPre = p.WB*p.SK*p.WB';
    case 'D'
        
        p.FItransPre = zeros(p.Nlm);
        for s = 1:p.Ns
            p.FItransPreS{s} = p.WBs{s}*p.Dsbb{s}*p.WBs{s}';
            p.FItransPre = p.FItransPre + p.FItransPreS{s};
        end
        
        p.FItransPrez = p.WB*p.Dbb*p.WB';
        %norm(p.FItransPre - p.FItransPrez)
    case 'SchurD'
        p.FItransPre = p.WB*p.SD*p.WB';
end
end