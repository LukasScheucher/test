function [p] = FETIcoarse(p)

if strcmp(p.mode,'dynamic')
    % compute GenEO "bad" modes
    for s = 1:p.Ns
        [p.VqGen{s}, DlambdaGen{s}] = eig(p.SDs{s},p.Bbs{s}'*p.FItransPre*p.Bbs{s});
        % sort eigenvalues and corresp. eig.vec. in ascending order
        [p.VqGen{s}, DlambdaGen{s}] = sortem(p.VqGen{s},DlambdaGen{s});
         p.GeneoEigenValues{s}= diag(DlambdaGen{s});
        % use small eigvalues until jump in eigvalue
    end
    
    if isfield(p,'saveGeneoEigenValues') && p.saveGeneoEigenValues
    for s = 1:p.Ns
        if p.Ns == 4
            eigs4{s} = p.GeneoEigenValues{s};
        end
        if p.Ns == 16
            eigs16{s} = p.GeneoEigenValues{s};
        end
        if p.Ns == 9
            eigs9{s} = p.GeneoEigenValues{s};
        end
        if p.Ns == 36
            eigs36{s} = p.GeneoEigenValues{s};
        end
        if p.Ns == 144
            eigs144{s} = p.GeneoEigenValues{s};
        end
    end
    
    if p.Ns == 4
            save('daten4.mat','eigs4');
        end
        if p.Ns == 16
            save('daten16.mat','eigs16');
        end
        if p.Ns == 9
            save('daten9.mat','eigs9');
        end
        if p.Ns == 36
            save('daten36.mat','eigs36');
        end
        if p.Ns == 144
            save('daten144.mat','eigs144');
        end
    end
    % 1: use Rigid Body Modes
    if p.Coarse == 1
        p.GI = p.Bb*p.R;
    end
    
    % 2: use Geneo Modes
    % 3: use Geneo Modes + Rigid Body Modes
    if p.Coarse == 2 || p.Coarse == 3
        i = 0;
        for s = 1:p.Ns
            if iscell(p.GeneoModes)
                GeneoModeSelection = p.GeneoModes{s};
            else
                GeneoModeSelection = p.GeneoModes;
            end
            for GeneoMode = GeneoModeSelection
                i = i + 1;
                p.GI(:,i) = p.FItransPre*p.Bbs{s}*p.VqGen{s}(:,GeneoMode);
                %p.GI(:,i) = p.Bbs{s}*p.VqGen{s}(:,GeneoMode);
            end
        end
    end
    
    % 3: use Geneo Modes + Rigid Body Modes
    if p.Coarse == 3
        p.GI = [p.GI p.Bb*p.R];
    end

    if p.Coarse
        p.CoarseOp = p.GI'*p.FItrans*p.GI;
        p.Ptrans = eye(p.Nlm) - p.GI*(p.CoarseOp\(p.GI'*p.FItrans));
    end
    if isfield(p,'sizeGI') && p.sizeGI
        display(['Coarse-Grid size: ', num2str(size(p.GI,2))]);
    end
elseif strcmp(p.mode,'static')
    % compute GenEO "bad" modes
    for s = 1:p.Ns
        [p.VqGen{s}, DlambdaGen{s}] = eig(p.SKs{s},p.Bbs{s}'*p.FIPre*p.Bbs{s});
        % sort eigenvalues and corresp. eig.vec. in ascending order
        [p.VqGen{s}, DlambdaGen{s}] = sortem(p.VqGen{s},DlambdaGen{s});
        
        % delete rigid body modes (zero eigenvalues)
        p.VqGen{s}(:,1:p.Nrbm_s(s)) = [];
        DlambdaGen{s}(:,1:p.Nrbm_s(s)) = [];
        DlambdaGen{s}(1:p.Nrbm_s(s),:) = [];
        
        p.GeneoEigenValues{s} = diag(DlambdaGen{s});
        % use small eigvalues until jump in eigenvalue
    end
if isfield(p,'saveGeneoEigenValues') && p.saveGeneoEigenValues    
    for s = 1:p.Ns
        if p.Ns == 4
            Seigs4{s} = p.GeneoEigenValues{s};
        end
        if p.Ns == 16
            Seigs16{s} = p.GeneoEigenValues{s};
        end
        if p.Ns == 9
            Seigs9{s} = p.GeneoEigenValues{s};
        end
        if p.Ns == 36
            Seigs36{s} = p.GeneoEigenValues{s};
        end
    end
    
    if p.Ns == 4
            save('datenS4.mat','Seigs4');
        end
        if p.Ns == 16
            save('datenS16.mat','Seigs16');
        end
        if p.Ns == 9
            save('datenS9.mat','Seigs9');
        end
        if p.Ns == 36
            save('datenS36.mat','Seigs36');
        end
end
        
    switch p.SwitchQ
        case 'eye'
            p.Q = eye(p.Nlm);
        case 'pre'
            p.Q = p.FIPre;
    end

    % anyway: Rigid Body Modes as coarse grid
    p.GI = p.Bb*p.R;
    p.NCoarseOp = p.GI'*p.Q*p.GI;
    p.PN = eye(p.Nlm) - p.Q*p.GI*(p.NCoarseOp\(p.GI'));
    
    % 2: use Geneo Modes
    if p.Coarse == 2
        i = 0;
        for s = 1:p.Ns
            if iscell(p.GeneoModes)
                GeneoModeSelection = p.GeneoModes{s};
            else
                GeneoModeSelection = p.GeneoModes;
            end
            for GeneoMode = GeneoModeSelection
                i = i + 1;
                %p.GI0(:,i) = p.PN*p.FIPre*p.Bbs{s}*p.VqGen{s}(:,GeneoMode);
                p.GI0(:,i) = p.PN*p.FIPre*p.Bbs{s}*p.VqGen{s}(:,GeneoMode);
                % added PN here --> improved rcond(N) from 10e-18 to
                % 10e-14
                % p.PN = 10 iterations, p.PN' = 26 iterations
            end
        end
        size(p.GI0)
        
        p.Ngm = i;
        
        p.FI0 = p.GI0'*(p.PN'*p.FI*p.PN)*p.GI0;
        p.P0 = eye(p.Nlm) - p.GI0*pinv(p.FI0)*p.GI0'*(p.PN'*p.FI*p.PN);
        p.P = p.PN*p.P0;
        
        p.PC = eye(p.Nlm) - p.GI0*pinv(p.GI0'*p.FI*p.GI0)*p.GI0'*p.FI;
        
        temp1 = p.GI'*p.GI0;
        temp2 = p.GI0'*p.GI;
        p.N = [p.GI0'*p.FI*p.GI0 p.GI0'*p.FI*p.Q*p.GI temp2; ...
            p.GI'*p.Q*p.FI*p.GI0 p.GI'*p.Q*p.FI*p.Q*p.GI p.GI'*p.Q*p.GI; ...
            temp1 p.GI'*p.Q*p.GI zeros(size(temp1,1),size(temp2,2))];
    end
    
    % use use the condition of exact matching corners as coarse grid
    % + use Geneo Modes
    % no RBMs!
    if p.Coarse == 3
        %CornerDOFs = sum(abs(p.B),1)==3;
        %CornerLMs = find(sum(abs(p.B(:,CornerDOFs)),2) ~= 0);
        p.PN = eye(p.Nlm);
        
        CornerLMsNrs = find(p.CornerLMs);
        
        col = 0;
        p.GI0 = zeros(p.Nlm,(length(CornerLMsNrs)+length(p.GeneoModes)));
        for CornerLM = CornerLMsNrs'
            col = col + 1;
            p.GI0(CornerLM,col) = 1;
            p.GI0(:,col) = p.PN*p.GI0(:,col);
        end
        
        for s = 1:p.Ns
            for GeneoMode = p.GeneoModes
                col = col + 1;
                p.GI0(:,col) = p.PN*p.FIPre*p.Bbs{s}*p.VqGen{s}(:,GeneoMode);
            end
        end
        
        p.PC = eye(p.Nlm) - p.GI0*pinv(p.GI0'*p.FI*p.GI0)*p.GI0'*p.FI;
    end
    
    % use use the condition of exact matching corners as coarse grid
    if p.Coarse == 4
        %CornerDOFs = sum(abs(p.B),1)==3;
        %CornerLMs = find(sum(abs(p.B(:,CornerDOFs)),2) ~= 0);
        p.PN = eye(p.Nlm);
        
        CornerLMsNrs = find(p.CornerLMs);
        
        col = 0;
        p.GI0 = zeros(p.Nlm,(length(CornerLMsNrs)+length(p.GeneoModes)));
        for CornerLM = CornerLMsNrs'
            col = col + 1;
            p.GI0(CornerLM,col) = 1;
            p.GI0(:,col) = p.PN*p.GI0(:,col);
        end
        
        p.PC = eye(p.Nlm) - p.GI0*pinv(p.GI0'*p.FI*p.GI0)*p.GI0'*p.FI;
    end
    
%     if .... ?
%         p.Ngm = 0;
%         p.GI0(1:size(p.Q,1),1) = 0;
%         
%         p.P = p.PN;
%         temp = p.GI'*p.Q*p.GI;
%         p.N = [p.GI'*p.Q*p.FI*p.Q*p.GI p.GI'*p.Q*p.GI; ...
%             p.GI'*p.Q*p.GI zeros(size(temp,1),size(temp,2))];
%     end
end

end