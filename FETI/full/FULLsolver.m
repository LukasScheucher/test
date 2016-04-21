function p = FULLsolver(p)%Full, Integration, Geometry)

    % full version: solver
    if strcmp(p.mode,'dynamic')
        Dfull = p.Full.M + p.Step/2.*p.Full.C + p.Step^2/4.*p.Full.K;
    else
        Dfull = p.Full.K;
    end

    % initial state of full version
    p.uFull = zeros(size(p.L_man,2),length(p.t));%p.NdofFull,length(p.t)); %uFull = zeros(Geometry.NdofFull,length(Integration.t));
    uFullP = p.uFull;
    uFullPP = p.uFull;
    uFullAst = p.uFull;
    uFullAstP = p.uFull;

    % initial acceleration of full version
    TimerInitFull = tic;
    [LMfull, UMfull, pMfull] = lu(p.Full.M,'vector');
    rhsInitFull = p.Full.f(:,1) - p.Full.C*uFullP(:,1) - p.Full.K*p.uFull(:,1);
    uFullPP(:,1) = UMfull\(LMfull\(rhsInitFull(pMfull,:)));
    TimeInitFull = toc(TimerInitFull);
    display(['Time for calc. of inital acceleration (full model): ' num2str(TimeInitFull)]);

    % LU decomposition of iteration matrix of full version
    TimerLUfull = tic;
    [LDfull, UDfull, pDfull] = lu(Dfull,'vector');
    TimeLUfull = toc(TimerLUfull);
    display(['Time for calc. of LU composition of FETI.D (full model): ' num2str(TimeLUfull)]);

    if strcmp(p.mode,'dynamic')
        % start Integration.Stepping
        for n = 1:(length(p.t)-1)
            % predictors
            uFullAst(:,n+1) = p.uFull(:,n) + p.Step.*uFullP(:,n) + p.Step^2/4.*uFullPP(:,n);
            uFullAstP(:,n+1) = uFullP(:,n) + p.Step/2.*uFullPP(:,n);

            % right hand side
            rhsFull = p.Full.f(:,n+1) - p.Full.C*uFullAstP(:,n+1) - p.Full.K*uFullAst(:,n+1);
            
            % solve
            uFullPP(:,n+1) = UDfull\(LDfull\(rhsFull(pDfull,:)));

            % integration
            uFullP(:,n+1) = uFullP(:,n) + p.Step/2.*(uFullPP(:,n)+uFullPP(:,n+1));
            p.uFull(:,n+1) = p.uFull(:,n) + p.Step.*uFullP(:,n) + p.Step^2/4.*(uFullPP(:,n)+uFullPP(:,n+1));
            %display(['time: ' num2str(Integration.t(n+1))]);
        end
    else
            % right hand side
            rhsFull = p.Full.f;
            disp(size(p.Full.f))
            % solve
            p.uFull = UDfull\(LDfull\(rhsFull(pDfull,:)));
    end
    
    p.PlotRes = 0;
end
