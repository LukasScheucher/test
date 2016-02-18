function params = FETIdynamicsolverS(params)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));

% tune external load
if FETISTuneLoad
    AddForceMax = 0;
    AddForceMin = 1;
    AddForce = (AddForceMax*rand(Nlm,1) + AddForceMin*ones(Nlm,1));
    for n = 1:(length(t)-1)
        for s = 1:Ns
            fs{s}(:,n) = fs{s}(:,n) + Bs{s}'*AddForce;
        end
        f(:,n) = f(:,n) + B'*AddForce;
    end
end

% init values
u = zeros(Ndof,length(t));
uP = u;
uPP = u;
uAst = u;
% initial acceleration
Temp1 = (B*(M\B'));
Temp2 = (M\(f(:,1)-C*uP(:,1)-K*u(:,1)));
lambda = pinv(Temp1)*(B*Temp2);

if Parallel
    for s = 1:Ns
        Ndof_end = sum(Ndof_s(1:s));
        Ndof_start = Ndof_end - Ndof_s(s) + 1;
        uPPs{s} = Ms{s}\(fs{s}(:,1)-Cs{s}*uP(Ndof_start:Ndof_end,1)-Ks{s}*u(Ndof_start:Ndof_end,1)-Bs{s}'*lambda);
    end
    for s = 1:Ns
        Ndof_end = sum(Ndof_s(1:s));
        Ndof_start = Ndof_end - Ndof_s(s) + 1;
        uPP(Ndof_start:Ndof_end,1) = uPPs{s};
    end
else
    uPP(:,1) = M\(f(:,1)-C*uP(:,1)-K*u(:,1)-B'*lambda);
end

iterations = 0;
NumDirections = 0;
itdivisor = 0;
a(:,1) = uPP(:,1);

tic
for n = 1:(length(t)-1)
    
    % predictors
    aPre = (GAalphaf.*uPP(:,n)-GAalpham.*a(:,n))./(1-GAalpham);
    uAst(:,n+1) = u(:,n) + Step.*uP(:,n) + Step^2*(0.5-GAbeta).*a(:,n) + Step^2*GAbeta.*aPre;
    uAstP(:,n+1) = uP(:,n) + Step*(1-GAgamma).*a(:,n) + Step*GAgamma.*aPre;

    dtrans(1:Nlm,1) = 0;
    if Parallel
        for s = 1:Ns
            Ndof_end = sum(Ndof_s(1:s));
            Ndof_start = Ndof_end - Ndof_s(s) + 1;
            dtransPre{s} = Bs{s}*(Ds{s}\(fs{s}(:,n+1)-Cs{s}*uAstP(Ndof_start:Ndof_end,n+1)-Ks{s}*uAst(Ndof_start:Ndof_end,n+1)));
        end
        for s = 1:Ns
            dtrans = dtrans + dtransPre{s};
        end
    else
        dtrans = B*(D\(f(:,n+1)-C*uAstP(:,n+1)-K*uAst(:,n+1)));
    end

    % initialize
    clear lambda r Zk W2k Q2k Deltak fac
    
    k = 1;

    if Coarse
        lambda(1:Nlm,k) = GI*(CoarseOp\(GI'*dtrans));
    else
        Ptrans = eye(Nlm);
        lambda(:,k) = zeros(Nlm,1);
    end
    
    r(:,1) = Ptrans'*(dtrans - FItrans*lambda(:,k));
    for s = 1:Ns
        Zk{1}(:,s) = FItransPreS{s}*r(:,1);
        %Zk{1}(:,s) = Zk{1}(:,s)./norm(Zk{1}(:,s));
        if FETISeigDirections
            LocalNorm = norm(Bbs{s}*Ksbb{s}*Bbs{s}'*Zk{1}(:,s));
            GlobalNorm = norm(Bbs{s}*Bbs{s}'*FItransPre*Zk{1}(:,s));
            Amps(s) = GlobalNorm/LocalNorm;
        end
        BBT{s} = Bs{s}*Bs{s}';
    end
    
    if FETISeigDirections2
        SummedDirection = sum(Zk{1}(:,:),2);
        for s = 1:Ns
            GlobalNorm = norm(BBT{s}*SummedDirection);
            LocalNorm = norm(Zk{1}(:,s));
            Amps(s) = GlobalNorm/LocalNorm;
        end
    end
    
    if FETISeigDirections || FETISeigDirections2
        [~, DirectionIndex] = min(Amps);
        KeptDirections(1:s) = false;
        KeptDirections(DirectionIndex) = true;
    end
    
    if FETISminmaxDirections
        DirectionsAmps = sum(abs(Zk{1}(:,:)),1);
        [~, DirectionIndex1] = max(DirectionsAmps);
        [~, DirectionIndex2] = min(DirectionsAmps);
        KeptDirections(1:s) = false;
        KeptDirections(DirectionIndex1) = true;
        KeptDirections(DirectionIndex2) = true;
    end
    
    if FETISAverageZeroDirections
        % DirectionsAmps = sum(abs(Zk{1}(:,:)),1);
        DirectionsAmps = sqrt(sum((Zk{1}).^2,1));
        KeptDirections = DirectionsAmps > FETISDirectionTol;
        if isempty(find(KeptDirections,1))
            [maxDirectionAmp, DirectionIndex] = max(DirectionsAmps);
            KeptDirections(DirectionIndex) = true;
            %display(['Warning: No directions outside zero tolerance. Using direction with Amp = ' num2str(maxDirectionAmp)]);
        end
    end
    
    if FETISAverageZeroDirections || FETISeigDirections || FETISminmaxDirections || FETISeigDirections2
        OtherDirections = sum(Zk{1}(:,~KeptDirections),2);
        Zk{1} = Zk{1}(:,KeptDirections);
        if size(Zk{1},2) < Ns
            Zk{1}(:,end+1) = OtherDirections;
        end
    end
    NumDirectionsk = size(Zk{1},2);
    
        
        
    
    W2k{1}(:,:) = Ptrans*Zk{1}(:,:);
    
    % iterate
    while (true)
        
        Q2k{k}(:,:) = FItrans*W2k{k}(:,:);

        Deltak{k}(:,:) = Q2k{k}(:,:)'*W2k{k}(:,:);
        Gamma = Zk{k}(:,:)'*r(:,k);
        fac{k}(:,1) = (Deltak{k}(:,:))\Gamma;

        lambda(:,k+1) = lambda(:,k) + W2k{k}(:,:)*fac{k}(:,1);
        
        ComputeExactResidual = 0;
        if ComputeExactResidual
            r(:,k+1) = Ptrans'*(dtrans - FItrans*lambda(:,k+1));
        else
            r(:,k+1) = r(:,k) - Ptrans'*Q2k{k}(:,:)*fac{k}(:,1);
        end
        
        SummedDirection = zeros(Nlm,1);
        for s = 1:Ns
            Zk{k+1}(:,s) = FItransPreS{s}*r(:,k+1);
            %Zk{k+1}(:,s) = Zk{k+1}(:,s)./norm(Zk{k+1}(:,s));
            
            SummedDirection = SummedDirection + Zk{k+1}(:,s);
            
            if FETISeigDirections
                LocalNorm = norm(Bbs{s}*Ksbb{s}*Bbs{s}'*Zk{k+1}(:,s));
                GlobalNorm = norm(Bbs{s}*Bbs{s}'*FItransPre*Zk{k+1}(:,s));
                Amps(s) = GlobalNorm/LocalNorm;
            end
        end
        
        if FETISeigDirections2
            %SummedDirection = sum(Zk{k+1}(:,:),2);
            for s = 1:Ns
                GlobalNorm = norm(BBT{s}*SummedDirection);
                LocalNorm = norm(Zk{k+1}(:,s));
                Amps(s) = GlobalNorm/LocalNorm;
            end
        end
        
        if FETISminmaxDirections
            DirectionsAmps = sum(abs(Zk{k+1}(:,:)),1);
            [~, DirectionIndex1] = max(DirectionsAmps);
            [~, DirectionIndex2] = min(DirectionsAmps);
            KeptDirections(1:s) = false;
            KeptDirections(DirectionIndex1) = true;
            KeptDirections(DirectionIndex2) = true;
        end
        
        if FETISeigDirections || FETISeigDirections2
            [~, DirectionIndex] = max(Amps);
            KeptDirections(1:s) = false;
            KeptDirections(DirectionIndex) = true;
        end
    
        if FETISAverageZeroDirections
            % DirectionsAmps = sum(abs(Zk{k+1}(:,:)),1);
            DirectionsAmps = sqrt(sum((Zk{k+1}).^2,1));
            KeptDirections = DirectionsAmps > FETISDirectionTol;
            if isempty(find(KeptDirections,1))
                [maxDirectionAmp, DirectionIndex] = max(DirectionsAmps);
                KeptDirections(DirectionIndex) = true;
                %display(['Warning: No directions outside zero tolerance. Using direction with Amp = ' num2str(maxDirectionAmp)]);
            end
        end
        
        if FETISAverageZeroDirections || FETISeigDirections || FETISminmaxDirections || FETISeigDirections2
            OtherDirections = sum(Zk{k+1}(:,~KeptDirections),2);
        	Zk{k+1} = Zk{k+1}(:,KeptDirections);
            if size(Zk{k+1},2) < Ns
                Zk{k+1}(:,end+1) = OtherDirections;
            end
        end
        NumDirectionsk = NumDirectionsk + size(Zk{k+1},2);
        
        W2k{k+1}(:,:) = Ptrans*Zk{k+1}(:,:);

        for j = 1:k % full orthogonalization
            Theta = Q2k{j}(:,:)'*W2k{k+1}(:,:);
            W2k{k+1}(:,:) = W2k{k+1}(:,:) - W2k{j}(:,:)*((Deltak{j}(:,:))\Theta);
        end
        
        k = k + 1;

        %display(['norm(r(' num2str(k) ')) = ' num2str(norm(r(:,k)))]);
        %display(['norm(r(' num2str(k) '))/norm(r(' num2str(1) ')) = ' num2str(norm(r(:,k))/norm(r(:,1)))]);

        res(k) = sqrt(r(:,k)'*sum(Zk{k}(:,:),2));
        
        % stopping criterion
        if (res(k) < tol)
            break;
        elseif (NumDirectionsk >= Nlm)
            error('No Convergence');
%             elseif k > 10 && ...
%                 ((norm(r(:,k-10))-norm(r(:,k)))/norm(r(:,k-10)) < 0.1)
%                 break;
        end
    end
    lastlambda = lambda(:,k);
    iterations = iterations + k;
    NumDirections = NumDirections + NumDirectionsk;
    itdivisor = itdivisor + 1;
    Residual = res;

    uPPs = [];
    if Parallel
        for s = 1:Ns
            Ndof_end = sum(Ndof_s(1:s));
            Ndof_start = Ndof_end - Ndof_s(s) + 1;
            uPPs{s} = Ds{s}\(fs{s}(:,n+1) - Cs{s}*uAstP(Ndof_start:Ndof_end,n+1) - Ks{s}*uAst(Ndof_start:Ndof_end,n+1) - Bs{s}'*lambda(:,k));
        end
        for s = 1:Ns
            Ndof_end = sum(Ndof_s(1:s));
            Ndof_start = Ndof_end - Ndof_s(s) + 1;
            uPP(Ndof_start:Ndof_end,n+1) = uPPs{s};
        end
    else
        uPP(:,n+1) = D\(f(:,n+1) - C*uAstP(:,n+1) - K*uAst(:,n+1) - B'*lambda(:,k));
    end

    
%     if n == 1 || n == 5
%         figure(Figh.its);
%         set(Figh.its,'units','normalized','outerposition',[0 0 0.3 0.6]);
%         %Residual = sqrt(sum((w).^2,1));
%         semilogy(Residual);
%         set(gca,'YGrid','on')
%         set(gca,'YMinorGrid','off');
%         %set(gca,'YTick',[1e-10 1e-8 1e-6 1e-4 1e-2 1 1e2 1e4]);
%         %set(gca,'XTick',[0 20 40 60 80 100 120 140 160]);
%         xlim([0 k+1]);
%         ylim([min(Residual)/7 max(Residual)*7]);
%         set(gcf,'color','w');
%         if ~NoWrite
%         plot2svg([Experiment '\Ns' num2str(Nsx) '_It_n' num2str(n) '_h' num2str(Step) '_' CoarseGrid '.svg']);
%         end
%     end


    a(:,n+1) = aPre + ((1-GAalphaf).*uPP(:,n+1))./(1-GAalpham);
    uP(:,n+1) = uP(:,n) + Step*(1-GAgamma).*a(:,n) + Step*GAgamma.*a(:,n+1);
    u(:,n+1) = u(:,n) + Step.*uP(:,n) + Step^2*(0.5-GAbeta).*a(:,n) + Step^2*GAbeta.*a(:,n+1);

    %[maxval, index] = max(B*uPP(:,n));

    if DisplayTimeSteps
        display(['time: ' num2str(t(n+1)) '; iteration count: ' num2str(k) '; Residual = ' num2str(res(k))]);
    end
    %display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
end
display(['iterations: ' num2str(iterations) '; avg. directions per iteration: ' num2str(ceil(NumDirections/iterations))]);
averageIterations = iterations/itdivisor;
time = toc
AverageTime = time/n;

params.u = u;
params.Residual = Residual;
params.StaticIterations = k;

if ~NoWrite
    fid = fopen([Experiment '\log.txt'], 'a');
    fprintf(fid, ['\nNs' num2str(Nsx) ' Nely' num2str(Step) ' ' CoarseGrid ':' num2str(averageIterations) '\n']);
    fclose(fid);
end

    if exist([Experiment '\Iterations.mat'], 'file')
        load([Experiment '\Iterations.mat']);
    end
    IterationsPerStep(CaseNr,Subcase) = averageIterations;
    TimePerStep(CaseNr,Subcase) = AverageTime;
    DirectionsPerStep(CaseNr,Subcase) = ceil(NumDirections/iterations);
    CasesToSave{CaseNr} = description;
    save([Experiment '\Iterations.mat'],'IterationsPerStep','TimePerStep','DirectionsPerStep','CasesToSave');

end