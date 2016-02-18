function params = FETIdynamicsolver(params)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));

%% Teste GenEo - Tenhumberg
TestGenEo = 0;
if (TestGenEo && Coarse > 1)
    [GI] = TestGI(GI, Odd, Nsx*Nsy, III);
    
    % Projektion entsprechend anpassen
    CoarseOp = GI'*FItrans*GI;
    Ptrans = eye(Nlm) - GI*(CoarseOp\(GI'*FItrans));
end

%%




% init values
params.ti_ti = []; % calculation-time for timestep i
params.it_ti = []; % iteration for timestep i
result_L = []; % für Vergleich mit FullSolver

u = zeros(Ndof,length(t)); % u: displacements
uP = u; % uP: derivative of u / velocities
uPP = u; % uPP: derivative of uP / accelerations
uAst = u; % uAst: Predictor for u in generalized alpha integration


% calculate initial lambda
Temp1 = (B*(M\B'));
Temp2 = (M\(f(:,1)-C*uP(:,1)-K*u(:,1)));
lambda = pinv(Temp1)*(B*Temp2);


% caculate initial accelerations uPP using the inital lambda
% here, parallel means that the calculations are done for each substructure
% separately, it does not mean parallelization in time
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
    % compute all substructures in one big matrix equation
    % this is NOT recommended, just academic
    uPP(:,1) = M\(f(:,1)-C*uP(:,1)-K*u(:,1)-B'*lambda);
end


iterations = 0;
itdivisor = 0;
a(:,1) = uPP(:,1);
lastlambda = [];
tic
% timeintegration starts here
for n = 1:(length(t)-1)

    % predictors for generalized alpha integration
    aPre = (GAalphaf.*uPP(:,n)-GAalpham.*a(:,n))./(1-GAalpham);
    uAst(:,n+1) = u(:,n) + Step.*uP(:,n) + Step^2*(0.5-GAbeta).*a(:,n) + Step^2*GAbeta.*aPre;
    uAstP(:,n+1) = uP(:,n) + Step*(1-GAgamma).*a(:,n) + Step*GAgamma.*aPre;

    % computation of d in the transient case
    % meaning of "Parallel" see above
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

    % initialize PCPG
    tic
    k = 1;
    if Coarse
        % use a coarse grid, coarse space defined in GI
        % select some admisible lambda (a lambda that meets the coarse
        % space requirements, that means it solves the problem within the
        % coarse space defined by the columns of GI)
        lambda(1:Nlm,k) = GI*(CoarseOp\(GI'*dtrans));
    else
        % do not use any coarse grid, projector is set to eye
        % and lambda can start with zero
        Ptrans = eye(Nlm);
        lambda(:,k) = zeros(Nlm,1);
    end
    
    r(:,1) = Ptrans'*(dtrans - FItrans*lambda(:,k)); % projection
    z(:,1) = FItransPre*r(:,1); % preconditioning
    w(:,1) = Ptrans*z(:,1); % reprojection
    % w is the new search direction
    
    % w(:,1) = z(:,1); % no reprojection: doesnt converge (the
    % preconditioner does not preserve the projection)

    
    % iterate
    while (true)
        
        % calculate step size alpha for the new search direction w
        q(:,k) = FItrans*w(:,k);
        alpha = (z(:,k)'*r(:,k)) / (q(:,k)'*w(:,k));
        
        % update lambda
        lambda(:,k+1) = lambda(:,k) + alpha.*w(:,k);
        % update residual (already projected)
        r(:,k+1) = r(:,k) - Ptrans'*(alpha.*q(:,k));
        
        % compute new search direction by preconditioning and reprojecting
        z(:,k+1) = FItransPre*r(:,k+1);
        w(:,k+1) = Ptrans*z(:,k+1);
        
        
        % w(:,k+1) = z(:,k+1); % no reprojection: doesnt converge (the
        % preconditioner does not preserve the projection)
        
        
        % orthogonalize the new search direction against all previous ones
        % in CG theory only orthogonalization against the last search
        % direction is necessary, practically, full reorthogonalization is
        % required due to numerical inaccuracies
        for j = 1:k
            xi = (q(:,j)'*w(:,k+1)) / (q(:,j)'*w(:,j));
            w(:,k+1) = w(:,k+1) - xi.*w(:,j);
        end

        
        % next iteration
        k = k + 1;
        
        % stopping criterion
        res(k) = sqrt(r(:,k)'*(z(:,k)));
        if isfield(params,'DisplayStepResidual') && params.DisplayStepResidual
            display(['Timstep ', num2str(n) ', Residual ', num2str(res(k))]);
        end
        if (res(k) < tol)
            break;
        elseif k == Nlm
            break;
%             elseif k > 10 && ...
%                 ((norm(r(:,k-10))-norm(r(:,k)))/norm(r(:,k-10)) < 0.1)
%                 break;
        end
    end
    % PCPG iteration end here
    
    % help variables
    params.ti_ti = [params.ti_ti, toc];
    params.it_ti = [params.it_ti, k];
    
    lastlambda = lambda(:,k);    
    result_L = [result_L, lastlambda]; 
    % keep track of number of iterations over all timesteps
    iterations = iterations + k;
    itdivisor = itdivisor + 1;
    Residual = res(k);

    
    % PCPG has converged, the new lambda is determined
    % compute the corresponding accelerations
    % meaning of "Parallel": see above
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
%         Residual = sqrt(sum((w).^2,1));
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



    % generalized alpha:
    % calculate velocities uP and displacements u using generalized alpha
    % integration:
    a(:,n+1) = aPre + ((1-GAalphaf).*uPP(:,n+1))./(1-GAalpham);
    uP(:,n+1) = uP(:,n) + Step*(1-GAgamma).*a(:,n) + Step*GAgamma.*a(:,n+1);
    u(:,n+1) = u(:,n) + Step.*uP(:,n) + Step^2*(0.5-GAbeta).*a(:,n) + Step^2*GAbeta.*a(:,n+1);

    
    %[maxval, index] = max(B*uPP(:,n));
    if DisplayTimeSteps
        display(['time: ' num2str(t(n+1)) '; iteration count: ' num2str(k) '; norm(r(' num2str(k) ')) = ' num2str(norm(r(:,1)))]);
    end
    %display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
end
% time integration ends here


display(['iterations: ' num2str(iterations)]);
display(['DoF = ' num2str(size(lambda,1)) ]);
averageIterations = iterations/itdivisor;
time = toc;
AverageTime = time/n;
display(['CPU time: ' num2str(time) 's']);
% load('Odd_cS');
% %str = ['Odd' num2str(Odd) 'cS' num2str(cS)];
% Odd_cS{Odd,cS} = params.it_ti;
% save('Odd_cS','Odd_cS')


% load('Odd_cS');
% str = ['Odd' num2str(Odd) 'cS' num2str(cS)];
% Odd_cS{Odd,cS} = params.it_ti;
% save('Odd_cS','Odd_cS')

params.u = u;
params.Residual = Residual;
params.StaticIterations = k;


% write some logfiles
if ~NoWrite
    fid = fopen([Experiment '\log.txt'], 'a');
    fprintf(fid, ['\nNs' num2str(Nsx) ' Nely' num2str(Step) ' ' CoarseGrid ':' num2str(averageIterations) '\n']);
    fclose(fid);
end

if exist([Experiment '\Iterations.mat'], 'file')
    load([Experiment '\Iterations.mat']);
end
IterationsPerStep(Case,Subcase) = averageIterations;
TimePerStep(Case,Subcase) = AverageTime;
DirectionsPerStep(Case,Subcase) = 1;
CasesToSave{Case} = description;
save([Experiment '\Iterations.mat'],'IterationsPerStep','TimePerStep','DirectionsPerStep','CasesToSave');


end