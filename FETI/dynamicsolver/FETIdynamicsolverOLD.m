function params = FETIdynamicsolverOLD(params)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));

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
itdivisor = 0;
a(:,1) = uPP(:,1);
lastlambda = [];
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
    clear k r lambda w wd z zd y xi p v
    k = 1;

    if Coarse
        lambda(1:Nlm,k) = GI*(CoarseOp\(GI'*dtrans));
        r(:,1) = dtrans - FItrans*lambda(:,k);
    else
        Ptrans = eye(Nlm);
        lambda(:,k) = zeros(Nlm,1);
        r(:,1) = dtrans - FItrans*lambda(:,k);
    end

    % iterate
    while (true)
        k = k + 1;

        % preconditioning and projection
        if Coarse
            w(:,k-1) = Ptrans'*r(:,1);
            %norm(w(:,k-1)-r(:,k-1))
            z(:,k-1) = FItransPre*w(:,k-1);
            y(:,k-1) = Ptrans*z(:,k-1);
            %norm(y(:,k-1)-z(:,k-1))
        else
            w(:,k-1) = r(:,1);
            z(:,k-1) = FItransPre*w(:,k-1);
            y(:,k-1) = z(:,k-1);
        end

        % conjugate gradient
        if (k > 2)
            xi(k) = (y(:,k-1)'*w(:,k-1)) / (y(:,k-2)'*w(:,k-2));
            p(:,k) = y(:,k-1) + xi(k).*p(:,k-1);
        else
            p(:,k) = y(:,k-1);
        end
        temp = FItrans*p(:,k);
        v(k) = (y(:,k-1)'*w(:,k-1)) / (p(:,k)'*temp);
        lambda(:,k) = lambda(:,k-1) + v(k)*p(:,k);
        r = r - v(k).*temp;
        %display(['norm(r(' num2str(k) ')) = ' num2str(norm(r(:,k)))]);
        %display(['norm(r(' num2str(k) '))/norm(r(' num2str(1) ')) = ' num2str(norm(r(:,k))/norm(r(:,1)))]);

        if Coarse
            w(:,k) = Ptrans'*r(:,1);
        else
            w(:,k) = r(:,1);
        end

        % stopping criterion
        if (norm(w(:,k)) < tol)
            break;
        elseif k == Nlm
            break;
%             elseif k > 10 && ...
%                 ((norm(r(:,k-10))-norm(r(:,k)))/norm(r(:,k-10)) < 0.1)
%                 break;
        end
    end
    lastlambda = lambda(:,k);
    iterations = iterations + k;
    itdivisor = itdivisor + 1;
    Residual = sqrt(sum((w).^2,1));

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


    a(:,n+1) = aPre + ((1-GAalphaf).*uPP(:,n+1))./(1-GAalpham);
    uP(:,n+1) = uP(:,n) + Step*(1-GAgamma).*a(:,n) + Step*GAgamma.*a(:,n+1);
    u(:,n+1) = u(:,n) + Step.*uP(:,n) + Step^2*(0.5-GAbeta).*a(:,n) + Step^2*GAbeta.*a(:,n+1);

    %[maxval, index] = max(B*uPP(:,n));

    display(['time: ' num2str(t(n+1)) '; iteration count: ' num2str(k) '; norm(r(' num2str(k) ')) = ' num2str(norm(r(:,1)))]);
    %display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
end
display(['iterations: ' num2str(iterations)]);
averageIterations = iterations/itdivisor;


params.u = u;
params.Residual = Residual;
params.StaticIterations = k;


fid = fopen([Experiment '\log.txt'], 'a');
fprintf(fid, ['\nNs' num2str(Nsx) ' Nely' num2str(Step) ' ' CoarseGrid ':' num2str(averageIterations) '\n']);
fclose(fid);

if exist([Experiment '\Iterations.mat'], 'file')
    load([Experiment '\Iterations.mat']);
end
Fieldname = ['Ns' num2str(Nsx)];
Timesteps{timestepNo} = num2str(Step);
IterationsPerTimestep.(Fieldname)(Coarse+1,timestepNo) = averageIterations;
if ~NoWrite
save([Experiment '\Iterations.mat'],'IterationsPerTimestep','Timesteps');
end

end