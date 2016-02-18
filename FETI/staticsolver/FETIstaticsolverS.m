function [p] = FETIstaticsolverS(p)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(p),struct2cell(p));

tic

% tune external load
if FETISTuneLoad
    AddForceMax = 0;
    AddForceMin = 1;
    AddForce = (AddForceMax*rand(Nlm,1) + AddForceMin*ones(Nlm,1));
    for s = 1:Ns
        fs{s}(:,1) = fs{s}(:,1) + Bs{s}'*AddForce;
    end
    f(:,1) = f(:,1) + B'*AddForce;
end


% initialize
k = 1;
vOne = ones(Ns,1);

lambdaN0 = Q*GI*((GI'*Q*GI)\e);
lambda(:,1) = zeros(Nlm,1);
r(:,1) = FETIpN((d - FI*lambdaN0(:,1)),1,p);

for s = 1:Ns
    Zk{1}(:,s) = FIPreS{s}*r(:,1);
    BBT{s} = Bs{s}*Bs{s}';
end

if FETISeigDirections1
    for s = 1:Ns
        LocalNorm = norm(Bbs{s}*Ksbb{s}*Bbs{s}'*Zk{1}(:,s));
        GlobalNorm = norm(Bbs{s}*Bbs{s}'*FIPre*Zk{1}(:,s));
        Amps(s) = GlobalNorm/LocalNorm;
    end
end

if FETISeigDirections2
    SummedDirection = sum(Zk{1}(:,:),2);
    for s = 1:Ns
        GlobalNorm = norm(BBT{s}*SummedDirection);
        LocalNorm = norm(Zk{1}(:,s));
        Amps(s) = GlobalNorm/LocalNorm;
    end
end

if FETISeigDirections1 || FETISeigDirections2
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

if FETISAverageZeroDirections || FETISeigDirections1 || FETISminmaxDirections || FETISeigDirections2
    OtherDirections = sum(Zk{1}(:,~KeptDirections),2);
    Zk{1} = Zk{1}(:,KeptDirections);
    if size(Zk{1},2) < Ns
        Zk{1}(:,end+1) = OtherDirections;
    end
end
NumDirectionsk = size(Zk{1},2);




W2k{1}(:,:) = FETIpN(Zk{1}(:,:),0,p);

lambdafull(:,1) = lambdaN0 + lambda(:,1);
alpha(:,1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,1)));

res(k) = sqrt(  r(:,k)' * sum( Zk{k}(:,:) , 2 )  );
realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );

% iterate

while (true)

    Q2k{k}(:,:) = FI*W2k{k}(:,:);
    Deltak{k}(:,:) = Q2k{k}(:,:)'*W2k{k}(:,:);
    Gamma = Zk{k}(:,:)'*r(:,k);
    lambda(:,k+1) = lambda(:,k) + W2k{k}(:,:)*pinv(Deltak{k}(:,:))*Gamma;
    r(:,k+1) = r(:,k) - FETIpN(Q2k{k}(:,:)*pinv(Deltak{k}(:,:))*Gamma,1,p);
    for s = 1:Ns
        Zk{k+1}(:,s) = FIPreS{s}*r(:,k+1);
    end
    
    
    
    
    
    SummedDirection = zeros(Nlm,1);
    for s = 1:Ns
        SummedDirection = SummedDirection + Zk{k+1}(:,s);
    end
    
    if FETISeigDirections1
        for s = 1:Ns
            LocalNorm = norm(Bbs{s}*Ksbb{s}*Bbs{s}'*Zk{k+1}(:,s));
            GlobalNorm = norm(Bbs{s}*Bbs{s}'*FIPre*Zk{k+1}(:,s));
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

    if FETISeigDirections1 || FETISeigDirections2
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

    if FETISAverageZeroDirections || FETISeigDirections1 || FETISminmaxDirections || FETISeigDirections2
        OtherDirections = sum(Zk{k+1}(:,~KeptDirections),2);
        Zk{k+1} = Zk{k+1}(:,KeptDirections);
        if size(Zk{k+1},2) < Ns
            Zk{k+1}(:,end+1) = OtherDirections;
        end
    end
    NumDirectionsk = NumDirectionsk + size(Zk{k+1},2);
    
    
    
    
    
    
    W2k{k+1}(:,:) = FETIpN(Zk{k+1}(:,:),0,p);
    
    for j = 1:k % full orthogonalization
        Theta = Q2k{j}(:,:)'*W2k{k+1}(:,:);
        W2k{k+1}(:,:) = W2k{k+1}(:,:) - W2k{j}(:,:)*pinv(Deltak{j}(:,:))*Theta;
    end
    
    lambdafull(:,k+1) = lambdaN0 + lambda(:,k+1);
    alpha(:,k+1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,k+1)));
    
    k = k + 1;
    
    % stopping criterion
    %display([num2str(norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) )) ' = norm( d - FI*lambdafull(' num2str(k) ') - GI*alpha(' num2str(k) ') )']);
    %display([num2str(norm( FETIpN(FETIpC( (d - FI*lambdafullUnprojected(:,k)) ,1,p),1,p) )) ' = PN''*PC''*(d - FI*lambdafullUnprojected(' num2str(k) '))']);
    %display([num2str(sqrt(r(:,k)'*z(:,k))) ' = sqrt(r(' num2str(k) ')''*z(' num2str(k) '))']);
    %display([num2str(norm(w(:,k))) ' = norm(w(' num2str(k) '))']);
    %display('\n');
    
    res(k) = sqrt(  r(:,k)' * sum(Zk{k}(:,:),2)  );
    realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );
    
    %display(['norm(r(' num2str(k) '))/norm(r(' num2str(1) ')) = ' num2str(norm(r(:,k))/norm(r(:,1)))]);
    if (res(k) < tol)
        break;
    elseif k == Nlm
        break;
    end
    
end
toc

p.plotRes = realRes;

p.Residual = res;

p.StaticIterations = k;
p.lambda = lambdafull;
p.alpha = alpha;

temp = num2str(sqrt(  r(:,k)' * sum(Zk{k}(:,:),2)  ));
display(['iteration count: ' num2str(k) '; sqrt(  r(:,' num2str(k) ')'' * sum(Zk{' num2str(k) '}(:,:),2)  ) = ' temp]);
display(['direction count: ' num2str(NumDirectionsk)]);
%display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
%display(['iterations: ' num2str(p.StaticIterations)]);
end