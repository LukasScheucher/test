function [p] = FETIstaticsolverStest(p)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(p),struct2cell(p));


% initialize
k = 1;
vOne = ones(Ns,1);

lambdaN0 = Q*GI*((GI'*Q*GI)\e);
lambda(:,1) = zeros(Nlm,1);
r(:,1) = FETIpN((d - FI*lambdaN0(:,1)),1,p);

for s = 1:Ns
    Zk{1}(:,s) = FIPreS{s}*r(:,1);
end
%Zk{1} = sum(Zk{1}(:,s),2);

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
    SearchDirection = W2k{k}(:,:)*pinv(Deltak{k}(:,:))*Gamma;
    lambda(:,k+1) = lambda(:,k) + SearchDirection;
    r(:,k+1) = r(:,k) - FETIpN(Q2k{k}(:,:)*pinv(Deltak{k}(:,:))*Gamma,1,p);
    for s = 1:Ns
        Zk{k+1}(:,s) = FIPreS{s}*r(:,k+1);
    end
    %Zk{k+1} = sum(Zk{k+1}(:,s),2);
    
    
    W2k{k+1}(:,:) = FETIpN(Zk{k+1}(:,:),0,p);

    
    for j = 1:k % full orthogonalization
        Theta = Q2k{j}(:,:)'*W2k{k+1}(:,:);
        W2k{k+1}(:,:) = W2k{k+1}(:,:) - W2k{j}(:,:)*pinv(Deltak{j}(:,:))*Theta;
    end
    
    
    
    % % %
    LocalNorm = norm(Bbs{s}*Ksbb{s}*Bbs{s}'*SearchDirection);
    GlobalNorm = norm(Bbs{s}*Bbs{s}'*FIPre*SearchDirection);
    Amps(k) = GlobalNorm/LocalNorm;
    SearchDirections(:,k) = SearchDirection;
    % % %
    
    
    
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

p.plotRes = realRes;

p.Residual = res;

p.StaticIterations = k;
p.lambda = lambdafull;
p.alpha = alpha;

temp = num2str(sqrt(  r(:,k)' * sum(Zk{k}(:,:),2)  ));
display(['iteration count: ' num2str(k) '; realRes(k) = ' num2str(realRes(k))]);
display(num2str(res(1)));
display(num2str(res(2)));
display(num2str(res(3)));
%display(['direction count: ' num2str(NumDirectionsk)]);
%display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
%display(['iterations: ' num2str(p.StaticIterations)]);



display(num2str(Nlm));

[LowSortedAmps, LowIndices] = sort(Amps,'ascend');
[HighSortedAmps, HighIndices] = sort(Amps,'descend');
LowSearchDirections = SearchDirections(:,LowIndices);
HighSearchDirections = SearchDirections(:,HighIndices);


for HighLow = 1:3
    
for numDir = 1:size(Amps,2)
    %display('\n');
    
    if (HighLow == 1)
        % use HighDirections
        Directions = SearchDirections(:,HighIndices(1:numDir));
    elseif (HighLow == 2)
        % use LowDirections
        Directions = SearchDirections(:,LowIndices(1:numDir));
    elseif (HighLow == 3)
        Directions = SearchDirections(:,1:numDir);
    end

    
clear i j Zk W2k k vOne lambdaN0 lambda r lambdafull alpha res realRes Q2k Deltak Gamma Theta

% initialize
k = 1;
vOne = ones(Ns,1);

lambdaN0 = Q*GI*((GI'*Q*GI)\e);
lambda(:,1) = zeros(Nlm,1);
r(:,1) = FETIpN((d - FI*lambdaN0(:,1)),1,p);

%for s = 1:Ns
%    Zk{1}(:,s) = FIPreS{s}*r(:,1);
%end
Zk{1} = Directions;

W2k{1}(:,:) = Directions;%FETIpN(Zk{1}(:,:),0,p);

lambdafull(:,1) = lambdaN0 + lambda(:,1);
alpha(:,1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,1)));

res(k) = sqrt(  r(:,k)' * sum( Zk{k}(:,:) , 2 )  );
realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );

% iterate

%while (true)

    Q2k{k}(:,:) = FI*W2k{k}(:,:);
    Deltak{k}(:,:) = Q2k{k}(:,:)'*W2k{k}(:,:);
    Gamma = Zk{k}(:,:)'*r(:,k);
    lambda(:,k+1) = lambda(:,k) + W2k{k}(:,:)*pinv(Deltak{k}(:,:))*Gamma;
    r(:,k+1) = r(:,k) - FETIpN(Q2k{k}(:,:)*pinv(Deltak{k}(:,:))*Gamma,1,p);
    for s = 1:Ns
        Zk{k+1}(:,s) = FIPreS{s}*r(:,k+1);
    end
    
    
    W2k{k+1}(:,:) = FETIpN(Zk{k+1}(:,:),0,p);

    
    %for j = 1:k % full orthogonalization
    %    Theta = Q2k{j}(:,:)'*W2k{k+1}(:,:);
    %    W2k{k+1}(:,:) = W2k{k+1}(:,:) - W2k{j}(:,:)*pinv(Deltak{j}(:,:))*Theta;
    %end
    
    
    lambdafull(:,k+1) = lambdaN0 + lambda(:,k+1);
    alpha(:,k+1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,k+1)));
    
    k = k + 1;
    
    res(k) = sqrt(  r(:,k)' * sum(Zk{k}(:,:),2)  );
    realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );
    
%end

p.plotRes = realRes;

p.Residual = res;

p.StaticIterations = k;
p.lambda = lambdafull;
p.alpha = alpha;

if HighLow == 1
    %display(['High Directions (' num2str(numDir) '):']);
    ResidualsHigh(numDir) = realRes(k);
elseif HighLow == 2
    %display(['Low Directions (' num2str(numDir) '):']);
    ResidualsLow(numDir) = realRes(k);
elseif HighLow == 3
    display(['Unsorted Directions (' num2str(numDir) '):']);
    ResidualsUnsorted(numDir) = realRes(k);
end
display(['iteration count: ' num2str(k) '; realRes = ' num2str(realRes(k))]);

end
end

figure(1);
plot([1:numDir],ResidualsHigh);
title('use directions with large factor')

figure(2);
plot([1:numDir],ResidualsLow);
title('use directions with small factor');

figure(3);
plot([1:numDir],ResidualsUnsorted);
title('use directions unsorted');

end