function [p] = FETIstaticsolverS2(p)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(p),struct2cell(p));

% initialize
k = 1;
vOne = ones(Ns,1);

lambdaN0 = Q*GI*((GI'*Q*GI)\e);
lambda(:,1) = zeros(Nlm,1);
r(:,1) = PN'*(d - FI*lambdaN0(:,1));
for s = 1:Ns
    Z(:,s,1) = FIPreS{s}*r(:,1);
end
W2(:,:,1) = PN*(Z(:,:,1));

lambdafull(:,1) = lambdaN0 + lambda(:,1);
alpha(:,1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,1)));

res(k) = sqrt(r(:,k)'*(Z(:,:,k)*vOne));
realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );

% iterate
tic
while (true)

    Q2(:,:,k) = FI*W2(:,:,k);
    Delta(:,:,k) = Q2(:,:,k)'*W2(:,:,k);
    Gamma = Z(:,:,k)'*r(:,k);
    lambda(:,k+1) = lambda(:,k) + W2(:,:,k)*pinv(Delta(:,:,k))*Gamma;
    r(:,k+1) = r(:,k) - PN'*(Q2(:,:,k)*pinv(Delta(:,:,k))*Gamma);
    for s = 1:Ns
        Z(:,s,k+1) = FIPreS{s}*r(:,k+1);
    end
    W2(:,:,k+1) = PN*(Z(:,:,k+1));
    
    for j = 1:k % full orthogonalization
        Theta = Q2(:,:,j)'*W2(:,:,k+1);
        W2(:,:,k+1) = W2(:,:,k+1) - W2(:,:,j)*pinv(Delta(:,:,j))*Theta;
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
    
    res(k) = sqrt(r(:,k)'*(Z(:,:,k)*vOne));
    realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );
    
    %display(['norm(r(' num2str(k) '))/norm(r(' num2str(1) ')) = ' num2str(norm(r(:,k))/norm(r(:,1)))]);
    if (res(k) < 1e-9)
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


display(['iteration count: ' num2str(k) '; norm(w(' num2str(k) ')) = ' num2str(norm(W2(:,k)))]);
%display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
display(['iterations: ' num2str(p.StaticIterations)]);
end