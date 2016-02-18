function [p] = FETIstaticsolverExp1(p)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(p),struct2cell(p));

% initialize
k = 1;

lambdaN0 = Q*GI*((GI'*Q*GI)\e);
lambdaC0 = GI0*((GI0'*FI*GI0)\(GI0'*(d-FI*lambdaN0)));
lambda(:,1) = zeros(Nlm,1);
r(:,1) = PN'*PC'*(d - FI*lambdaN0(:,1));
z(:,1) = FIPre*r(:,1);
w(:,1) = PN*z(:,1);

lambdafull(:,1) = lambdaN0 + lambdaC0 + PC*lambda(:,1);
alpha(:,1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,1)));

res(k) = sqrt(r(:,k)'*z(:,k));
realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );

% iterate
tic
while (true)

    q(:,k) = PC'*FI*w(:,k);
    delta = q(:,k)'*w(:,k);
    gamma = r(:,k)'*z(:,k);
    lambda(:,k+1) = lambda(:,k) + (gamma/delta).*w(:,k);
    r(:,k+1) = r(:,k) - (gamma/delta).*PN'*q(:,k);
    z(:,k+1) = FIPre*r(:,k+1);
    w(:,k+1) = PN*z(:,k+1);
    
    for j = k%1:k % full orthogonalization
        w(:,k+1) = w(:,k+1) - ( (q(:,j)'*w(:,k+1) )/( q(:,j)'*w(:,j) )).*w(:,j);
    end
    
    lambdafull(:,k+1) = lambdaN0 + lambdaC0 + PC*lambda(:,k+1);
    lambdafullUnprojected(:,k+1) = lambdaN0 + lambdaC0 + lambda(:,k+1);
    alpha(:,k+1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,k+1)));
    
    k = k + 1;

    % stopping criterion
    %display([num2str(norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) )) ' = norm( d - FI*lambdafull(' num2str(k) ') - GI*alpha(' num2str(k) ') )']);
    %display([num2str(norm( PN'*PC'*(d - FI*lambdafullUnprojected(:,k)) )) ' = PN''*PC''*(d - FI*lambdafullUnprojected(' num2str(k) '))']);
    %display([num2str(sqrt(r(:,k)'*z(:,k))) ' = sqrt(r(' num2str(k) ')''*z(' num2str(k) '))']);
    %display([num2str(norm(w(:,k))) ' = norm(w(' num2str(k) '))']);
    %display('\n');
    
    res(k) = sqrt(r(:,k)'*z(:,k));
    realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );

    %display(['norm(r(' num2str(k) '))/norm(r(' num2str(1) ')) = ' num2str(norm(r(:,k))/norm(r(:,1)))]);
    if (sqrt(r(:,k)'*z(:,k)) < 1e-7)
        break;
    elseif k == Nlm
        break;
    end
    
end
toc

p.Residual = realRes;

p.StaticIterations = k;
p.lambda = lambdafull;
p.alpha = alpha;

display(['iteration count: ' num2str(k) '; norm(w(' num2str(k) ')) = ' num2str(norm(w(:,k)))]);
%display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
display(['iterations: ' num2str(p.StaticIterations)]);
end