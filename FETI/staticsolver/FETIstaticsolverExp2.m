function [p] = FETIstaticsolverExp2(p)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(p),struct2cell(p));

tic

% initialize
k = 1;

if Coarse == 3
    lambdaN0 = zeros(Nlm,1);
else
    lambdaN0 = Q*GI*((GI'*Q*GI)\e);
end
if Coarse >= 2
    lambdaC0 = GI0*((GI0'*FI*GI0)\(GI0'*(d-FI*lambdaN0)));
else
    lambdaC0 = zeros(Nlm,1);
end
lambda(:,1) = zeros(Nlm,1);
r(:,1) = FETIpN(FETIpC((d - FI*lambdaN0(:,1)),1,p),1,p);
z(:,1) = FIPre*r(:,1);
w(:,1) = FETIpN(z(:,1),0,p);

lambdafull(:,1) = lambdaN0 + lambdaC0 + FETIpC(lambda(:,1),0,p);
alpha(:,1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,1)));

res(k) = sqrt(r(:,k)'*z(:,k));
realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );

% iterate

while (true)

    q(:,k) = FETIpC(FI*w(:,k),1,p);
    delta = q(:,k)'*w(:,k);
    gamma = r(:,k)'*z(:,k);
    lambda(:,k+1) = lambda(:,k) + (gamma/delta).*w(:,k);
    r(:,k+1) = r(:,k) - (gamma/delta).*FETIpN(q(:,k),1,p);
    z(:,k+1) = FIPre*r(:,k+1);
    w(:,k+1) = FETIpN(z(:,k+1),0,p);
    
    for j = k%1:k % full orthogonalization
        w(:,k+1) = w(:,k+1) - ( (q(:,j)'*w(:,k+1) )/( q(:,j)'*w(:,j) )).*w(:,j);
    end
    
    lambdafull(:,k+1) = lambdaN0 + lambdaC0 + FETIpC(lambda(:,k+1),0,p);
    lambdafullUnprojected(:,k+1) = lambdaN0 + lambdaC0 + lambda(:,k+1);
    if Coarse == 3
        alpha(:,k+1) = zeros(size(GI,2),1);
    else
        alpha(:,k+1) = (GI'*GI)\(GI'*(d-FI*lambdafull(:,k+1)));
    end
    
    k = k + 1;
    
    % stopping criterion
    %display([num2str(norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) )) ' = norm( d - FI*lambdafull(' num2str(k) ') - GI*alpha(' num2str(k) ') )']);
    %display([num2str(norm( FETIpN(FETIpC( (d - FI*lambdafullUnprojected(:,k)) ,1,p),1,p) )) ' = PN''*PC''*(d - FI*lambdafullUnprojected(' num2str(k) '))']);
    %display([num2str(sqrt(r(:,k)'*z(:,k))) ' = sqrt(r(' num2str(k) ')''*z(' num2str(k) '))']);
    %display([num2str(norm(w(:,k))) ' = norm(w(' num2str(k) '))']);
    %display('\n');
    
    res(k) = sqrt(r(:,k)'*z(:,k));
    realRes(k) = norm( d - FI*lambdafull(:,k) - GI*alpha(:,k) );
    
    %display(['norm(r(' num2str(k) '))/norm(r(' num2str(1) ')) = ' num2str(norm(r(:,k))/norm(r(:,1)))]);
    if (res(k) < tol)
        break;
    elseif k == max_iteration*Nlm
        disp('Max. number of iterations reached!')
        break;
    end
    
end

toc

p.Residual = res;

p.StaticIterations = k;
p.lambda = lambdafull;
p.alpha = alpha;


%display(['iteration count: ' num2str(k) '; norm(w(' num2str(k) ')) = ' num2str(norm(w(:,k)))]);
temp = num2str(sqrt(r(:,k)'*z(:,k)));
display(['iteration count: ' num2str(k) '; sqrt(r(:,' num2str(k) ')''*z(:,' num2str(k) ')) = ' temp]);
%display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
%display(['iterations: ' num2str(p.StaticIterations)]);
end