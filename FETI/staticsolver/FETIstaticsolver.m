function [params] = FETIstaticsolver(params)

cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));

% initialize
clear k r lambda w wd z zd y xi p v
k = 1;

% Np = pinv(N); no convergence with this!

if Coarse == 2
    startLMs = N\[GI0'*d; GI'*Q*d; e];
else
    startLMs = N\[GI'*Q*d; e];
end
if Ngm>0
    mu(:,k) = startLMs(1:Ngm);
else
    mu(1,k) = 0;
end
gamma(:,k) = startLMs(Ngm+1:Ngm+Nrbm);
alpha(:,k) = startLMs(Ngm+Nrbm+1:end);
lambda(:,k) = Q*GI*gamma(:,k) + GI0*mu(:,k);
w(:,k) = d - FI*lambda(:,k) - GI*alpha(:,k);
p(:,k) = FIPre*w(:,k);

% iterate
while (true)

    if k > 1
        for iOrth = 1:(k-1)
            p(:,k) = p(:,k) - ((Dw(:,iOrth)'*p(:,k))/(Dw(:,iOrth)'*Dlambda(:,iOrth))).*Dlambda(:,iOrth);
        end
    end

    if Coarse == 2
        LMs = N\[-GI0'*FI*p(:,k); -GI'*Q*FI*p(:,k); -GI'*p(:,k)];
    else
        LMs = N\[-GI'*Q*FI*p(:,k); -GI'*p(:,k)];
    end

    if Ngm>0
        mu(:,k) = LMs(1:Ngm);
    else
        mu(1,k) = 0;
    end
    gamma(:,k) = LMs(Ngm+1:Ngm+Nrbm);
    Dalpha(:,k) = LMs(Ngm+Nrbm+1:end);

    Dlambda(:,k) = p(:,k) + Q*GI*gamma(:,k) + GI0*mu(:,k);
    Dw(:,k) = -FI*Dlambda(:,k) - GI*Dalpha(:,k);

    eta(:,k) = -(Dlambda(:,k)'*w(:,k))/(Dlambda(:,k)'*Dw(:,k));
    lambda(:,k+1) = lambda(:,k) + eta(:,k)*Dlambda(:,k);
    w(:,k+1) = w(:,k) + eta(:,k)*Dw(:,k);
    alpha(:,k+1) = alpha(:,k) + eta(:,k)*Dalpha(:,k);
    
    %w(:,k+1) = d - FI*lambda(:,k+1) - GI*alpha(:,k+1);
    
    p(:,k+1) = FIPre*w(:,k);

    display(['norm(w(' num2str(k) ')) = ' num2str(norm(w(:,k)))]);
    %display(['norm(r(' num2str(k) '))/norm(r(' num2str(1) ')) = ' num2str(norm(r(:,k))/norm(r(:,1)))]);


    % stopping criterion
    if (norm(w(:,k+1)) < tol)
        break;
    elseif k+1 == Nlm
        break;
%             elseif k > 10 && ...
%                 ((norm(r(:,k-10))-norm(r(:,k)))/norm(r(:,k-10)) < 0.1)
%                 break;
    end

    k = k + 1;
end
params.StaticIterations = k+1;
params.alpha = alpha;
params.lambda = lambda;


display(['iteration count: ' num2str(k) '; norm(w(' num2str(k) ')) = ' num2str(norm(w(:,k)))]);
%display(['maxval = ' num2str(maxval) '; index = ' num2str(index)]);
display(['iterations: ' num2str(params.StaticIterations)]);
end