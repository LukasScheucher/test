function N = shapefunc(i,derivative,xi)
% Linear shape functions
tol=0.0000001;
switch i
    case 1
        if derivative==1
            N=-0.5;
        else
            N=0.5*(1-xi);
            if abs(xi-1)<=tol
                N=0;
            else
                if abs(xi+1)<=tol
                    N=1;
                end
            end
        end
    case 2
        if derivative==1
            N=0.5;
        else        
            N=0.5*(1+xi);
            if abs(xi+1)<=tol
                N=0;
            else
                if abs(xi-1)<=tol
                    N=1;
                end
            end
        end
end

end