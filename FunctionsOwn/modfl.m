function rest=modfl(x,y,tol)
    % Calculates the modulus of x and y for floatingpoint numbers
    %rest=mod(round(x/tol),round(y/tol))*tol;
    y_add=0;
    while abs(x-y_add)>tol && y_add-y<=x
        y_add=y_add+y;
    end
    if abs(x-y_add)<=tol
        rest=0;
    else
        rest=mod(x,y);
    end
end