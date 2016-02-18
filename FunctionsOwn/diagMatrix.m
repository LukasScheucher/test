function matrix = diagMatrix(matrices,n,m)

matrix = zeros(n,m); % sparse(n,m);
ni = 1;
mi = 1;
for i = 1:length(matrices)
    niend = ni+size(matrices{i},1) - 1;
    miend = mi+size(matrices{i},2) - 1;
    matrix(ni:niend,mi:miend) ...
        = matrices{i};
    ni = niend + 1;
    mi = miend + 1;
end

end